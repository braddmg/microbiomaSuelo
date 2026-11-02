# Lee la lista de taxones de interés. Acepta CSV con `,` o `;`, y UTF-8 o
# Latin-1 (lo habitual al guardar desde Excel en Windows).
leer_taxones <- function(archivo) {
  if (!file.exists(archivo)) stop("No existe la lista de taxones: ", archivo)

  lineas <- readLines(archivo, warn = FALSE)
  if (all(validUTF8(lineas))) {
    Encoding(lineas) <- "UTF-8"
  } else {
    lineas <- iconv(lineas, "latin1", "UTF-8")
  }
  lineas <- sub("^﻿", "", lineas)
  sep <- if (grepl(";", lineas[1], fixed = TRUE)) ";" else ","
  x <- utils::read.csv(text = lineas, sep = sep, stringsAsFactors = FALSE,
                       check.names = FALSE, colClasses = "character")

  faltan <- setdiff(c("Grupo", "Genero", "Especie", "Clasificacion"), names(x))
  if (length(faltan) > 0) {
    stop("Faltan columnas en la lista de taxones: ", paste(faltan, collapse = ", "))
  }
  if (!"Funcion" %in% names(x)) x$Funcion <- ""

  x$Genero <- trimws(gsub("*", "", x$Genero, fixed = TRUE))
  x$Especie <- trimws(gsub("*", "", x$Especie, fixed = TRUE))
  x <- x[nzchar(x$Genero), , drop = FALSE]
  # Especie admite "R. solani", "Rhizoctonia solani" o "solani"; vacía = solo género.
  x$Epiteto <- tolower(sub(".*\\s", "", x$Especie))

  clase <- tolower(trimws(x$Clasificacion))
  x$Clasificacion <- ifelse(startsWith(clase, "pat"), "Patógeno",
                            ifelse(startsWith(clase, "ben"), "Beneficioso", NA))
  if (anyNA(x$Clasificacion)) {
    stop("Clasificacion debe ser 'Patógeno' o 'Beneficioso'. Revise: ",
         paste(x$Genero[is.na(x$Clasificacion)], collapse = ", "))
  }
  x
}

# Coincide por partes del nombre: SILVA agrupa géneros en nombres compuestos
# como "Burkholderia-Caballeronia-Paraburkholderia".
coincide_genero <- function(generos, genero) {
  vapply(strsplit(tolower(generos), "[-/_ ]+"), function(p) tolower(genero) %in% p, logical(1))
}

#' Abundancia de taxones de interés (potenciales patógenos y beneficiosos)
#'
#' Busca primero cada especie de la lista y luego su género. Un género es
#' potencial patógeno o beneficioso si todas sus especies en la lista lo son,
#' y `No concluyente` si la lista incluye especies de ambas clases.
#'
#' @param ps Objeto phyloseq.
#' @param tipo_datos `"16S"` (usa filas `Bacteria`) o `"ITS"` (filas `Hongo`).
#' @param variable_analisis Columna de agrupación.
#' @param archivo_taxones CSV del usuario (`,` o `;`). Columnas: `Grupo`
#'   (`Bacteria` o `Hongo`), `Genero`, `Especie` (vacía para incluir solo el
#'   género), `Clasificacion` (`Patógeno` o `Beneficioso`) y `Funcion` (opcional).
#' @return Lista con `especies`, `generos`, `resumen` (formato largo) y
#'   `muestras`, o `NULL` si la lista no tiene taxones para el marcador.
#'   Las abundancias son porcentajes de las lecturas de cada muestra.
#' @export
taxones_interes <- function(ps, tipo_datos, variable_analisis, archivo_taxones) {
  lista <- leer_taxones(archivo_taxones)
  lista <- lista[grepl(if (tipo_datos == "16S") "^bact" else "^(hongo|fung)",
                       lista$Grupo, ignore.case = TRUE), , drop = FALSE]
  if (nrow(lista) == 0) return(NULL)

  meta <- metadata_ps(ps)
  grupo <- as.character(meta[phyloseq::sample_names(ps), variable_analisis])
  muestras <- phyloseq::sample_names(ps)[order(grupo)]

  conteos <- matriz_otu(ps)[, muestras, drop = FALSE]
  pct <- sweep(conteos, 2, pmax(colSums(conteos), 1), "/") * 100
  tax <- tabla_taxonomia(ps)
  tax <- tax[rownames(pct), , drop = FALSE]
  genero_asv <- ifelse(is.na(tax$Genus), "", tax$Genus)
  epiteto_asv <- tolower(sub(".*\\s", "", ifelse(is.na(tax$Species), "", tax$Species)))

  sumar <- function(i) as.data.frame(t(round(colSums(pct[i, , drop = FALSE]), 4)),
                                     check.names = FALSE)

  # ---- Especies ----
  por_especie <- lista[nzchar(lista$Epiteto), , drop = FALSE]
  especies <- do.call(rbind, lapply(seq_len(nrow(por_especie)), function(k) {
    f <- por_especie[k, ]
    i <- coincide_genero(genero_asv, f$Genero) & epiteto_asv == f$Epiteto
    data.frame(Genero = f$Genero, Especie = paste(f$Genero, f$Epiteto),
               Clasificacion = f$Clasificacion, Funcion = f$Funcion,
               Detectada = any(i), ASVs = sum(i), sumar(i), check.names = FALSE)
  }))

  # ---- Géneros ----
  nombres <- unique(lista$Genero)
  clase_asv <- rep(NA_character_, nrow(pct))
  generos <- do.call(rbind, lapply(nombres, function(g) {
    filas <- lista[lista$Genero == g, , drop = FALSE]
    clases <- unique(filas$Clasificacion)
    categoria <- if (length(clases) > 1) "No concluyente" else paste("Potencial", tolower(clases))
    i <- coincide_genero(genero_asv, g)
    con_especie <- i & nzchar(epiteto_asv)
    listada <- con_especie & epiteto_asv %in% filas$Epiteto[nzchar(filas$Epiteto)]

    conteo <- sort(table(epiteto_asv[con_especie]), decreasing = TRUE)
    nombres_sp <- names(conteo)
    if (length(nombres_sp) > 0) {
      marca <- ifelse(nombres_sp %in% filas$Epiteto, "", " (no listada)")
      detalle <- paste0(g, " ", nombres_sp, marca, collapse = "; ")
    } else {
      detalle <- ""
    }

    data.frame(Genero = g,
               Nombre_en_base = paste(unique(genero_asv[i]), collapse = "; "),
               Categoria = categoria,
               Evidencia = if (any(listada)) "Especie listada" else if (any(i)) "Solo género" else "No detectado",
               Detectado = any(i), ASVs = sum(i), ASVs_con_especie = sum(con_especie),
               Especies_detectadas = detalle, sumar(i), check.names = FALSE)
  }))

  # Por ASV, para no contar dos veces una ASV de un género compuesto.
  for (k in seq_len(nrow(generos))) {
    i <- coincide_genero(genero_asv, generos$Genero[k])
    clase_asv[i] <- ifelse(is.na(clase_asv[i]) | clase_asv[i] == generos$Categoria[k],
                           generos$Categoria[k], "No concluyente")
  }
  categorias <- c("Potencial patógeno", "Potencial beneficioso", "No concluyente")
  resumen <- do.call(rbind, lapply(categorias, function(ct) {
    i <- !is.na(clase_asv) & clase_asv == ct
    data.frame(SampleID = muestras, Grupo = grupo[match(muestras, phyloseq::sample_names(ps))],
               Categoria = factor(ct, levels = categorias),
               Porcentaje = unname(colSums(pct[i, , drop = FALSE])))
  }))

  list(especies = especies, generos = generos, resumen = resumen, muestras = muestras)
}

#' Gráfico de burbujas de taxones de interés
#'
#' Una fila por género detectado y por especie detectada, y una columna por
#' muestra. El tamaño está en escala logarítmica para que se vean a la vez
#' abundancias de 0.01 % y de 20 %. Las celdas vacías son abundancia 0.
#'
#' @param interes Resultado de [taxones_interes()].
#' @return Un objeto ggplot, o `NULL` si no se detectó ningún taxón.
#' @export
grafico_taxones_interes <- function(interes) {
  m <- interes$muestras
  g <- interes$generos[interes$generos$Detectado, , drop = FALSE]
  e <- interes$especies[interes$especies$Detectada, , drop = FALSE]
  filas <- data.frame(Taxon = g$Genero, Categoria = g$Categoria, g[m], check.names = FALSE)
  if (nrow(e) > 0) {
    filas <- rbind(filas, data.frame(Taxon = e$Especie,
                                     Categoria = paste("Potencial", tolower(e$Clasificacion)),
                                     e[m], check.names = FALSE))
  }
  if (nrow(filas) == 0) return(NULL)

  orden <- filas$Taxon[order(rowMeans(filas[m]))]
  df <- data.frame(
    Taxon = factor(rep(filas$Taxon, times = length(m)), levels = orden),
    Categoria = factor(rep(filas$Categoria, times = length(m)),
                       levels = c("Potencial patógeno", "Potencial beneficioso", "No concluyente")),
    SampleID = factor(rep(m, each = nrow(filas)), levels = m),
    Porcentaje = unlist(filas[m], use.names = FALSE)
  )
  df$Grupo <- interes$resumen$Grupo[match(df$SampleID, interes$resumen$SampleID)]
  df <- df[df$Porcentaje > 0, , drop = FALSE]

  pastel <- RColorBrewer::brewer.pal(9, "Pastel1")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$SampleID, y = .data$Taxon,
                                   size = .data$Porcentaje, fill = .data$Categoria)) +
    ggplot2::geom_point(shape = 21, colour = "grey30") +
    ggplot2::facet_grid(Categoria ~ Grupo, scales = "free", space = "free",
                        labeller = ggplot2::label_wrap_gen(12)) +
    ggplot2::scale_size_continuous(transform = "log10", range = c(1.5, 11),
                                   breaks = c(0.01, 0.1, 1, 10),
                                   labels = c("0.01", "0.1", "1", "10")) +
    ggplot2::scale_fill_manual(values = c("Potencial patógeno" = pastel[1],
                                          "Potencial beneficioso" = pastel[3],
                                          "No concluyente" = pastel[5]),
                               guide = "none") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(face = "italic"),
                   strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::labs(x = NULL, y = NULL, size = "% de lecturas")
}

alto_taxones_interes <- function(p) {
  max(3, 1.6 + 0.32 * nlevels(droplevels(p$data$Taxon)))
}
