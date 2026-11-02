# Gráficos Krona generados directamente desde R. No requiere KronaTools ni
# Perl: se escribe el mismo HTML que produce `ktImportText` en modo
# standalone, con el JavaScript de Krona embebido (inst/krona).

escapar_xml <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

# Linaje de cada ASV hasta el último rango clasificado. Las lecturas no
# clasificadas a un nivel quedan en el nodo padre (Krona las muestra como
# porción sin asignar).
linajes_krona <- function(ps) {
  tax <- tabla_taxonomia(ps)
  tax$Species <- etiqueta_taxon(tax, "Species")
  rangos <- intersect(RANGOS, names(tax))
  lapply(seq_len(nrow(tax)), function(i) {
    x <- unlist(tax[i, rangos], use.names = FALSE)
    corte <- which(is.na(x) | x == "")[1]
    if (!is.na(corte)) x <- x[seq_len(corte - 1)]
    x
  })
}

nodos_krona <- function(nombre, filas, nivel, linajes, conteos) {
  vals <- colSums(conteos[filas, , drop = FALSE])
  profundos <- filas[lengths(linajes[filas]) >= nivel]
  hijos <- if (length(profundos) > 0) {
    split(profundos, vapply(linajes[profundos], `[[`, character(1), nivel))
  } else {
    list()
  }
  hijos <- hijos[order(-vapply(hijos, function(f) sum(conteos[f, ]), numeric(1)))]
  c(
    sprintf('<node name="%s"><magnitude>%s</magnitude>', escapar_xml(nombre),
            paste0("<val>", round(vals), "</val>", collapse = "")),
    unlist(lapply(names(hijos), function(h) {
      nodos_krona(h, hijos[[h]], nivel + 1, linajes, conteos)
    })),
    "</node>"
  )
}

escribir_krona <- function(linajes, conteos, archivo, titulo, raiz = "Total") {
  plantilla <- system.file("krona", "krona_header.html", package = "microbiomaSuelo")
  header <- readLines(plantilla, warn = FALSE, encoding = "UTF-8")
  header <- sub("{{TITULO}}", escapar_xml(titulo), header, fixed = TRUE)

  cuerpo <- c(
    '<krona collapse="false" key="true">',
    '<attributes magnitude="magnitude">',
    '<attribute display="Lecturas">magnitude</attribute>',
    "</attributes>",
    "<datasets>",
    paste0("<dataset>", escapar_xml(colnames(conteos)), "</dataset>"),
    "</datasets>",
    nodos_krona(raiz, seq_len(nrow(conteos)), 1, linajes, conteos),
    "</krona>",
    "</div>",
    "</body>",
    "</html>"
  )

  con <- file(archivo, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(c(header, cuerpo), con)
  invisible(archivo)
}

#' Gráficos Krona por muestra
#'
#' Genera un HTML interactivo por muestra y uno combinado (`todas_las_muestras.html`)
#' en el que se cambia de muestra con el selector de la izquierda. Los
#' archivos funcionan sin conexión a internet.
#'
#' @param ps Objeto phyloseq.
#' @param carpeta_krona Carpeta de salida.
#' @param tipo_datos Texto para el nodo raíz (p. ej. "16S").
#' @return data.frame con `SampleID`, `Lecturas` y `Archivo` (nombre del HTML).
#' @export
crear_krona <- function(ps, carpeta_krona, tipo_datos = "") {
  crear_directorio(carpeta_krona)
  conteos <- matriz_otu(ps)
  linajes <- linajes_krona(ps)
  raiz <- trimws(paste("Total", tipo_datos))

  archivos <- paste0(gsub("[^A-Za-z0-9._-]", "_", colnames(conteos)), ".html")
  for (j in seq_len(ncol(conteos))) {
    filas <- which(conteos[, j] > 0)
    escribir_krona(
      linajes[filas], conteos[filas, j, drop = FALSE],
      file.path(carpeta_krona, archivos[j]),
      titulo = paste("Krona -", colnames(conteos)[j]), raiz = raiz
    )
  }
  escribir_krona(linajes, conteos,
                 file.path(carpeta_krona, "todas_las_muestras.html"),
                 titulo = "Krona - todas las muestras", raiz = raiz)

  data.frame(SampleID = colnames(conteos), Lecturas = colSums(conteos),
             Archivo = archivos, row.names = NULL)
}
