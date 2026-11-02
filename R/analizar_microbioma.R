#' Análisis completo de amplicones 16S/ITS de muestras de suelo
#'
#' Procesa todos los archivos `*.extendedFrags.fastq.gz` (fragmentos ya
#' extendidos con FLASH) de una carpeta y genera tablas, figuras, gráficos
#' Krona y un reporte HTML (`Reporte_final.html`).
#'
#' El ID de cada muestra es el nombre del archivo antes de
#' `.extendedFrags.fastq.gz`. En la metadata, la columna `columna_id` puede
#' contener ese ID con o sin el sufijo `.extendedFrags`.
#'
#' `tipo_datos` define la ruta del análisis:
#'
#' | | 16S | ITS |
#' |---|---|---|
#' | Primers por defecto | 341F / 805R | ITS1 / ITS2 |
#' | Largo mínimo (`min_len`) | 100 pb | 50 pb |
#' | Base de referencia | SILVA (+ `addSpecies`) | UNITE (`.tgz` o FASTA) |
#' | Filtro taxonómico | sin reino, eucariotas, cloroplastos, mitocondrias | sin reino |
#'
#' @param tipo_datos `"16S"` o `"ITS"`.
#' @param carpeta_entrada Carpeta con los archivos `*.extendedFrags.fastq.gz`.
#' @param archivo_metadata CSV (separado por `,` o `;`) con una fila por muestra.
#' @param variable_analisis Columna de la metadata con los grupos a comparar.
#' @param ruta_taxonomia Base de referencia para `dada2::assignTaxonomy()`.
#'   Acepta `.fa`, `.fa.gz` o el `.tgz` de UNITE tal como se descarga.
#' @param ruta_especies Base para `dada2::addSpecies()` (opcional, solo 16S).
#' @param columna_id Columna de la metadata con el ID de muestra.
#' @param carpeta_salida Carpeta donde se escriben los resultados.
#' @param primer_fwd,primer_rev Primers (5'->3', IUPAC) a recortar. `NULL` usa
#'   los del marcador; `NA` desactiva el recorte de ese primer. Si un primer no
#'   aparece en las lecturas (ya fue removido), no se recorta nada.
#' @param trunc_len,max_ee,trunc_q Parámetros de `dada2::filterAndTrim()`.
#' @param min_len Largo mínimo tras el filtrado; `NULL` usa el del marcador.
#' @param min_boot Bootstrap mínimo para asignar un rango taxonómico.
#' @param min_lecturas_asv Se eliminan las ASVs con menos de estas lecturas,
#'   sumando todas las muestras.
#' @param pool Agrupamiento de muestras en `dada2::dada()` (`FALSE`, `TRUE` o `"pseudo"`).
#' @param nbases Bases usadas para aprender el modelo de error
#'   (`dada2::learnErrors()`). Reducirlo (p. ej. `5e7`) acorta mucho esta etapa
#'   con fragmentos largos, a costa de un modelo algo menos preciso.
#' @param multithread Usar varios núcleos.
#' @param permutations Permutaciones para PERMANOVA y dispersión.
#' @param crear_krona_plots Generar gráficos Krona.
#' @param archivo_taxones CSV del usuario con los taxones de interés (ver
#'   [taxones_interes()]). Con `NULL` esa sección no se incluye en el reporte.
#' @param reanudar Si `TRUE`, reutiliza los resultados de DADA2 y taxonomía de
#'   una ejecución previa en la misma carpeta, siempre que las entradas y los
#'   parámetros sean idénticos.
#' @param abrir_reporte Abrir el reporte en el navegador al terminar.
#'
#' @return Lista (invisible) con el objeto phyloseq, tablas y ruta del reporte.
#' @export
#'
#' @examples
#' \dontrun{
#' analizar_microbioma(
#'   tipo_datos = "16S",
#'   carpeta_entrada = "Data 16s",
#'   archivo_metadata = "mdata_16S.csv",
#'   variable_analisis = "Tratamiento",
#'   ruta_taxonomia = "Databases/silva_nr99_v138.2_toGenus_trainset.fa.gz",
#'   ruta_especies = "Databases/silva_v138.2_assignSpecies.fa.gz",
#'   carpeta_salida = "Resultados_16S"
#' )
#'
#' analizar_microbioma(
#'   tipo_datos = "ITS",
#'   carpeta_entrada = "Data ITS",
#'   archivo_metadata = "mdata_ITS.csv",
#'   variable_analisis = "Tratamiento",
#'   ruta_taxonomia = "Databases/sh_general_release_19.02.2025.tgz",
#'   carpeta_salida = "Resultados_ITS"
#' )
#' }
analizar_microbioma <- function(
    tipo_datos = c("16S", "ITS"),
    carpeta_entrada,
    archivo_metadata,
    variable_analisis,
    ruta_taxonomia,
    ruta_especies = NULL,
    columna_id = "SampleID",
    carpeta_salida = file.path(getwd(), paste0("Resultados_", tipo_datos[1])),
    primer_fwd = NULL,
    primer_rev = NULL,
    trunc_len = 0,
    max_ee = 2,
    trunc_q = 2,
    min_len = NULL,
    min_boot = 50,
    min_lecturas_asv = 5,
    pool = FALSE,
    nbases = 1e8,
    multithread = TRUE,
    permutations = 999,
    crear_krona_plots = TRUE,
    archivo_taxones = NULL,
    reanudar = TRUE,
    abrir_reporte = interactive()
) {
  inicio <- Sys.time()

  # ---- Ruta según el marcador ----
  tipo_datos <- match.arg(tipo_datos)
  ruta <- CONFIG_MARCADOR[[tipo_datos]]
  if (is.null(primer_fwd)) primer_fwd <- PRIMERS[[tipo_datos]][["fwd"]]
  if (is.null(primer_rev)) primer_rev <- PRIMERS[[tipo_datos]][["rev"]]
  if (is.null(min_len)) min_len <- ruta$min_len

  if (grepl(ruta$base_incorrecta, basename(ruta_taxonomia), ignore.case = TRUE)) {
    stop("La base '", basename(ruta_taxonomia), "' no corresponde a datos ",
         tipo_datos, ". Revise tipo_datos o ruta_taxonomia.")
  }
  if (!grepl(ruta$base_esperada, basename(ruta_taxonomia), ignore.case = TRUE)) {
    warning("No se reconoce '", basename(ruta_taxonomia), "' como base para ",
            tipo_datos, "; verifique que sea la correcta.")
  }
  if (tipo_datos == "ITS" && !is.null(ruta_especies)) {
    warning("ruta_especies solo se usa con 16S; se ignora para ITS.")
    ruta_especies <- NULL
  }
  message("Análisis ", tipo_datos, ": primers ", primer_fwd, " / ", primer_rev,
          ", largo mínimo ", min_len, " pb",
          if (ruta$filtro_organelos) ", sin cloroplastos ni mitocondrias", ".")

  entradas <- validar_entradas(
    carpeta_entrada = carpeta_entrada,
    tipo_datos = tipo_datos,
    archivo_metadata = archivo_metadata,
    variable_analisis = variable_analisis,
    columna_id = columna_id
  )
  tipo_datos <- entradas$tipo_datos
  metadata <- entradas$metadata

  archivos <- detectar_fastq(carpeta_entrada)
  ids_sin_metadata <- setdiff(names(archivos), metadata[[columna_id]])
  if (length(ids_sin_metadata) > 0) {
    stop(
      "No se encontraron en la metadata estos IDs: ",
      paste(ids_sin_metadata, collapse = ", "),
      "\nIDs en la metadata: ", paste(metadata[[columna_id]], collapse = ", ")
    )
  }

  ruta_taxonomia <- preparar_referencia(ruta_taxonomia)
  ruta_especies <- preparar_referencia(ruta_especies)
  if (!is.null(archivo_taxones)) leer_taxones(archivo_taxones)  # falla aquí y no tras horas de DADA2

  carpeta_salida <- crear_directorio(carpeta_salida)
  dirs <- list(
    qc = crear_directorio(file.path(carpeta_salida, "01_Control_calidad")),
    dada2 = crear_directorio(file.path(carpeta_salida, "02_DADA2")),
    tax = crear_directorio(file.path(carpeta_salida, "03_Taxonomia")),
    abund = crear_directorio(file.path(carpeta_salida, "04_Abundancia")),
    krona = crear_directorio(file.path(carpeta_salida, "05_Krona")),
    heatmaps = crear_directorio(file.path(carpeta_salida, "06_Heatmaps")),
    diversity = crear_directorio(file.path(carpeta_salida, "07_Diversidad")),
    tables = crear_directorio(file.path(carpeta_salida, "08_Tablas")),
    reporte = crear_directorio(file.path(carpeta_salida, "09_Reporte"))
  )

  config <- list(
    tipo_datos = tipo_datos,
    carpeta_entrada = normalizePath(carpeta_entrada, winslash = "/"),
    archivo_metadata = normalizePath(archivo_metadata, winslash = "/"),
    variable_analisis = variable_analisis,
    columna_id = columna_id,
    ruta_taxonomia = ruta_taxonomia,
    ruta_especies = ruta_especies,
    primer_fwd = primer_fwd,
    primer_rev = primer_rev,
    trunc_len = trunc_len,
    max_ee = max_ee,
    trunc_q = trunc_q,
    min_len = min_len,
    min_boot = min_boot,
    min_lecturas_asv = min_lecturas_asv,
    pool = pool,
    permutations = permutations,
    rarefaccion = FALSE,
    version_paquete = as.character(utils::packageVersion("microbiomaSuelo")),
    version_dada2 = as.character(utils::packageVersion("dada2")),
    fecha_ejecucion = format(inicio)
  )
  saveRDS(config, file.path(carpeta_salida, "configuracion_analisis.rds"))
  capture.output(str(config), file = file.path(carpeta_salida, "configuracion_analisis.txt"))

  # ---- DADA2 (con caché) ----
  clave_dada2 <- list(archivos = basename(archivos), tamanos = file.size(archivos),
                      primer_fwd = primer_fwd, primer_rev = primer_rev,
                      trunc_len = trunc_len, max_ee = max_ee, trunc_q = trunc_q,
                      min_len = min_len, pool = pool)
  cache_dada2 <- file.path(dirs$dada2, "cache_dada2.rds")
  dada_res <- NULL
  if (reanudar && file.exists(cache_dada2)) {
    previo <- readRDS(cache_dada2)
    if (identical(previo$clave, clave_dada2)) {
      message("Reutilizando resultados de DADA2 de una ejecución previa.")
      dada_res <- previo$resultado
    }
  }
  if (is.null(dada_res)) {
    dada_res <- procesar_dada2(
      archivos = archivos,
      carpeta_dada2 = dirs$dada2,
      carpeta_qc = dirs$qc,
      primer_fwd = primer_fwd,
      primer_rev = primer_rev,
      trunc_len = trunc_len,
      max_ee = max_ee,
      trunc_q = trunc_q,
      min_len = min_len,
      pool = pool,
      nbases = nbases,
      multithread = multithread
    )
    saveRDS(list(clave = clave_dada2, resultado = dada_res), cache_dada2)
  }

  # ---- Taxonomía (con caché) ----
  seqs <- colnames(dada_res$seqtab_nochim)
  clave_tax <- list(seqs = seqs, ref = basename(ruta_taxonomia),
                    ref_tamano = file.size(ruta_taxonomia),
                    especies = if (!is.null(ruta_especies)) basename(ruta_especies),
                    min_boot = min_boot)
  cache_tax <- file.path(dirs$tax, "cache_taxonomia.rds")
  tax <- NULL
  if (reanudar && file.exists(cache_tax)) {
    previo <- readRDS(cache_tax)
    if (identical(previo$clave, clave_tax)) {
      message("Reutilizando la asignación taxonómica previa.")
      tax <- previo$tax
    }
  }
  if (is.null(tax)) {
    tax <- asignar_taxonomia(seqs, tipo_datos, ruta_taxonomia, ruta_especies,
                             min_boot = min_boot, multithread = multithread)
    saveRDS(list(clave = clave_tax, tax = tax), cache_tax)
  }

  # ---- phyloseq y tablas ----
  ps <- construir_phyloseq(dada_res$seqtab_nochim, tax, metadata,
                           columna_id = columna_id, tipo_datos = tipo_datos,
                           min_lecturas_asv = min_lecturas_asv)
  saveRDS(ps, file.path(carpeta_salida, "phyloseq_final.rds"))
  readr::write_csv(
    data.frame(ASV = phyloseq::taxa_names(ps), tabla_taxonomia(ps), check.names = FALSE),
    file.path(dirs$tax, "taxonomia.csv")
  )

  tracking <- dada_res$tracking
  tracking$final <- 0
  tracking$final[match(phyloseq::sample_names(ps), tracking$SampleID)] <-
    phyloseq::sample_sums(ps)
  tracking$pct_retenido <- round(100 * tracking$final / tracking$input, 1)
  readr::write_csv(tracking, file.path(dirs$qc, "tracking_lecturas.csv"))

  exportar_tablas(ps, dirs$tables, dirs$abund)

  # ---- Figuras y estadística ----
  message("Figuras y diversidad...")
  composicion <- list(
    Phylum = grafico_composicion(ps, "Phylum", variable_analisis, top_n = 10),
    Genus = grafico_composicion(ps, "Genus", variable_analisis, top_n = 15)
  )
  for (r in names(composicion)) {
    ggplot2::ggsave(file.path(dirs$abund, paste0("composicion_", tolower(r), ".png")),
                    composicion[[r]], width = 10, height = 6, dpi = 150)
  }

  heatmaps <- list(
    Genus = heatmap_top(ps, "Genus", variable_analisis, 20,
                        file.path(dirs$heatmaps, "heatmap_top20_genus.png")),
    Species = heatmap_top(ps, "Species", variable_analisis, 20,
                          file.path(dirs$heatmaps, "heatmap_top20_species.png"))
  )

  interes <- NULL
  interes_plot <- NULL
  if (!is.null(archivo_taxones)) {
    interes <- taxones_interes(ps, tipo_datos, variable_analisis, archivo_taxones)
  }
  if (!is.null(interes)) {
    dirs$interes <- crear_directorio(file.path(carpeta_salida, "10_Taxones_interes"))
    interes_plot <- grafico_taxones_interes(interes)
    readr::write_csv(interes$especies, file.path(dirs$interes, "especies_interes.csv"))
    readr::write_csv(interes$generos, file.path(dirs$interes, "generos_interes.csv"))
    readr::write_csv(interes$resumen, file.path(dirs$interes, "resumen_por_categoria.csv"))
    if (!is.null(interes_plot)) {
      ggplot2::ggsave(file.path(dirs$interes, "abundancia_taxones_interes.png"), interes_plot,
                      width = 9, height = alto_taxones_interes(interes_plot), dpi = 150)
    }
  }

  alfa <- calcular_alfa(ps, variable_analisis)
  alfa_plot <- grafico_alfa(alfa, variable_analisis)
  alfa_test <- pruebas_alfa(alfa, variable_analisis)
  readr::write_csv(alfa, file.path(dirs$diversity, "alfa_diversidad_por_muestra.csv"))
  if (!is.null(alfa_test)) {
    readr::write_csv(alfa_test, file.path(dirs$diversity, "alfa_kruskal_wallis.csv"))
  }
  ggplot2::ggsave(file.path(dirs$diversity, "alfa_diversidad.png"), alfa_plot,
                  width = 10, height = 4.5, dpi = 150)

  beta <- analisis_beta(ps, variable_analisis, permutations = permutations)
  ggplot2::ggsave(file.path(dirs$diversity, "pcoa_bray_curtis.png"), beta$plot,
                  width = 8, height = 6, dpi = 150)
  readr::write_csv(beta$coordinates, file.path(dirs$diversity, "pcoa_coordenadas.csv"))
  if (!is.null(beta$permanova)) {
    utils::write.csv(beta$permanova, file.path(dirs$diversity, "permanova_bray_curtis.csv"))
    utils::write.csv(beta$dispersion, file.path(dirs$diversity, "homogeneidad_dispersion.csv"))
  }

  krona <- NULL
  if (isTRUE(crear_krona_plots)) {
    message("Gráficos Krona...")
    krona <- crear_krona(ps, dirs$krona, tipo_datos)
  }

  # ---- Reporte ----
  message("Reporte HTML...")
  reporte <- crear_reporte_html(
    resultados = list(
      config = config,
      metadata = metadata,
      tracking = tracking,
      primers = dada_res$primers,
      pct_quimeras = dada_res$pct_quimeras,
      n_asv = phyloseq::ntaxa(ps),
      clasificacion = colMeans(!is.na(tabla_taxonomia(ps))) * 100,
      composicion = composicion,
      heatmaps = heatmaps,
      interes = interes,
      interes_plot = interes_plot,
      alfa = alfa,
      alfa_plot = alfa_plot,
      alfa_test = alfa_test,
      beta = beta[c("plot", "permanova", "dispersion")],
      krona = krona,
      img_calidad = file.path(dirs$qc, "perfiles_calidad.png")
    ),
    carpeta_salida = carpeta_salida,
    carpeta_reporte = dirs$reporte
  )

  message("Análisis completado en ",
          format(round(difftime(Sys.time(), inicio, units = "mins"), 1)), ".")
  message("Reporte: ", reporte)
  if (isTRUE(abrir_reporte)) utils::browseURL(reporte)

  invisible(list(
    phyloseq = ps,
    tracking = tracking,
    alpha = alfa,
    beta = beta,
    krona = krona,
    output_dir = carpeta_salida,
    report = reporte
  ))
}
