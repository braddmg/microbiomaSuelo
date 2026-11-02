# dada() solo aprovecha ~4 hilos por muestra (gran parte del algoritmo es
# secuencial). Con pool = FALSE las muestras son independientes, así que se
# procesan varias a la vez en procesos separados de ~4 hilos cada uno. El
# resultado es idéntico al de procesarlas una por una.
dada_muestras <- function(archivos, err, pool = FALSE, multithread = TRUE) {
  nucleos <- if (isTRUE(multithread)) parallel::detectCores()
             else if (is.numeric(multithread)) multithread else 1
  procesos <- min(length(archivos), max(1, nucleos %/% 4))

  if (!identical(pool, FALSE) || procesos < 2) {
    out <- dada2::dada(archivos, err = err, multithread = multithread, pool = pool)
    if (inherits(out, "dada")) out <- list(out)
    names(out) <- names(archivos)
    return(out)
  }

  hilos <- max(1, nucleos %/% procesos)
  message("  ", length(archivos), " muestras en ", procesos,
          " procesos paralelos de ", hilos, " hilos")
  cl <- parallel::makeCluster(procesos)
  on.exit(parallel::stopCluster(cl))
  out <- parallel::parLapplyLB(
    cl, archivos,
    function(f, err, hilos) dada2::dada(f, err = err, multithread = hilos),
    err = err, hilos = hilos
  )
  names(out) <- names(archivos)
  out
}

# DADA2 sobre fragmentos ya extendidos (single-end, sin mergePairs()).
procesar_dada2 <- function(
    archivos,
    carpeta_dada2,
    carpeta_qc,
    primer_fwd = NA,
    primer_rev = NA,
    trunc_len = 0,
    max_ee = 2,
    trunc_q = 2,
    min_len = 100,
    rm_phix = TRUE,
    pool = FALSE,
    nbases = 1e8,
    multithread = TRUE
) {
  ids <- names(archivos)

  message("Perfiles de calidad...")
  p_calidad <- dada2::plotQualityProfile(
    archivos[seq_len(min(length(archivos), 12))], n = 1e5
  )
  ggplot2::ggsave(
    file.path(carpeta_qc, "perfiles_calidad.png"), p_calidad,
    width = 12, height = 3 * ceiling(min(length(archivos), 12) / 3),
    dpi = 120, limitsize = FALSE
  )

  message("Remoción de primers...")
  primers <- remover_primers(
    archivos, file.path(carpeta_dada2, "sin_primers"),
    primer_fwd = primer_fwd, primer_rev = primer_rev
  )
  readr::write_csv(primers$resumen, file.path(carpeta_qc, "deteccion_primers.csv"))

  message("Filtrado...")
  filt_dir <- crear_directorio(file.path(carpeta_dada2, "filtered"))
  filt_files <- file.path(filt_dir, paste0(ids, "_filtered.fastq.gz"))
  names(filt_files) <- ids

  filt_stats <- dada2::filterAndTrim(
    fwd = primers$archivos,
    filt = filt_files,
    truncLen = trunc_len,
    maxEE = max_ee,
    truncQ = trunc_q,
    minLen = min_len,
    rm.phix = rm_phix,
    compress = TRUE,
    multithread = multithread,
    verbose = TRUE
  )
  # filterAndTrim nombra las filas con el archivo, no con el ID.
  filt_stats <- data.frame(
    SampleID = ids,
    input = primers$resumen$lecturas,
    filtered = filt_stats[, 2]
  )
  readr::write_csv(filt_stats, file.path(carpeta_dada2, "estadisticas_filtrado.csv"))

  buenos <- filt_files[file.exists(filt_files) & filt_stats$filtered > 0]
  if (length(buenos) < 2) {
    stop(
      "Menos de dos muestras tienen lecturas tras el filtrado. ",
      "Revise calidad y parámetros."
    )
  }

  message("Modelo de error...")
  err <- dada2::learnErrors(
    buenos, nbases = nbases, multithread = multithread, randomize = TRUE,
    verbose = TRUE
  )
  p_errores <- dada2::plotErrors(err, nominalQ = TRUE)
  ggplot2::ggsave(file.path(carpeta_qc, "modelo_error.png"), p_errores,
                  width = 10, height = 8, dpi = 120)

  message("Inferencia de ASVs...")
  dada_out <- dada_muestras(buenos, err, pool = pool, multithread = multithread)

  seqtab <- dada2::makeSequenceTable(dada_out)

  message("Eliminación de quimeras...")
  seqtab_nochim <- dada2::removeBimeraDenovo(
    seqtab, method = "consensus", multithread = multithread, verbose = TRUE
  )

  get_n <- function(x) sum(dada2::getUniques(x))
  track <- filt_stats
  track$denoised <- 0
  track$nonchim <- 0
  track$denoised[match(names(dada_out), ids)] <- vapply(dada_out, get_n, numeric(1))
  track$nonchim[match(rownames(seqtab_nochim), ids)] <- rowSums(seqtab_nochim)

  largos <- nchar(dada2::getSequences(seqtab_nochim))
  readr::write_csv(
    as.data.frame(table(largo_pb = largos), responseName = "n_ASV"),
    file.path(carpeta_dada2, "largo_ASVs.csv")
  )
  saveRDS(seqtab_nochim, file.path(carpeta_dada2, "seqtab_nochim.rds"))

  list(
    seqtab_nochim = seqtab_nochim,
    tracking = track,
    primers = primers$resumen,
    pct_quimeras = 100 * (1 - sum(seqtab_nochim) / sum(seqtab))
  )
}
