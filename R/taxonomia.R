# Resuelve la ruta de la base de referencia:
#  - si "x.fa" no existe pero "x.fa.gz" sí, usa la comprimida;
#  - si es un .tgz/.tar.gz (p. ej. UNITE), extrae el FASTA (no el "_dev").
preparar_referencia <- function(ruta) {
  if (is.null(ruta)) return(NULL)
  if (!file.exists(ruta) && file.exists(paste0(ruta, ".gz"))) {
    ruta <- paste0(ruta, ".gz")
  }
  if (!file.exists(ruta)) stop("No existe la base de referencia: ", ruta)

  if (grepl("\\.(tgz|tar\\.gz)$", ruta)) {
    contenido <- utils::untar(ruta, list = TRUE)
    fastas <- grep("\\.(fa|fasta|fna)(\\.gz)?$", contenido, value = TRUE)
    if (length(fastas) == 0) stop("No hay archivos FASTA dentro de: ", ruta)
    elegido <- c(fastas[!grepl("_dev\\.", fastas)], fastas)[1]
    destino <- file.path(dirname(ruta), sub("\\.(tgz|tar\\.gz)$", "", basename(ruta)))
    if (!file.exists(file.path(destino, elegido))) {
      message("Extrayendo ", elegido, "...")
      utils::untar(ruta, files = elegido, exdir = destino)
    }
    ruta <- file.path(destino, elegido)
  }
  normalizePath(ruta, winslash = "/")
}

normalizar_taxonomia <- function(tax, tipo_datos) {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)

  for (r in RANGOS) {
    if (!r %in% names(tax)) tax[[r]] <- NA_character_
    x <- as.character(tax[[r]])
    # UNITE: prefijos k__, p__ ... y rangos "unidentified".
    x <- sub("^[a-z]__", "", x)
    x[x == "" | grepl("^unidentified|_sp$", x, ignore.case = TRUE)] <- NA_character_
    if (r == "Species") x <- gsub("_", " ", x)
    tax[[r]] <- x
  }

  tax[, RANGOS]
}

asignar_taxonomia <- function(seqs, tipo_datos, ruta_taxonomia,
                              ruta_especies = NULL, min_boot = 50,
                              multithread = TRUE) {
  message("Asignación taxonómica (puede tardar varios minutos)...")
  tax <- dada2::assignTaxonomy(
    seqs,
    refFasta = ruta_taxonomia,
    minBoot = min_boot,
    tryRC = TRUE,
    multithread = multithread
  )

  if (!is.null(ruta_especies)) {
    tax <- dada2::addSpecies(
      tax, refFasta = ruta_especies, allowMultiple = FALSE, tryRC = TRUE
    )
  }

  tax <- normalizar_taxonomia(tax, tipo_datos)
  rownames(tax) <- seqs
  tax
}
