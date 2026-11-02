# Primers por defecto según el marcador.
#   16S V3-V4: 341F / 805R      ITS1: ITS1 / ITS2
PRIMERS <- list(
  `16S` = c(fwd = "CCTACGGGNGGCWGCAG", rev = "GACTACHVGGGTATCTAATCC"),
  ITS = c(fwd = "TCCGTAGGTGAACCTGCGG", rev = "GCTGCGTTCTTCATCGATGC")
)

# Ruta de análisis según el marcador.
CONFIG_MARCADOR <- list(
  `16S` = list(min_len = 100, filtro_organelos = TRUE,
               base_esperada = "silva|rdp|gtdb|greengenes",
               base_incorrecta = "unite|sh_general"),
  ITS = list(min_len = 50, filtro_organelos = FALSE,
             base_esperada = "unite|sh_general",
             base_incorrecta = "silva|rdp|gtdb|greengenes")
)

# Posición de la primera coincidencia de `patron` en cada lectura (0 = sin
# coincidencia). Admite códigos IUPAC y `max_mismatch` sustituciones.
buscar_primer <- function(lecturas, patron, max_mismatch) {
  m <- Biostrings::vmatchPattern(patron, lecturas, max.mismatch = max_mismatch,
                                 fixed = FALSE)
  hay <- lengths(Biostrings::startIndex(m)) > 0
  ini <- fin <- integer(length(lecturas))
  ini[hay] <- vapply(Biostrings::startIndex(m)[hay], `[`, integer(1), 1)
  fin[hay] <- vapply(Biostrings::endIndex(m)[hay], `[`, integer(1), 1)
  list(hay = hay, inicio = ini, fin = fin)
}

#' Remover primers de fragmentos extendidos
#'
#' Busca el primer forward al inicio de cada lectura y el complemento
#' reverso del primer reverse al final, y recorta todo lo que esté antes del
#' forward y desde el reverse en adelante. Las lecturas en las que no se
#' encuentra un primer se conservan sin recortar en ese extremo.
#'
#' @param archivos Vector nombrado (IDs de muestra) de archivos FASTQ.
#' @param carpeta Carpeta donde escribir los FASTQ recortados.
#' @param primer_fwd,primer_rev Secuencias 5'->3' (IUPAC). `NA` para omitir.
#' @param ventana Bases de cada extremo donde se busca el primer.
#' @param max_mismatch Sustituciones permitidas.
#' @param min_deteccion Si un primer aparece en menos de esta proporción de
#'   lecturas de una muestra, no se recorta en esa muestra.
#' @return Lista con `archivos` (rutas a usar después) y `resumen` (data.frame
#'   con el porcentaje de lecturas con cada primer por muestra).
#' @export
remover_primers <- function(archivos, carpeta, primer_fwd, primer_rev,
                            ventana = 40, max_mismatch = 2,
                            min_deteccion = 0.05) {
  crear_directorio(carpeta)
  rc_rev <- if (!is.na(primer_rev)) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer_rev)))
  }
  salida <- archivos
  resumen <- vector("list", length(archivos))

  for (i in seq_along(archivos)) {
    fq <- ShortRead::readFastq(archivos[i])
    s <- ShortRead::sread(fq)
    w <- Biostrings::width(s)
    inicio <- rep(1L, length(s))
    fin <- w
    pct_fwd <- pct_rev <- NA_real_

    if (!is.na(primer_fwd)) {
      b <- buscar_primer(Biostrings::subseq(s, 1, pmin(w, ventana)), primer_fwd, max_mismatch)
      pct_fwd <- mean(b$hay)
      if (pct_fwd >= min_deteccion) inicio[b$hay] <- b$fin[b$hay] + 1L
    }
    if (!is.na(primer_rev)) {
      desde <- pmax(1L, w - ventana + 1L)
      b <- buscar_primer(Biostrings::subseq(s, desde, w), rc_rev, max_mismatch)
      pct_rev <- mean(b$hay)
      if (pct_rev >= min_deteccion) fin[b$hay] <- desde[b$hay] + b$inicio[b$hay] - 2L
    }

    recortar <- any(inicio > 1L | fin < w)
    if (recortar) {
      ok <- fin - inicio + 1L >= 20L
      fq <- ShortRead::narrow(fq[ok], start = inicio[ok], end = fin[ok])
      salida[i] <- file.path(carpeta, paste0(names(archivos)[i], "_sin_primers.fastq.gz"))
      if (file.exists(salida[i])) file.remove(salida[i])
      ShortRead::writeFastq(fq, salida[i], compress = TRUE)
    }

    resumen[[i]] <- data.frame(
      SampleID = names(archivos)[i],
      lecturas = length(s),
      pct_primer_fwd = round(100 * pct_fwd, 1),
      pct_primer_rev = round(100 * pct_rev, 1),
      recortado = recortar
    )
  }

  resumen <- do.call(rbind, resumen)
  message("Primers detectados (% de lecturas): fwd ",
          paste(resumen$pct_primer_fwd, collapse = ", "), " | rev ",
          paste(resumen$pct_primer_rev, collapse = ", "))
  list(archivos = salida, resumen = resumen)
}
