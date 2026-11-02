# Los archivos de entrada son fragmentos ya extendidos con FLASH.
SUFIJO_FASTQ <- "\\.extendedFrags\\.fastq\\.gz$"

#' ID de muestra a partir del nombre del archivo
#'
#' El ID es todo lo que precede a `.extendedFrags.fastq.gz`.
#'
#' @param archivo Ruta(s) a archivos FASTQ.
#' @return Vector de caracteres con los IDs.
#' @examples
#' extraer_id_muestra("Data 16s/R1SG3006.extendedFrags.fastq.gz")
#' @export
extraer_id_muestra <- function(archivo) {
  sub(SUFIJO_FASTQ, "", basename(archivo))
}

# Permite que la metadata use "R1SG3006" o "R1SG3006.extendedFrags".
limpiar_id <- function(x) {
  sub("\\.extendedFrags(\\.fastq\\.gz)?$", "", trimws(as.character(x)))
}

detectar_fastq <- function(carpeta) {
  archivos <- list.files(carpeta, pattern = SUFIJO_FASTQ, full.names = TRUE)

  if (length(archivos) == 0) {
    stop("No se encontraron archivos *.extendedFrags.fastq.gz en: ", carpeta)
  }

  todos <- list.files(carpeta, pattern = "\\.(fastq|fq)(\\.gz)?$",
                      ignore.case = TRUE)
  ignorados <- setdiff(todos, basename(archivos))
  if (length(ignorados) > 0) {
    warning(
      "Se ignoran FASTQ que no terminan en .extendedFrags.fastq.gz ",
      "(p. ej. notCombined de FLASH): ", paste(ignorados, collapse = ", ")
    )
  }

  names(archivos) <- extraer_id_muestra(archivos)
  if (anyDuplicated(names(archivos))) {
    stop("Hay IDs de muestra duplicados: ",
         paste(unique(names(archivos)[duplicated(names(archivos))]),
               collapse = ", "))
  }
  archivos
}

leer_metadata <- function(archivo) {
  # Excel en español suele exportar CSV separados por ";".
  primera <- readLines(archivo, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (lengths(regmatches(primera, gregexpr(";", primera))) >
               lengths(regmatches(primera, gregexpr(",", primera)))) ";" else ","
  readr::read_delim(
    archivo, delim = delim, show_col_types = FALSE,
    name_repair = "minimal", trim_ws = TRUE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    as.data.frame(check.names = FALSE)
}

validar_entradas <- function(carpeta_entrada, tipo_datos,
                             archivo_metadata, variable_analisis,
                             columna_id = "SampleID") {
  if (!dir.exists(carpeta_entrada)) {
    stop("No existe la carpeta de entrada: ", carpeta_entrada)
  }
  if (!file.exists(archivo_metadata)) {
    stop("No existe el archivo de metadata: ", archivo_metadata)
  }

  tipo_datos <- match.arg(tipo_datos, choices = c("16S", "ITS"))
  metadata <- leer_metadata(archivo_metadata)

  if (anyDuplicated(names(metadata))) {
    stop("La metadata tiene nombres de columnas duplicados.")
  }
  if (!columna_id %in% names(metadata)) {
    stop(
      "La metadata debe contener la columna identificadora '",
      columna_id, "'. Columnas disponibles: ",
      paste(names(metadata), collapse = ", ")
    )
  }
  if (!variable_analisis %in% names(metadata)) {
    stop(
      "No existe la variable '", variable_analisis,
      "' en la metadata. Columnas disponibles: ",
      paste(names(metadata), collapse = ", ")
    )
  }

  metadata[[columna_id]] <- limpiar_id(metadata[[columna_id]])
  if (any(is.na(metadata[[columna_id]]) | metadata[[columna_id]] == "")) {
    stop("Hay IDs de muestra vacíos en la metadata.")
  }
  if (anyDuplicated(metadata[[columna_id]])) {
    stop("La columna ", columna_id, " contiene IDs duplicados.")
  }
  if (all(is.na(metadata[[variable_analisis]]))) {
    stop("La variable de análisis contiene solamente valores NA.")
  }
  metadata[[variable_analisis]] <- factor(
    metadata[[variable_analisis]],
    levels = unique(stats::na.omit(metadata[[variable_analisis]]))
  )

  list(
    tipo_datos = tipo_datos,
    metadata = metadata,
    columna_id = columna_id,
    variable_analisis = variable_analisis
  )
}
