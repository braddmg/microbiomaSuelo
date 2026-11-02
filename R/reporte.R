crear_reporte_html <- function(resultados, carpeta_salida, carpeta_reporte) {
  asegurar_pandoc()

  datos <- file.path(carpeta_reporte, "datos_reporte.rds")
  saveRDS(resultados, datos)

  # Se renderiza en una carpeta temporal para no escribir en la librería
  # de R ni dejar archivos intermedios en la carpeta de resultados.
  tmp <- tempfile("reporte_")
  dir.create(tmp)
  plantilla <- system.file("rmarkdown", "reporte.Rmd", package = "microbiomaSuelo")
  file.copy(plantilla, tmp)

  rmarkdown::render(
    input = file.path(tmp, "reporte.Rmd"),
    output_file = "Reporte_final.html",
    output_dir = carpeta_salida,
    intermediates_dir = tmp,
    knit_root_dir = tmp,
    params = list(datos = normalizePath(datos, winslash = "/")),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  file.path(carpeta_salida, "Reporte_final.html")
}

#' Copiar el tutorial paso a paso
#'
#' Copia el R Markdown que explica cada paso del análisis y reproduce las
#' figuras del reporte. Hay una versión para 16S y otra para ITS.
#'
#' @param tipo_datos `"16S"` o `"ITS"`.
#' @param destino Carpeta destino (por defecto, el directorio de trabajo).
#' @return Ruta del archivo copiado (invisible).
#' @export
copiar_tutorial <- function(tipo_datos = c("16S", "ITS"), destino = getwd()) {
  tipo_datos <- match.arg(tipo_datos)
  origen <- system.file("tutorial", paste0("tutorial_", tipo_datos, ".Rmd"),
                        package = "microbiomaSuelo")
  archivo <- file.path(destino, basename(origen))
  if (file.exists(archivo)) stop("Ya existe: ", archivo)
  file.copy(origen, archivo)
  message("Tutorial copiado en: ", archivo)
  invisible(archivo)
}
