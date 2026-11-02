#' @importFrom ggplot2 .data
NULL

RANGOS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Pastel1 y luego Pastel2, sin sus grises (reservados a "Otros" y "No
# clasificado"). Si se piden más de 15 colores, se interpolan.
paleta_pastel <- function(n) {
  colores <- c(RColorBrewer::brewer.pal(9, "Pastel1")[-9],
               RColorBrewer::brewer.pal(8, "Pastel2")[-8])
  if (n > length(colores)) return(grDevices::colorRampPalette(colores)(n))
  colores[seq_len(n)]
}

crear_directorio <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

matriz_otu <- function(ps) {
  m <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

metadata_ps <- function(ps) {
  methods::as(phyloseq::sample_data(ps), "data.frame")
}

tabla_taxonomia <- function(ps) {
  as.data.frame(methods::as(phyloseq::tax_table(ps), "matrix"),
                stringsAsFactors = FALSE)
}

# Etiqueta legible de un rango. Para "Species" se arma el binomio
# (SILVA/addSpecies solo devuelve el epíteto; UNITE ya trae el binomio).
etiqueta_taxon <- function(tax, rank) {
  x <- tax[[rank]]
  if (rank == "Species" && "Genus" %in% names(tax)) {
    g <- tax$Genus
    binomio <- !is.na(x) & !is.na(g) & !startsWith(x, g)
    x[binomio] <- paste(g[binomio], x[binomio])
  }
  x
}

# Busca pandoc (necesario para el reporte) cuando R no corre dentro de RStudio.
asegurar_pandoc <- function() {
  if (rmarkdown::pandoc_available()) return(invisible(TRUE))
  candidatos <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/quarto/bin/tools",
    "C:/Program Files/Quarto/bin/tools",
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
    "/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"
  )
  for (d in candidatos) {
    if (dir.exists(d)) {
      Sys.setenv(RSTUDIO_PANDOC = d)
      rmarkdown::find_pandoc(cache = FALSE)
      if (rmarkdown::pandoc_available()) return(invisible(TRUE))
    }
  }
  stop(
    "No se encontró pandoc, necesario para generar el reporte HTML. ",
    "Ejecute el análisis desde RStudio o instale pandoc (https://pandoc.org)."
  )
}
