#' Barras de composición taxonómica por muestra
#'
#' @param ps Objeto phyloseq.
#' @param rank Rango taxonómico.
#' @param variable_analisis Columna de la metadata para agrupar las muestras.
#' @param top_n Número de taxones a mostrar; el resto se agrupa en "Otros".
#' @return Un objeto ggplot.
#' @export
grafico_composicion <- function(ps, rank, variable_analisis, top_n = 10) {
  m <- tabla_rango(ps, rank, relativa = TRUE) * 100
  clasificados <- setdiff(rownames(m), "No clasificado")
  top <- clasificados[seq_len(min(top_n, length(clasificados)))]
  grupo <- ifelse(rownames(m) %in% top, rownames(m),
                  ifelse(rownames(m) == "No clasificado", "No clasificado", "Otros"))
  m <- rowsum(m, grupo)

  meta <- metadata_ps(ps)
  df <- data.frame(
    Taxon = rep(rownames(m), times = ncol(m)),
    Muestra = rep(colnames(m), each = nrow(m)),
    Porcentaje = as.vector(m)
  )
  df$Grupo <- meta[df$Muestra, variable_analisis]
  df$Taxon <- factor(df$Taxon, levels = c(top, "Otros", "No clasificado"))

  colores <- c(paleta_pastel(length(top)), "grey85", "grey65")
  names(colores) <- levels(df$Taxon)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Muestra, y = .data$Porcentaje,
                                   fill = .data$Taxon)) +
    ggplot2::geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(~ Grupo, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = colores, drop = TRUE) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = "Abundancia relativa (%)", fill = rank,
                  title = paste0("Composición a nivel de ", rank,
                                 " (top ", length(top), ")"))
}

#' Heatmap de los taxones más abundantes
#'
#' @inheritParams grafico_composicion
#' @param archivo Ruta PNG donde guardar la figura (opcional).
#' @return Objeto `pheatmap` (invisible), o `NULL` si no hay taxones.
#' @export
heatmap_top <- function(ps, rank, variable_analisis, top_n = 20, archivo = NULL) {
  m <- tabla_rango(ps, rank, relativa = TRUE) * 100
  m <- m[rownames(m) != "No clasificado", , drop = FALSE]
  if (nrow(m) == 0) {
    warning("No hay taxones clasificados a nivel de ", rank)
    return(invisible(NULL))
  }
  m <- m[seq_len(min(top_n, nrow(m))), , drop = FALSE]

  meta <- metadata_ps(ps)
  anotacion <- data.frame(meta[colnames(m), variable_analisis, drop = FALSE])
  niveles <- unique(as.character(anotacion[[1]]))
  colores_grupo <- stats::setNames(list(stats::setNames(paleta_pastel(length(niveles)), niveles)),
                                   variable_analisis)

  h <- pheatmap::pheatmap(
    m,
    color = grDevices::colorRampPalette(c("white", paleta_pastel(1)))(100),
    cluster_rows = nrow(m) > 1,
    cluster_cols = ncol(m) > 2,
    annotation_col = anotacion,
    annotation_colors = colores_grupo,
    display_numbers = ncol(m) <= 12,
    number_format = "%.1f",
    fontsize_row = 9,
    fontsize_col = 9,
    main = paste0("Top ", nrow(m), " ", rank, " (% de lecturas)"),
    silent = TRUE
  )

  if (!is.null(archivo)) {
    grDevices::png(archivo, width = 1800, height = 300 + 55 * nrow(m), res = 180)
    grid::grid.draw(h$gtable)
    grDevices::dev.off()
    readr::write_csv(data.frame(Taxon = rownames(m), m, check.names = FALSE),
                     sub("\\.png$", ".csv", archivo))
  }
  invisible(h)
}

#' Diversidad alfa por muestra
#'
#' Riqueza observada, Shannon y Simpson sin rarefacción.
#'
#' @inheritParams grafico_composicion
#' @return data.frame con una fila por muestra.
#' @export
calcular_alfa <- function(ps, variable_analisis) {
  # DADA2 no produce singletons; phyloseq advierte de eso sin que afecte a
  # estas métricas, así que se silencia solo ese aviso.
  alfa <- withCallingHandlers(
    phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson")),
    warning = function(w) {
      if (grepl("singletons", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
  # estimate_richness puede alterar los nombres; el orden sí se conserva.
  alfa <- data.frame(SampleID = phyloseq::sample_names(ps), alfa,
                     Lecturas = phyloseq::sample_sums(ps), row.names = NULL)
  alfa[[variable_analisis]] <- metadata_ps(ps)[alfa$SampleID, variable_analisis]
  alfa
}

#' Gráfico de diversidad alfa
#'
#' @param alfa Resultado de [calcular_alfa()].
#' @param variable_analisis Columna de agrupación.
#' @return Un objeto ggplot.
#' @export
grafico_alfa <- function(alfa, variable_analisis) {
  metricas <- c("Observed", "Shannon", "Simpson")
  df <- data.frame(
    SampleID = rep(alfa$SampleID, length(metricas)),
    Grupo = rep(alfa[[variable_analisis]], length(metricas)),
    Metrica = factor(rep(metricas, each = nrow(alfa)), levels = metricas),
    Valor = unlist(alfa[metricas], use.names = FALSE)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Grupo, y = .data$Valor,
                                   fill = .data$Grupo)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, seed = 1),
                        size = 2.5) +
    ggplot2::facet_wrap(~ Metrica, scales = "free_y") +
    ggplot2::scale_fill_manual(values = paleta_pastel(nlevels(factor(df$Grupo)))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = variable_analisis, y = NULL, title = "Diversidad alfa")
}

# Kruskal-Wallis con p-valor por permutación: la aproximación chi-cuadrado
# no es válida con pocas réplicas (p. ej. 3 vs 3 da 0.0495 cuando el mínimo
# exacto posible es 0.1). Con dos grupos y pocas combinaciones se enumeran
# todas (p exacto); si no, se usan `permutaciones` al azar.
pruebas_alfa <- function(alfa, variable_analisis, permutaciones = 9999) {
  grupo <- droplevels(factor(alfa[[variable_analisis]]))
  if (nlevels(grupo) < 2) return(NULL)
  n <- length(grupo)
  exacto <- nlevels(grupo) == 2 && choose(n, sum(grupo == levels(grupo)[1])) <= 1e5
  etiquetas <- if (exacto) {
    combn(n, sum(grupo == levels(grupo)[1]), simplify = FALSE)
  } else {
    replicate(permutaciones, sample(grupo), simplify = FALSE)
  }
  # Estadístico sin corrección por empates: es monótono con H, basta para comparar.
  h <- function(r, g) sum(tapply(r, g, sum)^2 / tapply(r, g, length))
  as_grupo <- function(e) {
    if (exacto) factor(seq_len(n) %in% e) else e
  }

  do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(m) {
    kw <- stats::kruskal.test(alfa[[m]] ~ grupo)
    r <- rank(alfa[[m]])
    obs <- h(r, grupo)
    perm <- vapply(etiquetas, function(e) h(r, as_grupo(e)), numeric(1))
    p <- if (exacto) mean(perm >= obs - 1e-9) else (sum(perm >= obs - 1e-9) + 1) / (length(perm) + 1)
    data.frame(Metrica = m, H = unname(kw$statistic), gl = unname(kw$parameter),
               p_valor = p, metodo = if (exacto) "exacto" else "permutaciones")
  }))
}

#' Diversidad beta: Bray-Curtis, PCoA, PERMANOVA y dispersión
#'
#' @inheritParams grafico_composicion
#' @param permutations Número de permutaciones.
#' @return Lista con `plot`, `coordinates`, `permanova` y `dispersion`.
#' @export
analisis_beta <- function(ps, variable_analisis, permutations = 999) {
  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  dist_bc <- phyloseq::distance(ps_rel, method = "bray")
  ord <- phyloseq::ordinate(ps_rel, method = "PCoA", distance = dist_bc)

  eig <- ord$values$Eigenvalues
  pct <- 100 * eig[1:2] / sum(eig[eig > 0])

  meta <- metadata_ps(ps)
  ord_df <- data.frame(SampleID = rownames(ord$vectors),
                       PCoA1 = ord$vectors[, 1], PCoA2 = ord$vectors[, 2])
  ord_df[[variable_analisis]] <- meta[ord_df$SampleID, variable_analisis]

  grupo <- droplevels(factor(meta[labels(dist_bc), variable_analisis]))
  n_grupo <- table(grupo)

  # Relleno pastel con borde oscuro: sin borde, los puntos casi no se ven sobre blanco.
  colores <- paleta_pastel(nlevels(grupo))
  p <- ggplot2::ggplot(ord_df, ggplot2::aes(x = .data$PCoA1, y = .data$PCoA2,
                                            fill = .data[[variable_analisis]])) +
    ggplot2::geom_point(shape = 21, size = 4, colour = "grey30") +
    ggplot2::scale_fill_manual(values = colores) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "PCoA de Bray-Curtis",
                  x = sprintf("PCoA1 (%.1f%%)", pct[1]),
                  y = sprintf("PCoA2 (%.1f%%)", pct[2]),
                  fill = variable_analisis)
  # Las elipses requieren al menos 4 muestras por grupo.
  if (length(n_grupo) > 0 && all(n_grupo >= 4)) {
    p <- p +
      ggplot2::stat_ellipse(ggplot2::aes(colour = .data[[variable_analisis]]),
                            type = "t", level = 0.95, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = colores)
  }

  permanova <- NULL
  dispersion <- NULL
  keep <- !is.na(grupo)
  if (nlevels(grupo) >= 2 && sum(keep) >= 3) {
    d <- stats::as.dist(as.matrix(dist_bc)[keep, keep])
    g <- droplevels(grupo[keep])
    permanova <- as.data.frame(vegan::adonis2(d ~ g, permutations = permutations))
    rownames(permanova)[1] <- variable_analisis
    bd <- vegan::betadisper(d, g)
    dispersion <- as.data.frame(vegan::permutest(bd, permutations = permutations)$tab)
  }

  list(distance = dist_bc, coordinates = ord_df, plot = p,
       permanova = permanova, dispersion = dispersion)
}
