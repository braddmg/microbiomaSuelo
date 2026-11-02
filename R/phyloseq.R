construir_phyloseq <- function(seqtab_nochim, tax, metadata,
                               columna_id = "SampleID",
                               tipo_datos = "16S",
                               min_lecturas_asv = 5) {
  ids_seq <- rownames(seqtab_nochim)
  ids_meta <- as.character(metadata[[columna_id]])

  faltan <- setdiff(ids_seq, ids_meta)
  if (length(faltan) > 0) {
    stop("Estas muestras del FASTQ no están en metadata: ",
         paste(faltan, collapse = ", "))
  }
  sin_fastq <- setdiff(ids_meta, ids_seq)
  if (length(sin_fastq) > 0) {
    warning("Muestras de la metadata sin lecturas tras DADA2: ",
            paste(sin_fastq, collapse = ", "))
  }

  metadata <- metadata[match(ids_seq, ids_meta), , drop = FALSE]
  rownames(metadata) <- metadata[[columna_id]]

  seqs <- colnames(seqtab_nochim)
  asv_ids <- paste0("ASV", seq_along(seqs))
  colnames(seqtab_nochim) <- asv_ids
  taxmat <- as.matrix(tax[seqs, RANGOS, drop = FALSE])
  rownames(taxmat) <- asv_ids
  refseq <- Biostrings::DNAStringSet(seqs)
  names(refseq) <- asv_ids

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(t(seqtab_nochim), taxa_are_rows = TRUE),
    phyloseq::tax_table(taxmat),
    phyloseq::sample_data(metadata),
    refseq
  )

  # ASVs sin reino asignado (y en 16S, cloroplastos, mitocondrias y
  # eucariotas) se eliminan antes de los análisis.
  taxdf <- tabla_taxonomia(ps)
  keep <- !is.na(taxdf$Kingdom)
  if (tipo_datos == "16S") {
    keep <- keep &
      taxdf$Kingdom != "Eukaryota" &
      (is.na(taxdf$Order) | taxdf$Order != "Chloroplast") &
      (is.na(taxdf$Family) | taxdf$Family != "Mitochondria")
  }
  ps <- phyloseq::prune_taxa(keep, ps)
  n_previo <- phyloseq::ntaxa(ps)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) >= max(min_lecturas_asv, 1), ps)
  message("ASVs con menos de ", min_lecturas_asv, " lecturas en total eliminadas: ",
          n_previo - phyloseq::ntaxa(ps), " de ", n_previo, ".")
  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 0, ps)

  if (phyloseq::nsamples(ps) < 2) {
    stop("Quedaron menos de dos muestras con lecturas utilizables.")
  }

  ps
}

#' Tabla de conteos agregada por rango taxonómico
#'
#' Las ASVs sin clasificar en `rank` se agrupan como "No clasificado", de
#' modo que las abundancias relativas se calculan sobre el total de lecturas.
#'
#' @param ps Objeto phyloseq.
#' @param rank Rango taxonómico (p. ej. "Phylum", "Genus").
#' @param relativa Si `TRUE` devuelve proporciones por muestra.
#' @return Matriz taxones x muestras.
#' @export
tabla_rango <- function(ps, rank, relativa = FALSE) {
  tax <- tabla_taxonomia(ps)
  etiqueta <- etiqueta_taxon(tax, rank)
  etiqueta[is.na(etiqueta)] <- "No clasificado"
  m <- rowsum(matriz_otu(ps), etiqueta)
  if (relativa) m <- sweep(m, 2, pmax(colSums(m), 1), "/")
  m[order(rowSums(m), decreasing = TRUE), , drop = FALSE]
}

exportar_tablas <- function(ps, carpeta_tablas, carpeta_abund) {
  otu <- matriz_otu(ps)
  tax <- tabla_taxonomia(ps)
  base <- data.frame(ASV = rownames(tax), tax,
                     Secuencia = as.character(phyloseq::refseq(ps)),
                     check.names = FALSE)

  rel <- sweep(otu, 2, colSums(otu), "/")
  readr::write_csv(cbind(base, otu), file.path(carpeta_tablas, "ASV_counts_taxonomy.csv"))
  readr::write_csv(cbind(base, rel), file.path(carpeta_tablas, "ASV_relative_abundance.csv"))
  Biostrings::writeXStringSet(phyloseq::refseq(ps), file.path(carpeta_tablas, "ASV_secuencias.fasta"))

  for (r in RANGOS[-1]) {
    m <- tabla_rango(ps, r, relativa = TRUE)
    readr::write_csv(
      data.frame(Taxon = rownames(m), m * 100, check.names = FALSE),
      file.path(carpeta_abund, paste0("abundancia_relativa_", tolower(r), "_porcentaje.csv"))
    )
  }
  invisible(NULL)
}
