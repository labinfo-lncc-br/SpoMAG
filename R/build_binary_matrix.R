#' Build binary presence/absence matrix of sporulation genes
#'
#' Transforms the output of `sporulation_gene_name()` into a wide-format matrix
#' indicating the presence (1) or absence (0) of each sporulation-associated gene per genome.
#'
#' @param df A data.frame from `sporulation_gene_name()` with columns `genome_ID` and `spo_gene_name`.
#'
#' @return A wide-format binary matrix with genomes in rows and genes in columns.
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @export
build_binary_matrix <- function(df) {
  required_genes <- c(
    "spo0A", "sigH", "spoIIE", "spoIIIE", "spoIIIJ", "pth", "spoVG", "spoVS", "divIC",
    "divIB", "divIVA", "ftsA", "ftsE", "ftsH", "ftsX", "ftsY", "ftsZ", "jag", "minC",
    "minD", "spo0B", "spo0F", "ald", "obg", "ftsL", "ymcA", "ylbF", "yaaT", "sda",
    "sigE", "sigF", "sigG", "spoIIAA", "spoIIAB", "spoIIGA", "parA", "soj", "parB", "spoIID",
    "spoIIM", "spoIIP", "spoIIQ", "spoIIIAA", "spoIIIAB", "spoIIIAC", "spoiiiAD", "spoIIIAE", "spoIIIAF", "spoIIIAG",
    "spoIIIAH", "spoIIB", "yunB", "sigK", "spoIIID", "spoIVFB", "spoVB", "spoVK", "ftsW", "bofA",
    "alr", "dacB", "spmA", "spmB", "yisY", "ylmC", "ytaF", "ytvI", "ctpB", "spoIVFA",
    "ykvI", "yqfU", "ydca", "ydcc", "yhbh", "spoIIR", "spoIVB", "spoVT", "dacF", "ytfJ",
    "yloC", "rsfA", "fin_yabk", "ymfJ", "yqhG", "ywzB", "spoVAC", "spoVAD", "nfo", "mglA_sspA",
    "sspB", "spoVAA", "spoVAF", "sspH", "tlp", "spoVAEB", "spoVAB", "sspC", "sspI", "sspF",
    "spoVFA", "spoVFB", "spoIVA", "cotJC", "cotSA", "gerM", "safA", "ydhD", "yhaX", "cotJB",
    "yhjR", "spoVID", "cotE", "hsps", "spoVD", "ylbJ", "cwlC_cwlD", "lytH", "yabP", "yabQ",
    "yqfC", "yqfD", "gerA", "gpr", "CSD", "gdh", "ypeB", "gerC", "lgt", "gerD",
    "gerE"
  )

  df_bin <- df %>%
    mutate(present = 1) %>%
    select(genome_ID, spo_gene_name, present) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = spo_gene_name,
      values_from = present,
      values_fill = list(present = 0)
    )

  colnames(df_bin) <- gsub("[^a-zA-Z0-9_]", "_", colnames(df_bin))

  missing_genes <- setdiff(required_genes, colnames(df_bin))
  for (gene in missing_genes) {
    df_bin[[gene]] <- 0
  }

  df_bin <- df_bin[, c("genome_ID", sort(setdiff(colnames(df_bin), "genome_ID")))]

  return(df_bin)
}
