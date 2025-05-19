#' Build binary presence/absence matrix of sporulation genes
#'
#' Transforms the output of `sporulation_gene_name()` into a wide-format matrix
#' indicating the presence (1) or absence (0) of each sporulation-associated gene per genome.
#'
#' @param df A data.frame from `sporulation_gene_name()` with columns `codigo_genoma` and `consensus_name_this_study`.
#'
#' @return A wide-format binary matrix with genomes in rows and genes in columns.
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @export
build_binary_matrix <- function(df) {
  # Lista de todos os genes obrigatórios
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

  # Criar matriz binária inicial
  df_bin <- df %>%
    mutate(presenca = 1) %>%
    select(codigo_genoma, consensus_name_this_study, presenca) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = consensus_name_this_study,
      values_from = presenca,
      values_fill = list(presenca = 0)
    )

  # Corrigir nomes inválidos
  colnames(df_bin) <- gsub("[^a-zA-Z0-9_]", "_", colnames(df_bin))

  # Garantir presença de todos os genes
  missing_genes <- setdiff(required_genes, colnames(df_bin))
  for (gene in missing_genes) {
    df_bin[[gene]] <- 0
  }

  # Reordenar colunas (mantendo codigo_genoma primeiro)
  df_bin <- df_bin[, c("codigo_genoma", sort(setdiff(colnames(df_bin), "codigo_genoma")))]

  return(df_bin)
}
