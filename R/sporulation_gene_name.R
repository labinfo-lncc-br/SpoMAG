#' Identify Sporulation-Associated Genes
#'
#' This function identifies sporulation-associated genes in a genome annotation data frame.
#' It searches for gene names and KEGG Orthology identifiers related to sporulation steps and returns
#' a data frame with annotated sporulation genes and a consensus name.
#'
#' @param df A data frame containing MAG annotation with the columns 'Preferred_name', 'KEGG_ko', and 'genome_ID'.
#'
#' @return A data frame of sporulation-associated genes with standardized names and spo_process tags.
#' @import dplyr
#' @importFrom tidyr pivot_wider
#'
#' @export
sporulation_gene_name <- function(df){
  ##### Sporulation onset and checkpoints #####
  message("Looking for genes related to Sporulation Onset and Checkpoints")

  ### Spo0A ###
  spo0a_gene <- df %>% filter(grepl('spo0a', tolower(Preferred_name)))
  tryCatch({
    spo0a_ko <- df %>% filter(grepl('07699', tolower(KEGG_ko)))
    spo0a_ko <- spo0a_ko[spo0a_ko$Preferred_name == "-",]
    spo0a_gene <- rbind(spo0a_gene, spo0a_ko)
  }, error = function(e) {})
  tryCatch({
    spo0a_gene$spo_gene_name <- "spo0A"
  }, error = function(e) {})
  ### *** spo0a_gene


  tryCatch({
    ### sigH ###
    sigh_gene <- df %>% filter(grepl('sigh', tolower(Preferred_name)))
    sigh_gene$spo_gene_name <- "sigH"
    ### *** sigh_gene
  }, error = function(e) {})

  ### spoIIE ###
  spoiie_gene <- df %>% filter(grepl('spoiie', tolower(Preferred_name)))
  tryCatch({
    spoiie_ko <- df %>% filter(grepl('06382', tolower(KEGG_ko)))
    spoiie_ko <- spoiie_ko[spoiie_ko$Preferred_name == "-",]
    spoiie_gene <- rbind(spoiie_gene, spoiie_ko)
  }, error = function(e) {})
  tryCatch({
    spoiie_gene$spo_gene_name <- "spoIIE"
  }, error = function(e) {})
  ### *** spoiie_gene



  ### spoIIIE/ftsk ###
  spoiiie_gene <- df %>% filter(grepl('spoiiie', tolower(Preferred_name)))
  tryCatch({
    spoiiie_ko <- df %>% filter(grepl('03466', tolower(KEGG_ko)))
    spoiiie_ko_hifen <- spoiiie_ko[spoiiie_ko$Preferred_name == "-",]
    spoiiie_ko_essc <- spoiiie_ko[spoiiie_ko$Preferred_name == "essC",]
    spoiiie_ko_ftsk <- spoiiie_ko[spoiiie_ko$Preferred_name == "ftsK",]
    spoiiie_ko_hifen <- spoiiie_ko_hifen[spoiiie_ko_hifen$EC == "-",]
    spoiiie_gene <- rbind(spoiiie_gene, spoiiie_ko_hifen)
    spoiiie_gene <- rbind(spoiiie_gene, spoiiie_ko_essc)
    spoiiie_gene <- rbind(spoiiie_gene, spoiiie_ko_eccca)
    spoiiie_gene <- rbind(spoiiie_gene, spoiiie_ko_ftsk)
    ### *** spoiiie_gene
  }, error = function(e) {})
  tryCatch({
    spoiiie_gene$spo_gene_name <- "spoIIIE"
  }, error = function(e) {})


  ### spoIIIJ/yidC ###
  yidc_gene <- df %>% filter(grepl('yidc', tolower(Preferred_name)))
  tryCatch({
    yidc_ko <- df %>% filter(grepl('03217', tolower(KEGG_ko)))
    yidc_ko <- yidc_ko[yidc_ko$Preferred_name == "-",]
    yidc_gene <- rbind(yidc_gene, yidc_ko)
  }, error = function(e) {})
  tryCatch({
    yidc_gene$spo_gene_name <- "spoIIIJ"
  }, error = function(e) {})
  ### *** yidc_gene



  ### pth ###
  pth_gene <- df %>% filter(grepl('pth', tolower(Preferred_name)))
  tryCatch({
    pth_ko <- df %>% filter(grepl('01056', tolower(KEGG_ko)))
    pth_ko <- pth_ko[pth_ko$Preferred_name == "-",]
    pth_gene <- rbind(pth_gene, pth_ko)
  }, error = function(e) {})
  tryCatch({
    pth_gene$spo_gene_name <- "pth"
  }, error = function(e) {})
  ### *** pth_gene


  ### spoVG ###
  spovg_gene <- df %>% filter(grepl('spovg', tolower(Preferred_name)))
  tryCatch({
    spovg_ko <- df %>% filter(grepl('06412', tolower(KEGG_ko)))
    spovg_com_nome <- spovg_ko[spovg_ko$Preferred_name == "-",]
    spovg_gene <- rbind(spovg_gene, spovg_com_nome)
  }, error = function(e) {})
  tryCatch({
    spovg_gene$spo_gene_name <- "spoVG"
  }, error = function(e) {})
  ### *** spovg_gene

  ### spoVS ###
  spovs_gene <- df %>% filter(grepl('spovs', tolower(Preferred_name)))
  tryCatch({
    spovs_ko <- df %>% filter(grepl('06416', tolower(KEGG_ko)))
    spovs_ko <- spovs_ko[spovs_ko$Preferred_name == "-",]
    spovs_gene <- rbind(spovs_gene, spovs_ko)
  }, error = function(e) {})
  tryCatch({
    spovs_gene$spo_gene_name <- "spoVS"
  }, error = function(e) {})
  ### *** spovs_gene

  ### divIB ###
  divib_gene <- df %>% filter(grepl('divib', tolower(Preferred_name)))
  tryCatch({
    divib_ko <- df %>% filter(grepl('03589', tolower(KEGG_ko)))
    divib_ko_hifen <- divib_ko[divib_ko$Preferred_name == "-",]
    divib_ko_ftsq <- divib_ko[divib_ko$Preferred_name == "ftsQ",]
    divib_gene <- rbind(divib_gene, divib_ko_hifen)
    divib_gene <- rbind(divib_gene, divib_ko_ftsq)
  }, error = function(e) {})
  tryCatch({
    divib_gene$spo_gene_name <- "divIB"
  }, error = function(e) {})
  ### *** divib_gene

  ### divIC ###
  divic_gene <- df %>% filter(grepl('divic', tolower(Preferred_name)))
  tryCatch({
    divic_ko <- df %>% filter(grepl('13052', tolower(KEGG_ko)))
    divic_ko <- divic_ko[divic_ko$Preferred_name == "-",]
    divic_gene <- rbind(divic_gene, divic_ko)
  }, error = function(e) {})
  tryCatch({
    divic_gene$spo_gene_name <- "divIC"
  }, error = function(e) {})
  ### *** divic_gene

  ### divIVA ###
  diviva_gene <- df %>% filter(grepl('diviva', tolower(Preferred_name)))
  tryCatch({
    diviva_ko <- df %>% filter(grepl('04074', tolower(KEGG_ko)))
    diviva_ko <- diviva_ko[diviva_ko$Preferred_name == "-",]
    diviva_gene <- rbind(diviva_gene, diviva_ko)
  }, error = function(e) {})
  tryCatch({
    diviva_gene$spo_gene_name <- "divIVA"
  }, error = function(e) {})
  ### *** diviva_gene

  tryCatch({
    ### ftsA ###
    ftsa_gene <- df %>% filter(grepl('ftsa', tolower(Preferred_name)))
    ftsa_gene$spo_gene_name <- "ftsA"
    ### *** ftsa_gene
  }, error = function(e) {})

  ### ftsE ###
  ftse_gene <- df %>% filter(grepl('ftse', tolower(Preferred_name)))
  tryCatch({
    ftse_ko <- df %>% filter(grepl('09812', tolower(KEGG_ko)))
    ftse_ko <- ftse_ko[ftse_ko$Preferred_name == "-",]
    ftse_gene <- rbind(ftse_gene, ftse_ko)
  }, error = function(e) {})
  tryCatch({
    ftse_gene$spo_gene_name <- "ftsE"
  }, error = function(e) {})
  ### *** ftse_gene

  ### ftsH ###
  ftsh_gene <- df %>% filter(grepl('ftsh', tolower(Preferred_name)))
  tryCatch({
    ftsh_ko <- df %>% filter(grepl('03798', tolower(KEGG_ko)))
    ftsh_ko <- ftsh_ko[ftsh_ko$Preferred_name == "-",]
    ftsh_gene <- rbind(ftsh_gene, ftsh_ko)
  }, error = function(e) {})
  tryCatch({
    ftsh_gene$spo_gene_name <- "ftsH"
  }, error = function(e) {})
  ### *** ftsh_gene

  ### ftsX ###
  ftsX_gene <- df %>% filter(grepl('ftsx', tolower(Preferred_name)))
  tryCatch({
    ftsX_ko <- df %>% filter(grepl('09811', tolower(KEGG_ko)))
    ftsX_ko$spo_gene_name <- "ftsX"
    ftsX_gene <- ftsX_ko
  }, error = function(e) {})
  tryCatch({
    ftsX_gene$spo_gene_name <- "ftsX"
  }, error = function(e) {})
  ### *** ftsX_gene

  tryCatch({
    ### ftsY ###
    ftsy_gene <- df %>% filter(grepl('ftsy', tolower(Preferred_name)))
    ftsy_gene$spo_gene_name <- "ftsY"
    ### *** ftsy_gene
  }, error = function(e) {})

  tryCatch({
    ### ftsZ ###
    ftsz_gene <- df %>% filter(grepl('ftsz', tolower(Preferred_name)))
    ftsz_gene$spo_gene_name <- "ftsZ"
    ### *** ftsz_gene
  }, error = function(e) {})

  ### jag ###
  jag_gene <- df %>% filter(grepl('jag', tolower(Preferred_name)))
  tryCatch({
    jag_ko <- df %>% filter(grepl('06346', tolower(KEGG_ko)))
    jag_ko <- jag_ko[jag_ko$Preferred_name == "-",]
    jag_gene <- rbind(jag_gene, jag_ko)
  }, error = function(e) {})
  ### *** jag_gene
  tryCatch({
    jag_gene$spo_gene_name <- "jag"
  }, error = function(e) {})

  ### minC ###
  minc_gene <- df %>% filter(grepl('minc', tolower(Preferred_name)))
  tryCatch({
    minc_ko <- df %>% filter(grepl('03610', tolower(KEGG_ko)))
    minc_gene <-minc_ko
  }, error = function(e) {})
  tryCatch({
    minc_gene$spo_gene_name <- "minC"
  }, error = function(e) {})
  ### *** minc_gene

  tryCatch({
    ### mind ###
    mind_gene <- df %>% filter(grepl('mind', tolower(Preferred_name)))
    mind_gene$spo_gene_name <- "minD"
    ### *** mind_gene
  }, error = function(e) {})

  tryCatch({
    ### obg ###
    obg_gene <- df %>% filter(grepl('obg', tolower(Preferred_name)))
    obg_gene$spo_gene_name <- "obg"
    ### *** obg_gene
  }, error = function(e) {})

  ### spo0B ###
  spo0b_gene <- df %>% filter(grepl('spo0b', tolower(Preferred_name)))
  tryCatch({
    spo0b_ko <- df %>% filter(grepl('06375', tolower(KEGG_ko)))
    spo0b_gene <- spo0b_ko
  }, error = function(e) {})
  tryCatch({
    spo0b_gene$spo_gene_name <- "spo0B"
  }, error = function(e) {})
  ### *** spo0b_gene

  ### spo0F ###
  spo0f_gene <- df %>% filter(grepl('spo0f', tolower(Preferred_name)))
  tryCatch({
    spo0f_ko <- df %>% filter(grepl('02490', tolower(KEGG_ko)))
    spo0f_ko <- spo0f_ko[spo0f_ko$Preferred_name == "-",]
    spo0f_gene <- rbind(spo0f_gene, spo0f_ko)
  }, error = function(e) {})
  tryCatch({
    spo0f_gene$spo_gene_name <- "spo0F"
  }, error = function(e) {})
  ### *** spo0f_gene

  ### ald ###
  ald_gene <- df %>% filter(grepl('alda', tolower(Preferred_name)))
  tryCatch({
    ald_ko <- df %>% filter(grepl('00259', tolower(KEGG_ko)))
    ald_ko <- ald_ko[ald_ko$Preferred_name == "-",]
    ald_gene <- rbind(ald_gene, ald_ko)
  }, error = function(e) {})
  tryCatch({
    ald_gene$spo_gene_name <- "ald"
  }, error = function(e) {})
  ### *** ald_gene

  tryCatch({
    ### ftsL ###
    ftsl_gene <- df %>% filter(grepl('ftsl', tolower(Preferred_name)))
    ftsl_gene$spo_gene_name <- "ftsL"
    ### *** ftsl_gene
  }, error = function(e) {})

  tryCatch({
    ### ymcA ###
    ymca_gene <- df %>% filter(grepl('ymca', tolower(Preferred_name)))
    ymca_gene$spo_gene_name <- "ymcA"
    ### *** ymca_gene
  }, error = function(e) {})

  tryCatch({
    ### ylbF ###
    ylbf_gene <- df %>% filter(grepl('ylbf', tolower(Preferred_name)))
    ylbf_gene$spo_gene_name <- "ylbF"
    ### *** ylbf_gene
  }, error = function(e) {})

  tryCatch({
    ### yaaT ###
    yaat_gene <- df %>% filter(grepl('yaat', tolower(Preferred_name)))
    yaat_gene$spo_gene_name <- "yaaT"
    ### *** yaat_gene
  }, error = function(e) {})

  tryCatch({
    ### sda ###
    sda_gene <- df[df$Preferred_name == 'sda',]
    sda_gene$spo_gene_name <- "sda"
    ### *** sda_gene
  }, error = function(e) {})

  tryCatch({
    sporulation_onset_and_checkpoints <- rbind(spo0a_gene, sigh_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spoiie_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spoiiie_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, yidc_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, pth_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spovg_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spovs_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, divic_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, divib_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, diviva_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsa_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftse_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsh_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsX_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsy_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsz_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, jag_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, minc_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, mind_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spo0b_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, spo0f_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ald_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, obg_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ftsl_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ymca_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, ylbf_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, yaat_gene)
    sporulation_onset_and_checkpoints <- rbind(sporulation_onset_and_checkpoints, sda_gene)
    sporulation_onset_and_checkpoints$spo_process <- "sporulation_onset_and_checkpoints"
  }, error = function(e) {})

  ##### Spo0A regulon #####
  message("Looking for genes related to Spo0A Regulon")

  tryCatch({
    ### sigE ###
    sigE_gene <- df %>% filter(grepl('sige', tolower(Preferred_name)))
    sigE_gene$spo_gene_name <- "sigE"
    ### *** sigE_gene
  }, error = function(e) {})

  tryCatch({
    ### sigF ###
    sigF_gene <- df %>% filter(grepl('sigf', tolower(Preferred_name)))
    sigF_gene$spo_gene_name <- "sigF"
    ### *** sigF_gene
  }, error = function(e) {})

  tryCatch({
    ### sigG ###
    sigg_gene <- df %>% filter(grepl('sigg', tolower(Preferred_name)))
    sigg_gene$spo_gene_name <- "sigG"
    ### *** sigg_gene
  }, error = function(e) {})

  tryCatch({
    ### spoiiaa ###
    spoiiaa_gene <- df %>% filter(grepl('spoiiaa', tolower(Preferred_name)))
    spoiiaa_gene$spo_gene_name <- "spoIIAA"
    ### *** spoiiaa_gene
  }, error = function(e) {})

  tryCatch({
    ### spoiiab ###
    spoiiab_gene <- df %>% filter(grepl('spoiiab', tolower(Preferred_name)))
    spoiiab_gene$spo_gene_name <- "spoIIAB"
    ### *** spoiiab_gene
  }, error = function(e) {})

  ### spoiiga ###
  spoiiga_gene <- df %>% filter(grepl('spoiiga', tolower(Preferred_name)))
  tryCatch({
    spoiiga_ko <- df %>% filter(grepl('06383', tolower(KEGG_ko)))
    spoiiga_ko <- spoiiga_ko[spoiiga_ko$Preferred_name == "-",]
    spoiiga_gene <- rbind(spoiiga_gene, spoiiga_ko)
  }, error = function(e) {})
  tryCatch({
    spoiiga_gene$spo_gene_name <- "spoIIGA"
  }, error = function(e) {})
  ### *** spoiiga_gene

  tryCatch({
    ### parA ###
    para_gene <- df %>% filter(grepl('para', tolower(Preferred_name)))
    para_gene$spo_gene_name <- "parA"
    ### *** para_gene
  }, error = function(e) {})

  tryCatch({
    ### soj ###
    soj_gene <- df %>% filter(grepl('soj', tolower(Preferred_name)))
    soj_gene$spo_gene_name <- "soj"
    ### *** soj_gene
  }, error = function(e) {})

  ### parB ###
  parb_gene <- df %>% filter(grepl('parb', tolower(Preferred_name)))
  tryCatch({
    parb_ko <- df %>% filter(grepl('03497', tolower(KEGG_ko)))
    parb_gene <- parb_ko
    ### *** parb_gene
  }, error = function(e) {})
  tryCatch({
    parb_gene$spo_gene_name <- "parB"
  }, error = function(e) {})

  tryCatch({
    spo0a_regulon <- rbind(sigE_gene, sigF_gene)
    spo0a_regulon <- rbind(spo0a_regulon, sigg_gene)
    spo0a_regulon <- rbind(spo0a_regulon, spoiiaa_gene)
    spo0a_regulon <- rbind(spo0a_regulon, spoiiab_gene)
    spo0a_regulon <- rbind(spo0a_regulon, spoiiga_gene)
    spo0a_regulon <- rbind(spo0a_regulon, para_gene)
    spo0a_regulon <- rbind(spo0a_regulon, soj_gene)
    spo0a_regulon <- rbind(spo0a_regulon, parb_gene)
    spo0a_regulon$spo_process <- "spo0a_regulon"
  }, error = function(e) {})

  ##### Engulfment #####
  message("Looking for genes related to Engulfment")

  ### spoIID ###
  spoiid_gene <- df %>% filter(grepl('spoiid', tolower(Preferred_name)))
  tryCatch({
    spoiid_ko <- df %>% filter(grepl('06381', tolower(KEGG_ko)))
    spoiid_gene <- spoiid_ko
  }, error = function(e) {})
  tryCatch({
    spoiid_gene$spo_gene_name <- "spoIID"
  }, error = function(e) {})
  ### *** spoiid_gene

  ### spoiim ###
  spoiim_gene <- df %>% filter(grepl('spoiim', tolower(Preferred_name)))
  tryCatch({
    spoiim_ko <- df %>% filter(grepl('06384', tolower(KEGG_ko)))
    spoiim_ko <- spoiim_ko[spoiim_ko$Description != "-",]
    spoiim_gene <- spoiim_ko
  }, error = function(e) {})
  tryCatch({
    spoiim_gene$spo_gene_name <- "spoIIM"
  }, error = function(e) {})
  ### *** spoiim_gene

  tryCatch({
    ### spoiip ###
    spoiip_gene <- df %>% filter(grepl('spoiip', tolower(Preferred_name)))
    spoiip_gene$spo_gene_name <- "spoIIP"
    ### *** spoiip_gene
  }, error = function(e) {})

  ### spoiiq ###
  spoiiq_gene <- df %>% filter(grepl('spoiiq', tolower(Preferred_name)))
  tryCatch({
    spoiiq_ko <- df %>% filter(grepl('06386', tolower(KEGG_ko)))
    spoiiq_gene <- spoiiq_ko
  }, error = function(e) {})
  tryCatch({
    spoiiq_gene$spo_gene_name <- "spoIIQ"
  }, error = function(e) {})
  ### *** spoiiq_gene

  ### spoiiiaa ###
  spoiiiaa_gene <- df %>% filter(grepl('spoiiiaa', tolower(Preferred_name)))
  tryCatch({
    spoiiiaa_ko <- df %>% filter(grepl('06390', tolower(KEGG_ko)))
    spoiiiaa_ko <- spoiiiaa_ko[spoiiiaa_ko$Preferred_name == "-",]
    spoiiiaa_gene <- rbind(spoiiiaa_gene,spoiiiaa_ko)
  }, error = function(e) {})
  tryCatch({
    spoiiiaa_gene$spo_gene_name <- "spoIIIAA"
  }, error = function(e) {})
  ### *** spoiiiaa_gene

  ### spoiiiab ###
  spoiiiab_gene <- df %>% filter(grepl('spoiiiab', tolower(Preferred_name)))
  tryCatch({
    spoiiiab_ko <- df %>% filter(grepl('06391', tolower(KEGG_ko)))
    spoiiiab_gene <- spoiiiab_ko
  }, error = function(e) {})
  tryCatch({
    spoiiiab_gene$spo_gene_name <- "spoIIIAB"
  }, error = function(e) {})
  ### *** spoiiiab_gene

  ### spoiiiac ###
  spoiiiac_gene <- df %>% filter(grepl('spoiiiac', tolower(Preferred_name)))
  tryCatch({
    spoiiiac_ko <- df %>% filter(grepl('06392', tolower(KEGG_ko)))
    spoiiiac_ko <- spoiiiac_ko[spoiiiac_ko$Preferred_name == "-",]
    spoiiiac_gene <- rbind(spoiiiac_gene,spoiiiac_ko)
  }, error = function(e) {})
  tryCatch({
    spoiiiac_gene$spo_gene_name <- "spoIIIAC"
  }, error = function(e) {})
  ### *** spoiiiac_gene

  ### spoiiiad ###
  spoiiiad_gene <- df %>% filter(grepl('spoiiiad', tolower(Preferred_name)))
  tryCatch({
    spoiiiad_ko <- df %>% filter(grepl('06393', tolower(KEGG_ko)))
    spoiiiad_gene <- spoiiiad_ko
  }, error = function(e) {})
  tryCatch({
    spoiiiad_gene$spo_gene_name <- "spoiiiAD"
  }, error = function(e) {})
  ### *** spoiiiad_gene

  ### spoiiiae ###
  spoiiiae_gene <- df %>% filter(grepl('spoiiiae', tolower(Preferred_name)))
  tryCatch({
    spoiiiae_ko <- df %>% filter(grepl('06394', tolower(KEGG_ko)))
    spoiiiae_ko <- spoiiiae_ko[spoiiiae_ko$Preferred_name == "-",]
    spoiiiae_gene <- rbind(spoiiiae_gene,spoiiiae_ko)
  }, error = function(e) {})
  tryCatch({
    spoiiiae_gene$spo_gene_name <- "spoIIIAE"
  }, error = function(e) {})
  ### *** spoiiiae_gene

  ### spoiiiaf ###
  spoiiiaf_gene <- df %>% filter(grepl('spoiiiaf', tolower(Preferred_name)))
  tryCatch({
    spoiiiaf_ko <- df %>% filter(grepl('06395', tolower(KEGG_ko)))
    spoiiiaf_gene <- spoiiiaf_ko
  }, error = function(e) {})
  tryCatch({
    spoiiiaf_gene$spo_gene_name <- "spoIIIAF"
  }, error = function(e) {})
  ### *** spoiiiaf_gene

  ### spoiiiag ###
  spoiiiag_gene <- df %>% filter(grepl('spoiiiag', tolower(Preferred_name)))
  tryCatch({
    spoiiiag_ko <- df %>% filter(grepl('06396', tolower(KEGG_ko)))
    spoiiiag_ko <- spoiiiag_ko[spoiiiag_ko$Preferred_name == "-",]
    spoiiiag_gene <- rbind(spoiiiag_gene,spoiiiag_ko)
  }, error = function(e) {})
  tryCatch({
    spoiiiag_gene$spo_gene_name <- "spoIIIAG"
  }, error = function(e) {})
  ### *** spoiiiag_gene

  ### spoiiiah ###
  spoiiiah_gene <- df %>% filter(grepl('spoiiiah', tolower(Preferred_name)))
  tryCatch({
    spoiiiah_ko <- df %>% filter(grepl('06397', tolower(KEGG_ko)))
    spoiiiah_gene <- spoiiiah_ko
  }, error = function(e) {})
  tryCatch({
    spoiiiah_gene$spo_gene_name <- "spoIIIAH"
  }, error = function(e) {})
  ### *** spoiiiah_gene

  ### spoIIB ###
  spoiib_gene <- df %>% filter(grepl('spoiib', tolower(Preferred_name)))
  tryCatch({
    spoiib_ko <- df %>% filter(grepl('06380', tolower(KEGG_ko)))
    spoiib_gene <- spoiib_ko
  }, error = function(e) {})
  tryCatch({
    spoiib_gene$spo_gene_name <- "spoIIB"
  }, error = function(e) {})
  ### *** spoiib_gene

  tryCatch({
    ### yunb ###
    yunb_gene <- df %>% filter(grepl('yunb', tolower(Preferred_name)))
    yunb_gene$spo_gene_name <- "yunB"
    ### *** yunb_gene
  }, error = function(e) {})

  tryCatch({
    Engulfment <- rbind(spoiid_gene, spoiim_gene)
    Engulfment <- rbind(Engulfment, spoiip_gene)
    Engulfment <- rbind(Engulfment, spoiiq_gene)
    Engulfment <- rbind(Engulfment, spoiiiaa_gene)
    Engulfment <- rbind(Engulfment, spoiiiab_gene)
    Engulfment <- rbind(Engulfment, spoiiiac_gene)
    Engulfment <- rbind(Engulfment, spoiiiad_gene)
    Engulfment <- rbind(Engulfment, spoiiiae_gene)
    Engulfment <- rbind(Engulfment, spoiiiaf_gene)
    Engulfment <- rbind(Engulfment, spoiiiag_gene)
    Engulfment <- rbind(Engulfment, spoiiiah_gene)
    Engulfment <- rbind(Engulfment, spoiib_gene)
    Engulfment <- rbind(Engulfment, yunb_gene)
    Engulfment$spo_process <- "Engulfment"
  }, error = function(e) {})

  ##### sigF regulon #####
  message("Looking for genes related to SigF Regulon")


  ### spoiir ###
  spoiir_gene <- df %>% filter(grepl('spoiir', tolower(Preferred_name)))
  tryCatch({
    spoiir_ko <- df %>% filter(grepl('06387', tolower(KEGG_ko)))
    spoiir_gene <- spoiir_ko
  }, error = function(e) {})
  tryCatch({
    spoiir_gene$spo_gene_name <- "spoIIR"
  }, error = function(e) {})
  ### *** spoiir_gene

  tryCatch({
    ### spoivb ###
    spoivb_gene <- df %>% filter(grepl('spoivb', tolower(Preferred_name)))
    spoivb_gene$spo_gene_name <- "spoIVB"
    ### *** spoivb_gene
  }, error = function(e) {})

  tryCatch({
    ### spovt ###
    spovt_gene <- df %>% filter(grepl('spovt', tolower(Preferred_name)))
    spovt_gene$spo_gene_name <- "spoVT"
    ### *** spovt_gene
  }, error = function(e) {})

  tryCatch({
    ### dacf ###
    dacf_gene <- df %>% filter(grepl('dacf', tolower(Preferred_name)))
    dacf_gene$spo_gene_name <- "dacF"
    ### *** dacf_gene
  }, error = function(e) {})

  tryCatch({
    ### ytfj ###
    ytfj_gene <- df %>% filter(grepl('ytfj', tolower(Preferred_name)))
    ytfj_gene$spo_gene_name <- "ytfJ"
    ### *** ytfj_gene
  }, error = function(e) {})

  tryCatch({
    ### yhcv ###
    yhcv_gene <- df %>% filter(grepl('yhcv', tolower(Preferred_name)))
    yhcv_gene$spo_gene_name <- "yhcV"
    ### *** yhcv_gene
  }, error = function(e) {})

  tryCatch({
    ### yloC ###
    yloc_gene <- df %>% filter(grepl('yloc', tolower(Preferred_name)))
    yloc_gene$spo_gene_name <- "yloC"
    ### *** yloc_gene
  }, error = function(e) {})

  tryCatch({
    ### bofC ###
    bofc_gene <- df %>% filter(grepl('bofc', tolower(Preferred_name)))
    bofc_gene$spo_gene_name <- "bofC"
    ### *** bofc_gene
  }, error = function(e) {})

  tryCatch({
    ### rsfA ###
    rsfa_gene <- df %>% filter(grepl('rsfa', tolower(Preferred_name)))
    rsfa_gene$spo_gene_name <- "rsfA"
    ### *** rsfa_gene
  }, error = function(e) {})

  tryCatch({
    ### fin ###
    fin_gene <- df %>% filter(grepl('fin', tolower(Preferred_name)))
    yabk_gene <- df %>% filter(grepl('yabk', tolower(Preferred_name)))
    yabk_gene$spo_gene_name <- "fin/yabk"
    fin_gene <- rbind(fin_gene, yabk_gene)
  }, error = function(e) {})
  tryCatch({
    fin_gene$spo_gene_name <- "fin/yabk"
  }, error = function(e) {})
  ### *** fin_gene

  tryCatch({
    ### ymfJ ###
    ymfj_gene <- df %>% filter(grepl('ymfj', tolower(Preferred_name)))
    ymfj_gene$spo_gene_name <- "ymfJ"
    ### *** ymfj_gene
  }, error = function(e) {})

  tryCatch({
    ### yqhG ###
    yqhg_gene <- df %>% filter(grepl('yqhg', tolower(Preferred_name)))
    yqhg_gene$spo_gene_name <- "yqhG"
    ### *** yqhg_gene
  }, error = function(e) {})

  tryCatch({
    ### ywzB ###
    ywzb_gene <- df %>% filter(grepl('ywzb', tolower(Preferred_name)))
    ywzb_gene$spo_gene_name <- "ywzB"
    ### *** ywzb_gene
  }, error = function(e) {})

  tryCatch({
    sigF_regulon <- rbind(spoiir_gene, spoivb_gene)
    sigF_regulon <- rbind(sigF_regulon, spovt_gene)
    sigF_regulon <- rbind(sigF_regulon, dacf_gene)
    sigF_regulon <- rbind(sigF_regulon, ytfj_gene)
    sigF_regulon <- rbind(sigF_regulon, yhcv_gene)
    sigF_regulon <- rbind(sigF_regulon, yloc_gene)
    sigF_regulon <- rbind(sigF_regulon, bofc_gene)
    sigF_regulon <- rbind(sigF_regulon, rsfa_gene)
    sigF_regulon <- rbind(sigF_regulon, fin_gene)
    sigF_regulon <- rbind(sigF_regulon, ymfj_gene)
    sigF_regulon <- rbind(sigF_regulon, yqhg_gene)
    sigF_regulon <- rbind(sigF_regulon, ywzb_gene)
    sigF_regulon$spo_process <- "sigF_regulon"
  }, error = function(e) {})

  ##### sigG regulon #####
  message("Looking for genes related to SigG Regulon")

  ### spovac ###
  spovac_gene <- df %>% filter(grepl('spovac', tolower(Preferred_name)))
  tryCatch({
    spovac_ko <- df %>% filter(grepl('06405', tolower(KEGG_ko)))
    spovac_ko <- spovac_ko[spovac_ko$Preferred_name == "-",]
    spovac_gene <- rbind(spovac_gene, spovac_ko)
  }, error = function(e) {})
  tryCatch({
    spovac_gene$spo_gene_name <- "spoVAC"
  }, error = function(e) {})
  ### *** spovac_gene

  ### spovad ###
  spovad_gene <- df %>% filter(grepl('spovad', tolower(Preferred_name)))
  tryCatch({
    spovad_ko <- df %>% filter(grepl('06406', tolower(KEGG_ko)))
    spovad_ko <- spovad_ko[spovad_ko$Preferred_name == "-",]
    spovad_gene <- rbind(spovad_gene, spovad_ko)
  }, error = function(e) {})
  tryCatch({
    spovad_gene$spo_gene_name <- "spoVAD"
  }, error = function(e) {})
  ### *** spovad_gene

  ### spoVAEB ###
  spovaeb_gene <- df %>% filter(grepl('spovaeb', tolower(Preferred_name)))
  tryCatch({
    spovaeb_ko <- df %>% filter(grepl('06407', tolower(KEGG_ko)))
    spovaeb_gene <- spovaeb_ko
  }, error = function(e) {})
  tryCatch({
    spovaeb_gene$spo_gene_name <- "spoVAEB"
  }, error = function(e) {})
  ### *** spovaeb_gene

  ### nfo ###
  nfo_gene <- df %>% filter(grepl('nfo', tolower(Preferred_name)))
  tryCatch({
    nfo_ko <- df %>% filter(grepl('01151', tolower(KEGG_ko)))
    nfo_gene <- nfo_ko
  }, error = function(e) {})
  tryCatch({
    nfo_gene$spo_gene_name <- "nfo"
  }, error = function(e) {})
  ### *** nfo_gene

  tryCatch({
    ### mgla ###
    mgla_gene <- df %>% filter(grepl('mgla', tolower(Preferred_name)))
    sspa_gene <- df %>% filter(grepl('sspa', tolower(Preferred_name)))
    mgla_gene <- rbind(mgla_gene, sspa_gene)
    mgla_gene$spo_gene_name <- "mglA / sspA"
    ### *** mgla_gene
  }, error = function(e) {})

  ### sspb ###
  sspb_gene <- df %>% filter(grepl('sspb', tolower(Preferred_name)))
  tryCatch({
    sspb_ko <- df %>% filter(grepl('03600', tolower(KEGG_ko)))
    sspb_ko <- sspb_ko[sspb_ko$Preferred_name == "-",]
    sspb_gene <- rbind(sspb_gene, sspb_ko)
  }, error = function(e) {})
  tryCatch({
    sspb_gene$spo_gene_name <- "sspB"
  }, error = function(e) {})
  ### *** sspb_gene

  tryCatch({
    ### sspC ###
    sspc_gene <- df %>% filter(grepl('sspc', tolower(Preferred_name)))
    sspc_gene$spo_gene_name <- "sspC"
    ### *** sspc_gene
  }, error = function(e) {})

  tryCatch({
    ### spovaa ###
    spovaa_gene <- df %>% filter(grepl('spovaa', tolower(Preferred_name)))
    spovaa_gene$spo_gene_name <- "spoVAA"
    ### *** spovaa_gene
  }, error = function(e) {})

  ### spovab ###
  spovab_gene <- df %>% filter(grepl('spovab', tolower(Preferred_name)))
  tryCatch({
    spovab_ko <- df %>% filter(grepl('06404', tolower(KEGG_ko)))
    spovab_gene <- spovab_ko
  }, error = function(e) {})
  tryCatch({
    spovab_gene$spo_gene_name <- "spoVAB"
  }, error = function(e) {})
  ### *** spovab_gene

  tryCatch({
    ### spovaf ###
    spovaf_gene <- df %>% filter(grepl('spovaf', tolower(Preferred_name)))
    spovaf_gene$spo_gene_name <- "spoVAF"
    ### *** spovaf_gene
  }, error = function(e) {})

  ### sspf ###
  sspf_gene <- df %>% filter(grepl('sspf', tolower(Preferred_name)))
  tryCatch({
    sspf_ko <- df %>% filter(grepl('06423', tolower(KEGG_ko)))
    sspf_gene <- sspf_ko
  }, error = function(e) {})
  tryCatch({
    sspf_gene$spo_gene_name <- "sspF"
  }, error = function(e) {})
  ### *** sspf_gene

  ### ssph ###
  ssph_gene <- df %>% filter(grepl('ssph', tolower(Preferred_name)))
  tryCatch({
    ssph_ko <- df %>% filter(grepl('06425', tolower(KEGG_ko)))
    ssph_gene <- ssph_ko
  }, error = function(e) {})
  tryCatch({
    ssph_gene$spo_gene_name <- "sspH"
  }, error = function(e) {})
  ### *** ssph_gene

  tryCatch({
    ### sspi ###
    sspi_gene <- df %>% filter(grepl('sspi', tolower(Preferred_name)))
    sspi_gene$spo_gene_name <- "sspI"
    ### *** sspi_gene
  }, error = function(e) {})

  ### tlp ###
  tlp_gene <- df %>% filter(grepl('tlp', tolower(Preferred_name)))
  tryCatch({
    tlp_ko <- df %>% filter(grepl('06434', tolower(KEGG_ko)))
    tlp_gene <- tlp_ko
  }, error = function(e) {})
  tryCatch({
    tlp_gene$spo_gene_name <- "tlp"
  }, error = function(e) {})
  ### *** tlp_gene

  tryCatch({
    sigG_regulon <- rbind(spovac_gene, spovad_gene)
    sigG_regulon <- rbind(sigG_regulon, nfo_gene)
    sigG_regulon <- rbind(sigG_regulon, mgla_gene)
    sigG_regulon <- rbind(sigG_regulon, sspb_gene)
    sigG_regulon <- rbind(sigG_regulon, spovaa_gene)
    sigG_regulon <- rbind(sigG_regulon, spovaf_gene)
    sigG_regulon <- rbind(sigG_regulon, ssph_gene)
    sigG_regulon <- rbind(sigG_regulon, tlp_gene)
    sigG_regulon <- rbind(sigG_regulon, spovaeb_gene)
    sigG_regulon <- rbind(sigG_regulon, spovab_gene)
    sigG_regulon <- rbind(sigG_regulon, sspc_gene)
    sigG_regulon <- rbind(sigG_regulon, sspi_gene)
    sigG_regulon <- rbind(sigG_regulon, sspf_gene)
    sigG_regulon$spo_process <- "sigG_regulon"
  }, error = function(e) {})

  ##### sigE regulon #####
  message("Looking for genes related to SigE Regulon")

  tryCatch({
    ### sigk ###
    sigk_gene <- df %>% filter(grepl('sigk', tolower(Preferred_name)))
    sigk_gene$spo_gene_name <- "sigK"
    ### *** sigk_gene
  }, error = function(e) {})

  ### spoiiid ###
  spoiiid_gene <- df %>% filter(grepl('spoiiid', tolower(Preferred_name)))
  tryCatch({
    spoiiid_ko <- df %>% filter(grepl('06283', tolower(KEGG_ko)))
    spoiiid_gene <- spoiiid_ko
  }, error = function(e) {})
  tryCatch({
    spoiiid_gene$spo_gene_name <- "spoIIID"
  }, error = function(e) {})
  ### *** spoiiid_gene

  ### spoivfb ###
  spoivfb_gene <- df %>% filter(grepl('spoivfb', tolower(Preferred_name)))
  tryCatch({
    spoivfb_ko <- df %>% filter(grepl('06402', tolower(KEGG_ko)))
    spoivfb_ko <- spoivfb_ko[spoivfb_ko$Preferred_name == "-",]
    spoivfb_gene <- rbind(spoivfb_gene, spoivfb_ko)
  }, error = function(e) {})
  tryCatch({
    spoivfb_gene$spo_gene_name <- "spoIVFB"
  }, error = function(e) {})
  ### *** spoivfb_gene

  ### spovb ###
  spovb_gene <- df %>% filter(grepl('spovb', tolower(Preferred_name)))
  tryCatch({
    spovb_ko <- df %>% filter(grepl('06409', tolower(KEGG_ko)))
    spovb_gene <- spovb_ko
  }, error = function(e) {})
  tryCatch({
    spovb_gene$spo_gene_name <- "spoVB"
  }, error = function(e) {})
  ### *** spovb_gene

  ### ftsw ###
  ftsw_gene <- df %>% filter(grepl('ftsw', tolower(Preferred_name)))
  tryCatch({
    ftsw_ko <- df %>% filter(grepl('03588', tolower(KEGG_ko)))
    ftsw_gene <- ftsw_ko
  }, error = function(e) {})
  tryCatch({
    ftsw_gene$spo_gene_name <- "ftsW"
  }, error = function(e) {})
  ### *** ftsw_gene

  ### spoVK ###
  spoVK_gene <- df %>% filter(grepl('spovk', tolower(Preferred_name)))
  tryCatch({
    spoVK_ko <- df %>% filter(grepl('06413', tolower(KEGG_ko)))
    spoVK_ko <- spoVK_ko[spoVK_ko$Preferred_name == "-",]
    spoVK_gene <- rbind(spoVK_gene,spoVK_ko)
  }, error = function(e) {})
  tryCatch({
    spoVK_gene$spo_gene_name <- "spoVK"
  }, error = function(e) {})
  ### *** spoVK_gene

  ### ctpB ###
  ctpb_gene <- df %>% filter(grepl('ctpb', tolower(Preferred_name)))
  tryCatch({
    ctpb_ko <- df %>% filter(grepl('03797', tolower(KEGG_ko)))
    ctpb_gene <- ctpb_ko
  }, error = function(e) {})
  tryCatch({
    ctpb_gene$spo_gene_name <- "ctpB"
  }, error = function(e) {})
  ### *** ctpb_gene

  tryCatch({
    ### spoIVFA ###
    spoivfa_gene <- df %>% filter(grepl('spoivfa', tolower(Preferred_name)))
    spoivfa_gene$spo_gene_name <- "spoIVFA"
    ### *** spoivfa_gene
  }, error = function(e) {})

  ### bofA ###
  bofa_gene <- df %>% filter(grepl('bofa', tolower(Preferred_name)))
  tryCatch({
    bofa_ko <- df %>% filter(grepl('06317', tolower(KEGG_ko)))
    bofa_gene <- bofa_ko
  }, error = function(e) {})
  tryCatch({
    bofa_gene$spo_gene_name <- "bofA"
  }, error = function(e) {})
  ### *** bofa_gene

  tryCatch({
    ### dacb ###
    dacb_gene <- df %>% filter(grepl('dacb', tolower(Preferred_name)))
    dacb_gene$spo_gene_name <- "dacB"
    ### *** dacb_gene
  }, error = function(e) {})

  tryCatch({
    ### alr ###
    alr_gene <- df %>% filter(grepl('alr', tolower(Preferred_name)))
    alr_gene$spo_gene_name <- "alr"
    ### *** alr_gene
  }, error = function(e) {})

  tryCatch({
    ### spmA ###
    spmA_gene <- df %>% filter(grepl('spma', tolower(Preferred_name)))
    spmA_gene$spo_gene_name <- "spmA"
    ### *** spmA_gene
  }, error = function(e) {})

  ### spmB ###
  spmB_gene <- df %>% filter(grepl('spmb', tolower(Preferred_name)))
  tryCatch({
    spmB_ko <- df %>% filter(grepl('06374', tolower(KEGG_ko)))
    spmB_gene <- spmB_ko
  }, error = function(e) {})
  tryCatch({
    spmB_gene$spo_gene_name <- "spmB"
  }, error = function(e) {})
  ### *** spmB_gene

  tryCatch({
    ### yisY ###
    yisY_gene <- df %>% filter(grepl('yisy', tolower(Preferred_name)))
    yisY_gene$spo_gene_name <- "yisY"
    ### *** yisY_gene
  }, error = function(e) {})

  tryCatch({
    ### yqfu ###
    yqfu_gene <- df %>% filter(grepl('yqfu', tolower(Preferred_name)))
    yqfu_gene$spo_gene_name <- "yqfU"
    ### *** yqfu_gene
  }, error = function(e) {})

  tryCatch({
    ### ylmc ###
    ylmc_gene <- df %>% filter(grepl('ylmc', tolower(Preferred_name)))
    ylmc_gene$spo_gene_name <- "ylmC"
    ### *** ylmc_gene
  }, error = function(e) {})

  tryCatch({
    ### ytaf ###
    ytaf_gene <- df %>% filter(grepl('ytaf', tolower(Preferred_name)))
    ytaf_gene$spo_gene_name <- "ytaF"
    ### *** ytaf_gene
  }, error = function(e) {})

  tryCatch({
    ### ytvi ###
    ytvi_gene <- df %>% filter(grepl('ytvi', tolower(Preferred_name)))
    ytvi_gene$spo_gene_name <- "ytvI"
    ### *** ytvi_gene
  }, error = function(e) {})

  tryCatch({
    ### ykvi ###
    ykvi_gene <- df %>% filter(grepl('ykvi', tolower(Preferred_name)))
    ykvi_gene$spo_gene_name <- "ykvI"
    ### *** ykvi_gene
  }, error = function(e) {})

  tryCatch({
    ### ydca ###
    ydca_gene <- df %>% filter(grepl('ydca', tolower(Preferred_name)))
    ydca_gene$spo_gene_name <- "ydca"
    ### *** ydca_gene
  }, error = function(e) {})

  tryCatch({
    ### ydcc ###
    ydcc_gene <- df %>% filter(grepl('ydcc', tolower(Preferred_name)))
    ydcc_gene$spo_gene_name <- "ydcc"
    ### *** ydcc_gene
  }, error = function(e) {})

  tryCatch({
    ### yhbh ###
    yhbh_gene <- df %>% filter(grepl('yhbh', tolower(Preferred_name)))
    yhbh_gene$spo_gene_name <- "yhbh"
    ### *** yhbh_gene
  }, error = function(e) {})

  tryCatch({
    sigE_regulon <- rbind(sigk_gene, spoiiid_gene)
    sigE_regulon <- rbind(sigE_regulon, spoivfb_gene)
    sigE_regulon <- rbind(sigE_regulon, spovb_gene)
    sigE_regulon <- rbind(sigE_regulon, spoVK_gene)
    sigE_regulon <- rbind(sigE_regulon, ftsw_gene)
    sigE_regulon <- rbind(sigE_regulon, bofa_gene)
    sigE_regulon <- rbind(sigE_regulon, alr_gene)
    sigE_regulon <- rbind(sigE_regulon, dacb_gene)
    sigE_regulon <- rbind(sigE_regulon, spmA_gene)
    sigE_regulon <- rbind(sigE_regulon, spmB_gene)
    sigE_regulon <- rbind(sigE_regulon, yisY_gene)
    sigE_regulon <- rbind(sigE_regulon, ylmc_gene)
    sigE_regulon <- rbind(sigE_regulon, ytaf_gene)
    sigE_regulon <- rbind(sigE_regulon, ytvi_gene)
    sigE_regulon <- rbind(sigE_regulon, ctpb_gene)
    sigE_regulon <- rbind(sigE_regulon, spoivfa_gene)
    sigE_regulon <- rbind(sigE_regulon, ykvi_gene)
    sigE_regulon <- rbind(sigE_regulon, yqfu_gene)
    sigE_regulon <- rbind(sigE_regulon, ydca_gene)
    sigE_regulon <- rbind(sigE_regulon, ydcc_gene)
    sigE_regulon <- rbind(sigE_regulon, yhbh_gene)
    sigE_regulon$spo_process <- "sigE_regulon"
  }, error = function(e) {})

  ##### sigK regulon #####
  message("Looking for genes related to SigK Regulon")

  ### spovfa ###
  spovfa_gene <- df %>% filter(grepl('spovfa', tolower(Preferred_name)))
  tryCatch({
    spovfa_ko <- df %>% filter(grepl('06410', tolower(KEGG_ko)))
    spovfa_gene <- spovfa_ko
  }, error = function(e) {})
  tryCatch({
    spovfa_gene$spo_gene_name <- "spoVFA"
  }, error = function(e) {})
  ### *** spovfa_gene

  tryCatch({
    spovfb_gene <- df %>% filter(grepl('spovfb', tolower(Preferred_name)))
    spovfb_gene$spo_gene_name <- "spoVFB"
    ### *** spovfb_gene
  }, error = function(e) {})

  tryCatch({
    ### ykud ###
    ykud_gene <- df %>% filter(grepl('ykud', tolower(Preferred_name)))
    ykud_gene$spo_gene_name <- "ykuD"
    ### *** ykud_gene
  }, error = function(e) {})

  tryCatch({
    sigK_regulon <- rbind(spovfa_gene, spovfb_gene)
    sigK_regulon <- rbind(sigK_regulon, ykud_gene)
    sigK_regulon$spo_process <- "sigK_regulon"
  }, error = function(e) {})


  ##### spore cortex #####
  message("Looking for genes related to Spore Cortex")

  tryCatch({
    ### spovd ###
    spovd_gene <- df %>% filter(grepl('spovd', tolower(Preferred_name)))
    spovd_gene$spo_gene_name <- "spoVD"
    ### *** spovd_gene
  }, error = function(e) {})

  tryCatch({
    ### ylbj ###
    ylbj_gene <- df %>% filter(grepl('ylbj', tolower(Preferred_name)))
    ylbj_gene$spo_gene_name <- "ylbJ"
    ### *** ylbj_gene
  }, error = function(e) {})

  tryCatch({
    ### cwlc ###
    cwlc_gene <- df %>% filter(grepl('cwlc', tolower(Preferred_name)))
    cwlc_gene$spo_gene_name <- "cwlC"
  }, error = function(e) {})

  tryCatch({
    cwld_gene <- df %>% filter(grepl('cwld', tolower(Preferred_name)))
    cwld_gene$spo_gene_name <- "cwlD"
  }, error = function(e) {})

  tryCatch({
    cwldc_gene <- rbind(cwlc_gene, cwld_gene)
    cwldc_gene$spo_gene_name <- "cwlC/cwlD"
    ### *** cwldc_gene
  }, error = function(e) {})

  tryCatch({
    ### lyth ###
    lyth_gene <- df %>% filter(grepl('lyth', tolower(Preferred_name)))
    lyth_gene$spo_gene_name <- "lytH"
    ### *** lyth_gene
  }, error = function(e) {})

  tryCatch({
    ### yabp ###
    yabp_gene <- df %>% filter(grepl('yabp', tolower(Preferred_name)))
    yabp_gene$spo_gene_name <- "yabP"
    ### *** yabp_gene
  }, error = function(e) {})

  tryCatch({
    ### yabq ###
    yabq_gene <- df %>% filter(grepl('yabq', tolower(Preferred_name)))
    yabq_gene$spo_gene_name <- "yabQ"
    ### *** yabq_gene
  }, error = function(e) {})

  tryCatch({
    ### yqfc ###
    yqfc_gene <- df %>% filter(grepl('yqfc', tolower(Preferred_name)))
    yqfc_gene$spo_gene_name <- "yqfC"
    ### *** yqfc_gene
  }, error = function(e) {})

  tryCatch({
    ### yqfD ###
    yqfD_gene <- df %>% filter(grepl('yqfd', tolower(Preferred_name)))
    yqfD_gene$spo_gene_name <- "yqfD"
    ### *** yqfD_gene
  }, error = function(e) {})

  ### cotd ###
  tryCatch({
    cotd_gene <- df %>% filter(grepl('cotd', tolower(Preferred_name)))
    cotd_ko <- df %>% filter(grepl('06327', tolower(KEGG_ko)))
    cotd_gene <- cotd_ko
    cotd_gene$spo_gene_name <- "cotD"
  }, error = function(e) {})
  ### *** cotd_gene

  tryCatch({
    spore_cortex <- rbind(spovd_gene, ylbj_gene)
    spore_cortex <- rbind(spore_cortex, cwldc_gene)
    spore_cortex <- rbind(spore_cortex, lyth_gene)
    spore_cortex <- rbind(spore_cortex, yabp_gene)
    spore_cortex <- rbind(spore_cortex, yabq_gene)
    spore_cortex <- rbind(spore_cortex, yqfc_gene)
    spore_cortex <- rbind(spore_cortex, yqfD_gene)
    spore_cortex <- rbind(spore_cortex, cotd_gene)
    spore_cortex$spo_process <- "spore_cortex"
  }, error = function(e) {})

  ##### spore coat #####
  message("Looking for genes related to Spore Coat")

  ### spoIVA ###
  spoIVA_gene <- df %>% filter(grepl('spoiva', tolower(Preferred_name)))
  tryCatch({
    spoIVA_ko <- df %>% filter(grepl('06398', tolower(KEGG_ko)))
    spoIVA_gene <- spoIVA_ko
  }, error = function(e) {})
  tryCatch({
    spoIVA_gene$spo_gene_name <- "spoIVA"
  }, error = function(e) {})
  ### *** spoIVA_gene

  ### cotjc ###
  cotjc_gene <- df %>% filter(grepl('cotjc', tolower(Preferred_name)))
  tryCatch({
    cotjc_ko <- df %>% filter(grepl('06334', tolower(KEGG_ko)))
    cotjc_gene <- cotjc_ko
  }, error = function(e) {})
  tryCatch({
    cotjc_gene$spo_gene_name <- "cotJC"
  }, error = function(e) {})
  ### *** cotjc_gene

  tryCatch({
    ### cotsa ###
    cotsa_gene <- df %>% filter(grepl('cotsa', tolower(Preferred_name)))
    cotsa_gene$spo_gene_name <- "cotSA"
    ### *** cotsa_gene
  }, error = function(e) {})

  ### gerM ###
  gerM_gene <- df %>% filter(grepl('germ', tolower(Preferred_name)))
  tryCatch({
    gerM_ko <- df %>% filter(grepl('06298', tolower(KEGG_ko)))
    gerM_gene <- gerM_ko
  }, error = function(e) {})
  tryCatch({
    gerM_gene$spo_gene_name <- "gerM"
  }, error = function(e) {})
  ### *** gerM_gene

  tryCatch({
    ### lipc ###
    lipc_gene <- df %>% filter(grepl('lipc', tolower(Preferred_name)))
    lipc_gene$spo_gene_name <- "lipC"
    ### *** lipc_gene
  }, error = function(e) {})

  tryCatch({
    ### safa ###
    safa_gene <- df %>% filter(grepl('safa', tolower(Preferred_name)))
    safa_gene$spo_gene_name <- "safA"
    ### *** safa_gene
  }, error = function(e) {})

  tryCatch({
    ### ydhd ###
    ydhd_gene <- df %>% filter(grepl('ydhd', tolower(Preferred_name)))
    ydhd_gene$spo_gene_name <- "ydhD"
    ### *** ydhd_gene
  }, error = function(e) {})

  tryCatch({
    ### yhax ###
    yhax_gene <- df %>% filter(grepl('yhax', tolower(Preferred_name)))
    yhax_gene$spo_gene_name <- "yhaX"
    ### *** yhax_gene
  }, error = function(e) {})

  ### spovid ###
  spovid_gene <- df %>% filter(grepl('spovid', tolower(Preferred_name)))
  tryCatch({
    spovid_ko <- df %>% filter(grepl('06417', tolower(KEGG_ko)))
    spovid_gene <- spovid_ko
  }, error = function(e) {})
  tryCatch({
    spovid_gene$spo_gene_name <- "spoVID"
  }, error = function(e) {})
  ### *** spovid_gene

  tryCatch({
    ### spovif ###
    spovif_gene <- df %>% filter(grepl('spovif', tolower(Preferred_name)))
    spovif_gene$spo_gene_name <- "spoVIF"
    ### *** spovif_gene
  }, error = function(e) {})

  ### cote ###
  cote_gene <- df %>% filter(grepl('cote', tolower(Preferred_name)))
  tryCatch({
    cote_ko <- df %>% filter(grepl('06328', tolower(KEGG_ko)))
    cote_gene <- cote_ko
  }, error = function(e) {})
  tryCatch({
    cote_gene$spo_gene_name <- "cotE"
  }, error = function(e) {})
  ### *** cote_gene

  tryCatch({
    ### cotja ###
    cotja_gene <- df %>% filter(grepl('cotja', tolower(Preferred_name)))
    cotja_gene$spo_gene_name <- "cotJA"
    ### *** cotja_gene
  }, error = function(e) {})

  ### cotjb ###
  cotjb_gene <- df %>% filter(grepl('cotjb', tolower(Preferred_name)))
  tryCatch({
    cotjb_ko <- df %>% filter(grepl('06333', tolower(KEGG_ko)))
    cotjb_gene <- cotjb_ko
  }, error = function(e) {})
  tryCatch({
    cotjb_gene$spo_gene_name <- "cotJB"
  }, error = function(e) {})
  ### *** cotjb_gene

  ### cotp ###
  cotp_gene <- df %>% filter(grepl('cotp', tolower(Preferred_name)))
  tryCatch({
    cotp_ko <- df %>% filter(grepl('13993', tolower(KEGG_ko)))
    cotp_gene <- cotp_ko
  }, error = function(e) {})
  tryCatch({
    cotp_gene$spo_gene_name <- "hsps"
  }, error = function(e) {})
  ### *** cotp_gene

  tryCatch({
    ### yhjr ###
    yhjr_gene <- df %>% filter(grepl('yhjr', tolower(Preferred_name)))
    yhjr_gene$spo_gene_name <- "yhjR"
    ### *** yhjr_gene
  }, error = function(e) {})

  tryCatch({
    spore_coat <- rbind(spoIVA_gene, cotjc_gene)
    spore_coat <- rbind(spore_coat, cotsa_gene)
    spore_coat <- rbind(spore_coat, gerM_gene)
    spore_coat <- rbind(spore_coat, safa_gene)
    spore_coat <- rbind(spore_coat, ydhd_gene)
    spore_coat <- rbind(spore_coat, yhax_gene)
    spore_coat <- rbind(spore_coat, cotjb_gene)
    spore_coat <- rbind(spore_coat, yhjr_gene)
    spore_coat <- rbind(spore_coat, lipc_gene)
    spore_coat <- rbind(spore_coat, spovid_gene)
    spore_coat <- rbind(spore_coat, spovif_gene)
    spore_coat <- rbind(spore_coat, cote_gene)
    spore_coat <- rbind(spore_coat, cotja_gene)
    spore_coat <- rbind(spore_coat, cotp_gene)
    spore_coat$spo_process <- "spore_coat"
  }, error = function(e) {})

  ##### germination #####
  message("Looking for genes related to Germination")

  tryCatch({
    ### gera ###
    gera_gene <- df %>% filter(grepl('gera', tolower(Preferred_name)))
    gera_gene$spo_gene_name <- "gerA"
    ### *** gera_gene
  }, error = function(e) {})

  ### gerc ###
  gerc_gene <- df %>% filter(grepl('gerc', tolower(Preferred_name)))
  tryCatch({
    gerc_ko <- df %>% filter(grepl('06308', tolower(KEGG_ko)))
    gerc_gene <- gerc_ko
  }, error = function(e) {})
  tryCatch({
    gerc_gene$spo_gene_name <- "gerC"
  }, error = function(e) {})
  ### *** gerc_gene

  ### lgt ###
  lgt_gene <- df %>% filter(grepl('lgt', tolower(Preferred_name)))
  tryCatch({
    lgt_ko <- df %>% filter(grepl('13292', tolower(KEGG_ko)))
    lgt_gene <- lgt_ko
  }, error = function(e) {})
  tryCatch({
    lgt_gene$spo_gene_name <- "lgt"
  }, error = function(e) {})
  ### *** lgt_gene

  ### gpr ###
  gpr_gene <- df %>% filter(grepl('gpr', tolower(Preferred_name)))
  tryCatch({
    gpr_ko <- df %>% filter(grepl('06012', tolower(KEGG_ko)))
    gpr_gene <- gpr_ko
  }, error = function(e) {})
  tryCatch({
    gpr_gene$spo_gene_name <- "gpr"
  }, error = function(e) {})
  ### *** gpr_gene

  tryCatch({
    ### cwlj ###
    cwlj_gene <- df %>% filter(grepl('cwlj', tolower(Preferred_name)))
    cwlj_gene$spo_gene_name <- "cwlJ"
  }, error = function(e) {})

  tryCatch({
    sleb_gene <- df %>% filter(grepl('sleb', tolower(Preferred_name)))
    sleb_gene$spo_gene_name <- "sleB"
  }, error = function(e) {})

  tryCatch({
    cwlj_gene <- rbind(cwlj_gene, sleb_gene)
    cwlj_gene$spo_gene_name <- "cwlJ/sleB"
    ### *** cwlj_gene
  }, error = function(e) {})

  ### cspa ###
  cspa_gene <- df %>% filter(grepl('cspa', tolower(Preferred_name)))
  tryCatch({
    cspa_ko <- df %>% filter(grepl('03704', tolower(KEGG_ko)))
    cspa_gene <- cspa_ko
  }, error = function(e) {})
  tryCatch({
    cspa_gene$spo_gene_name <- "CSD"
  }, error = function(e) {})
  ### *** cspa_gene

  tryCatch({
    ### gdh ###
    gdh_gene <- df[df$Preferred_name == 'gdh',]
    gdh_gene$spo_gene_name <- "gdh"
    ### *** gdh_gene
  }, error = function(e) {})

  tryCatch({
    ### gerd ###
    gerd_gene <- df %>% filter(grepl('gerd', tolower(Preferred_name)))
    gerd_gene$spo_gene_name <- "gerD"
    ### *** gerd_gene
  }, error = function(e) {})

  tryCatch({
    ### gere ###
    gere_gene <- df %>% filter(grepl('gere', tolower(Preferred_name)))
    gere_gene$spo_gene_name <- "gerE"
    ### *** gere_gene
  }, error = function(e) {})

  tryCatch({
    ### gerq ###
    gerq_gene <- df %>% filter(grepl('gerq', tolower(Preferred_name)))
    gerq_gene$spo_gene_name <- "gerQ"
    ### *** gerq_gene
  }, error = function(e) {})

  tryCatch({
    ### ypeb ###
    ypeb_gene <- df %>% filter(grepl('ypeb', tolower(Preferred_name)))
    ypeb_gene$spo_gene_name <- "ypeB"
    ### *** ypeb_gene
  }, error = function(e) {})

  tryCatch({
    germination <- rbind(gera_gene, gpr_gene)
    germination <- rbind(germination, cspa_gene)
    germination <- rbind(germination, gdh_gene)
    germination <- rbind(germination, ypeb_gene)
    germination <- rbind(germination, gerc_gene)
    germination <- rbind(germination, lgt_gene)
    germination <- rbind(germination, gerd_gene)
    germination <- rbind(germination, gere_gene)
    germination$spo_process <- "germination"
  }, error = function(e) {})

  tryCatch({
    df <- rbind(sporulation_onset_and_checkpoints, spo0a_regulon)
    df <- rbind(df, Engulfment)
    df <- rbind(df, sigE_regulon)
    df <- rbind(df, sigF_regulon)
    df <- rbind(df, sigG_regulon)
    df <- rbind(df, sigK_regulon)
    df <- rbind(df, spore_coat)
    df <- rbind(df, spore_cortex)
    df <- rbind(df, germination)
  }, error = function(e) {})

  df <- df %>% dplyr::filter(!is.na(genome_ID))

  return(df)
}
utils::globalVariables(c(
  "genome_ID", "spo_gene_name", "present", "rf_model",
  "svm_model", "meta_model", "Preferred_name", "KEGG_ko",
  "spoiiie_ko_eccca"
))
