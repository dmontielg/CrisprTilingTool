

#' Locate base modifications in a sgRNA window and apply base sustition
#'
#' @param df.target.trans is a dataframe transcript from the CDS human genome
#' @param gene.search gene name to look in the CDS human genome
#' @param trans transcript name to look in the CDS human genome
#' @param flank number of base pairs to spam before and after each exon
#' @param sgrna_len size of the Single Guided RNA to cover
#' @param window_start index of the window within the sgRNA
#' @param window_end  end of the index within the sgRNA
#' @param base_change base change to find and to convert e.g. A->G, C->T
#' @param seq_before Sequence to append before the sgRNA (protospacer)
#' @param seq_after Sequence to append after the sgRNA (protospacer)
#'
#' @return a dataframe object containing base changes and their positions
#' @export
#'
#' @examples
run_locate_modifications <- function(df.target.trans,
                                     gene.search,
                                     trans,flank,
                                     sgrna_len,
                                     window_start,
                                     window_end,
                                     base_change,
                                     seq_before,
                                     seq_after,
                                     Genome) {
  
  
  codons.list <- c("AGA","AGC","AGT","AGG",
                   "CGA","CGC","CGT","CGG",
                   "TGA","TGC","TGT","TGG",
                   "GGA","GGC","GGT","GGG")
  window_len <- window_end - window_start + 1
  n <- 3 # codon len
  base_current <- NULL
  base_mod <- NULL
  base_current_rev <- NULL
  base_mod_rev <- NULL
  
  for (base in stringr::str_split(base_change, "->")) {
    base_current <- base[1]
    base_current_rev <- as.character(
      Biostrings::complement(
        Biostrings::DNAString(base[1])))
    base_mod <- base[2]
    base_mod_rev <- as.character(
      Biostrings::complement(
        Biostrings::DNAString(base[2])))
  }
  
  trans_sequence <- NULL
  aa_trans_sequence <- NULL
  df_dna_info <- NULL
  trans_sequence <- as.character(
    get_dna_trasncript_sequence(df.target.trans, Genome))
  aa_trans_sequence <- as.character(Biostrings::translate(
    Biostrings::DNAString(trans_sequence)))
  
  characters_dna <- Biostrings::strsplit(repeat_characters(
    trans_sequence, 1), "")[[1]]
  df_dna_info <- as.data.frame(BiocGenerics::cbind(characters_dna, 
                                                   1:nchar(trans_sequence)))
  characters_aa <- Biostrings::strsplit(repeat_characters(
    aa_trans_sequence, 1), "")[[1]]
  df_aa_info <- as.data.frame(BiocGenerics::cbind(characters_aa,
                                                  1:nchar(aa_trans_sequence)))
  colnames(df_dna_info) <- c("dna","pos")
  df_dna_info$pos <- as.numeric(df_dna_info$pos)
  df.phase.exons.trans <- get_phase_transcript(df.target.trans, Genome)
  df.phase.exons.trans$complete <- unlist(
    complete_phase_sequences(df.phase.exons.trans))
  df.phase.exons.trans$exon.seq <- unlist(
    get_exon_phase_sequences(df.phase.exons.trans))
  df.aa.exons.index <- get_absolute_coordinates_aa(df.phase.exons.trans)
  
  # Reverse strand
  df.phase.exons.trans.rev <- assign_phase_reverse(df.phase.exons.trans)
  df_table_all <- list()
  
  for (exon_n in unique(df.target.trans$exon)) {
    df_table <- NULL
    df.target.trans.exon <- df.target.trans %>% dplyr::filter(exon == exon_n )
    # if exon_n is not the last 
    current.exon.seq <- Biostrings::getSeq(Genome, 
                                           df.target.trans.exon$chr, 
                                           start = df.target.trans.exon$start - flank, 
                                           end = df.target.trans.exon$end + flank,
                                           strand = df.target.trans.exon$strand)
    # Forward strand
    forward_out <- get_pam_matches_info(codons.list, 
                                        current.exon.seq, 
                                        "+", 
                                        sgrna_len,
                                        window_start,
                                        window_end)
    
    match_positions <- lapply(as.vector(unlist(forward_out$sgRNA)), function(x) 
      as.data.frame(stringr::str_locate_all(as.character(current.exon.seq),x)[[1]]))
    
    df_exon <- get_dataframe_sequence_exon_phase(df.phase.exons.trans, 
                                                 exon_n, 
                                                 flank)
    df_exon_info <- df_exon[[1]]
    df_exon_info_tmp <- df_exon[[2]]
    
    df_out_window_change_forward <- get_amino_acid_changes_forward(
      match_positions, 
      window_start, 
      window_end,
      df_exon_info, 
      df_exon_info_tmp, 
      base_current, 
      base_mod,
      flank)
    df_out_forward <- BiocGenerics::cbind(as.data.frame(forward_out),
                                          as.data.frame(df_out_window_change_forward))
    # Reverse strand
    current.exon.seq.rev <- Biostrings::reverseComplement(current.exon.seq)
    reverse_out <- get_pam_matches_info(codons.list, 
                                        current.exon.seq.rev, 
                                        "-", 
                                        sgrna_len,
                                        window_start,
                                        window_end)
    match_positions <- lapply(as.vector(unlist(reverse_out$sgRNA)), function(x) 
      as.data.frame(stringr::str_locate_all(
        as.character(current.exon.seq.rev), x)[[1]]))
    df_exon <- get_dataframe_sequence_exon_phase(df.phase.exons.trans.rev, 
                                                 exon_n, 
                                                 flank)
    df_exon_info_tmp <- df_exon[[2]]
    phase_rev <- df.phase.exons.trans.rev %>% filter(exon == exon_n)
    df_out_window_change_reverse <- get_amino_acid_changes_reverse(
      match_positions, 
      window_start, 
      window_end,
      df_exon_info, 
      df_exon_info_tmp, 
      base_current_rev, 
      base_mod_rev,
      phase_rev,
      flank)
    df_out_reverse <- BiocGenerics::cbind(as.data.frame(reverse_out), 
                                          as.data.frame(df_out_window_change_reverse))
    df_table <- BiocGenerics::rbind(df_out_forward, df_out_reverse)
    df_table$exon <- rep(exon_n, as.numeric(nrow(df_table)))
    df_table_all <- BiocGenerics::rbind(df_table_all, df_table)
  }
  
  df_table_all$out_aa_info <- gsub("\\*", "Ter", df_table_all$out_aa_info)
  df_table_all$tx <- unique(df.target.trans$tx)
  df_table_all$chr <- unique(df.target.trans$chr)
  df_table_all$gene <- unique(df.target.trans$gene)
  
  out_aa_info_abs <- get_map_coordinates(df_table_all, 
                                         df.phase.exons.trans, 
                                         df.aa.exons.index)
  
  df_table_all$out_aa_info_abs <- unlist(out_aa_info_abs)
  df_table_all$protospacer <- paste(seq_before,
                                    as.vector(unlist(df_table_all$sgRNA)),
                                    seq_after, sep = "-")
  df_table_all <- apply(df_table_all, 2 ,as.character)
  
  filename_out <- paste("Table", gene.search, "mod", base_current, base_mod,
                        "flank", flank, ".txt", sep = "_")
  return(df_table_all)
}

