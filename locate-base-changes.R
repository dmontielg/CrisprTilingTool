
get_aa_modification_info <- function(current_match, 
                                     df_exon_info, 
                                     base_current, 
                                     base_mod){ 
  # current_match: Contains the current window and the chromosomical position
  # df_exon_info: Contains the chromosomical position for the entire transcript
  # base_current base to look for within the window
  # base_mod base to modify once is found within the window
  
  out_dna_current  <- NULL
  out_dna_mod      <- NULL
  out_aa_current   <- NULL
  out_aa_mod       <- NULL
  current_seq      <- NULL
  current_seq_mod  <- NULL
  
  current_match_mod <- current_match
  col_name_mod <- paste("dna", base_current, base_mod, sep = "_")
  current_match_mod[col_name_mod] <- replace(current_match$dna, 
                                    current_match$dna == base_current, base_mod)
  # Verify is there is any modification
  flag_base_change <- as.vector(
                    current_match_mod[col_name_mod] == as.vector(current_match$dna)
                    )
  df_dna <- current_match_mod[, c("dna", "pos")]
  
  df_dna_mod <- current_match_mod[, c(col_name_mod, "pos")]
  
  out_aa_mod <- get_aa_info(current_match, 
                          flag_base_change, 
                          df_dna_mod, 
                          df_exon_info)
  
  out_dna_mod <- get_codon_mod_info(df_dna, df_dna_mod)
  out_dna_window <- paste(current_match$dna, collapse = "")
  
  return(list(out_dna_window = out_dna_window,
              out_dna_mod = out_dna_mod,
              out_aa_mod = out_aa_mod))
}


get_aa_info <- function(out_dna, 
                      flag_base_mod, 
                      df_current_match, 
                      df_exon_info){
  
  out_aa_info_str <- ""
  if (FALSE %in% flag_base_mod) { # if there is a base modification 
    
    out_dna_fill <- fill_codon_bases(
      out_dna,
      df_exon_info)
    
    out_dna_base_mod <- fill_codon_bases(
      df_current_match, 
      df_exon_info)
    Biostrings::translate(Biostrings::DNAString(
      paste(
        as.character(out_dna_fill$dna),collapse = "")
    ))
    Biostrings::translate(Biostrings::DNAString(
      paste(
        as.character(out_dna_base_mod$dna),collapse = "")
    ))
    flag <- as.vector(out_dna_base_mod$dna == out_dna_fill$dna)
    base_i <- 3
    out_base_mod <- NULL
    for (codon_i in 1:(length(flag)/3)){
      #codon_i <- 1
      base_start <- as.integer(base_i - 2)
      base_end <- base_i
      if (FALSE %in% flag[base_start : base_end ]) {
        # If False THERE is a base change in this codon
        current_codon_mod <- out_dna_base_mod[base_start : base_end ,]
        current_aa_mod <- Biostrings::translate(Biostrings::DNAString(paste(
          as.character(current_codon_mod$dna)
          , collapse = "")
        ), no.init.codon = T)
        
        #current_aa_mod_i <- current_codon_mod$pos[nrow(current_codon_mod) ] / 3
        repeated_sequence <- rep(1:as.integer(nrow(df_exon_info) / 3), each = 3)
        current_aa_mod_i <- as.numeric(rownames(df_exon_info[which(
          df_exon_info$pos == min(current_codon_mod$pos)),]))
        current_aa_mod_i <- repeated_sequence[current_aa_mod_i]
        current_aa <- out_dna_fill %>% 
                    dplyr::filter(pos >= min(current_codon_mod$pos) &
                             pos <= max(current_codon_mod$pos))
        
        current_aa <- Biostrings::translate(
                        Biostrings::DNAString(
                          paste(current_aa$dna, 
                                collapse = "")),
                                no.init.codon = T)
        out_aa <- paste(current_aa, current_aa_mod_i,current_aa_mod, sep = "")
        out_aa_info_str <- paste(out_aa_info_str, out_aa, sep = "")
      }
      base_i <- base_i + 3
    }
    return(list(out_aa_info = out_aa_info_str))
  }else{ # there is nothing to modified
    return(list(out_aa_info = NA))
  }
}


#' Fill missing bases to complete a codon before and after. 
#'
#' @param current_match dataframe of each base and position where the window is
#' @param df_exon_info dataframe of each base and position of the exon
#'
#' @return
#' @export
#'
#' @examples
fill_codon_bases <- function(current_match, df_exon_info){
  
  result <- NULL
  current_match_mod <- NULL
  val_check <- as.numeric(rownames(
                          df_exon_info[which(
                            df_exon_info$pos == min(current_match$pos)),]))
  # AAAC -> AAAC[CC]
  if(val_check %% 3 == 1){
    # AAACC -> AAACCC
    # If there is only one bp to modified 
    # A -> A[AA]
    if (nrow(current_match) == 1) {
      
      current_match_mod <- df_exon_info %>% 
        dplyr::filter(pos >= min(current_match$pos) & 
                        pos <= (max(current_match$pos) + 2))
      
    }else{
      
      if(nrow(current_match) %% 3 == 1){
        current_match_mod <- df_exon_info %>% 
          dplyr::filter(pos >= min(current_match$pos) & 
                          pos <= (max(current_match$pos) + 2))    
      }else if(nrow(current_match) %% 3 == 2){
        current_match_mod <- df_exon_info %>% 
          dplyr::filter(pos >= min(current_match$pos) & 
                          pos <= (max(current_match$pos) + 1))    
      }
    }
    
    result <- fill_missing_positions(current_match, current_match_mod)
  }else if(val_check %% 3 == 2){
    # AABBB -> AAABBB
    # If there is only one bp to modified 
    # A -> [A]A[A]
    if (nrow(current_match) == 1) {
        current_match_mod <- df_exon_info %>% 
                dplyr::filter(pos >= min(current_match$pos) - 1 & 
                              pos <= (max(current_match$pos)) + 1)
    }else{
        current_match_mod <- df_exon_info %>% 
                dplyr::filter(pos >= min(current_match$pos) - 1 & 
                              pos <= (max(current_match$pos)))
        if(nrow(current_match_mod) %% 3 == 1){
                current_match_mod <- df_exon_info %>% 
                dplyr::filter(pos >= min(current_match$pos) - 1 & 
                              pos <= (max(current_match$pos) + 2))
        }else if(nrow(current_match_mod) %% 3 == 2){
          current_match_mod <- df_exon_info %>% 
                dplyr::filter(pos >= min(current_match$pos) - 1 & 
                              pos <= (max(current_match$pos) + 1))
        }
    }
    
    result <- fill_missing_positions(current_match, current_match_mod)
    
  }else if(val_check %% 3 == 0){
    # ACCCB -> [AA]ACCCB[BB]
    # If there is only one bp to modified 
    # A -> [AA]A
    if (nrow(current_match) == 1) {
      current_match_mod <- df_exon_info %>% 
                      dplyr::filter(pos >= min(current_match$pos) - 2 & 
                                    pos <= (max(current_match$pos)))
    }else{
      current_match_mod <- df_exon_info %>% 
                      dplyr::filter(pos >= min(current_match$pos) - 2 & 
                                    pos <= (max(current_match$pos) + 2))
    }
    result <- fill_missing_positions(current_match, current_match_mod)
  }
  return(result)
}


fill_missing_positions <- function(df1, df2){
  # Context: Window is the region found within sgRNA
  # Identify the unique positions on two data frames windows
  # In this case df1 is current bases and df2 are modified base changes window
  colnames(df1) <- c("dna", "pos")
  pos1 <- df1$pos
  pos2 <- df2$pos
  # Identify the missing positions in df1 relative to df2
  missing_pos_in_df1 <- setdiff(pos2, pos1)
  # Add missing rows to df1
  if (length(missing_pos_in_df1) > 0) {
    missing_rows_in_df1 <- df2[df2$pos %in% missing_pos_in_df1, ]
    df1 <- BiocGenerics::rbind(df1, missing_rows_in_df1)
  }
  df1 <- df1[order(df1$pos), ]
  # Return the updated data frames as a list
  return(df1)
}


get_codon_mod_info <- function(out_dna, out_dna_mod){
  out_codon_base_chage <- NULL
  base_mod <- out_dna[,1] == out_dna_mod[,1]
  if (FALSE %in% base_mod) {
    out_codon_base_chage <- paste(out_dna_mod[,1], collapse = "")
  }else{
    out_codon_base_chage <- NA
  }
  return(list(out_codon_base_change = out_codon_base_chage))
}


get_amino_acid_changes_forward <- function(match_positions, 
                                           window_start, 
                                           window_end,
                                           df_exon_info, 
                                           df_exon_info_tmp, 
                                           base_current, 
                                           base_mod,
                                           flank){
  # base_current -> base to search within the window
  # base_mod -> base to modify to once is found in base_current within window
  
  df_out_window <- data.frame()
  for (i in 1:length(match_positions)){
    
    out_window <- NULL
    current_match <- match_positions[[i]]
    current_match$start <- match_positions[[i]]$start + window_start - 1
    current_match$end <- match_positions[[i]]$start + (window_end) - 1
    df_current_match <- df_exon_info_tmp %>% 
                            dplyr::filter(pos >= current_match$start &
                                          pos <= current_match$end)
    if (nrow(df_current_match) > 0) {
      out_window <- get_aa_modification_info(df_current_match, 
                                             df_exon_info,
                                             base_current, 
                                             base_mod)
      out_window$start <- min(df_current_match$pos) - flank
      out_window$end <- max(df_current_match$pos)  - flank
      out_window <- as.data.frame(out_window)
      df_out_window <- BiocGenerics::rbind(df_out_window, out_window )
    }else{
      out_window$out_dna_window <- NA
      out_window$out_codon_base_change <- NA
      out_window$out_aa_info <- NA
      out_window$start <- NA
      out_window$end <- NA
      out_window <- as.data.frame(out_window)
      df_out_window <- BiocGenerics::rbind(df_out_window, out_window )
    }
  }
  df_out_window$base_change <- paste(base_current, "->", base_mod, sep = "")
  return(df_out_window)
} 


get_amino_acid_changes_reverse <- function(match_positions, 
                                           window_start, 
                                           window_end,
                                           df_exon_info, 
                                           df_exon_info_tmp, 
                                           base_current, 
                                           base_mod,
                                           phase_rev,
                                           flank){
  df_out_window <- data.frame()
  for (i in 1:length(match_positions)){
    # convert the df_exon_info to the forward strand
    # This step is needed because it cmes from get_dataframe_sequence_exon_phase
    # this function contains information from the reverse complement exons
    out_window <- NULL
    current_match <- match_positions[[i]]
    current_match$start <- match_positions[[i]]$start + window_start - 1
    current_match$end <- match_positions[[i]]$start + (window_end) - 1
    # If at the end and breaks phase
    df_current_match <- df_exon_info_tmp %>% 
                dplyr::filter(pos >= current_match$start &
                              pos <= current_match$end)
    df_current_match_rev <- df_current_match
    if (nrow(df_current_match_rev) > 0) {
          df_current_match_rev$dna <- as.vector(
            Biostrings::reverseComplement(
              Biostrings::DNAString(paste(df_current_match_rev$dna, collapse = ""))
            ))
          if(phase_rev$start.phase == 0){
            s.match <- (max(df_exon_info$pos) - min(df_current_match_rev$pos) - 
                             nrow(df_current_match_rev) + flank) + 1
          }else if(phase_rev$start.phase == 1){
            s.match <- (max(df_exon_info$pos) - min(df_current_match_rev$pos) - 
                          nrow(df_current_match_rev) + flank)
          }else if(phase_rev$start.phase == 2){
            s.match <- (max(df_exon_info$pos) - min(df_current_match_rev$pos) - 
                          nrow(df_current_match_rev) + flank) - 1
          }
          to_ <- s.match + nrow(df_current_match_rev)
          df_current_match_rev$pos <- seq(
                          from = s.match + 1,
                          to = to_, 1)
          out_window <- get_aa_modification_info(df_current_match_rev, 
                                                 df_exon_info,
                                                 base_current, 
                                                 base_mod)
          
          out_window$start <- min(df_current_match_rev$pos) - flank
          out_window$end <- max(df_current_match_rev$pos)  - flank
          out_window$base_change <- paste(base_current, "->", base_mod, sep = "")
          out_window <- as.data.frame(out_window)
          df_out_window <- BiocGenerics::rbind(df_out_window, out_window )
    }else{
          out_window$out_dna_window <- NA
          out_window$out_codon_base_change <- NA
          out_window$out_aa_info <- NA
          out_window$start <- NA
          out_window$end <- NA
          out_window$base_change <- paste(base_current, "->", base_mod, sep = "")
          out_window <- as.data.frame(out_window)
          df_out_window <- BiocGenerics::rbind(df_out_window, out_window )
    }
  }
  return(df_out_window)
}

