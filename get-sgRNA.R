

# Extract the sequences where the match was found and 20bp before the match
#' Title
#'
#' @param dna_sequence 
#' @param match_positions 
#' @param sgrna_len 
#'
#' @return
#' @export
#'
#' @examples
extract_sequences <- function(dna_sequence, match_positions, sgrna_len){
  sequences <- character(nrow(match_positions))
  for (i in seq_len(nrow(match_positions))) {
    start_pos <- match_positions[i, 1] - sgrna_len
    if (start_pos < 1) start_pos <- 1
    end_pos <- match_positions[i, 2]
    sequences[i] <- substr(dna_sequence, start_pos, end_pos)
  }
  return(sequences)
}


#' Title
#'
#' @param sub_seq 
#' @param sgrna_len 
#'
#' @return
#' @export
#'
#' @examples
get_info <- function(sub_seq, sgrna_len, window_start, window_end){
  # From a match extract pam, sgRNA based on the 20 bps for example and window
  pam <- substr(sub_seq, 
                    nchar(sub_seq) - 2, 
                    nchar(sub_seq))
  sgRNA <- substr(sub_seq, 
                    nchar(sub_seq) - (sgrna_len + 2), 
                    nchar(sub_seq) - 3)
  window <- substr(sub_seq, 
                    nchar(sub_seq) - (sgrna_len + 2) + window_start - 1, 
                    nchar(sub_seq) - (sgrna_len + 2) + window_end - 1)
  
  return(list(pam = pam, sgRNA = sgRNA, window = window))
}


#' Title
#'
#' @param pattern 
#' @param dna_sequence 
#' @param strand.direction 
#' @param sgrna_len 
#'
#' @return
#' @export
#'
#' @examples
get_pam_matches_info <- function(pattern,
                              dna_sequence,
                              strand.direction,
                              sgrna_len,
                              window_start,
                              window_end)
{
  # Search for PAMs for a given dna sequence
  match_positions <- lapply(pattern, function(x) 
    as.data.frame(stringr::str_locate_all(as.character(dna_sequence), x)[[1]]))
  match_positions <- dplyr::bind_rows(match_positions)
  match_positions <- match_positions[order(match_positions[,1]),]
  match_positions <- match_positions %>%
    filter(start > sgrna_len)
  extracted_sequences <- extract_sequences(as.character(dna_sequence), 
                                           match_positions, sgrna_len)
  out <- lapply(extracted_sequences, get_info, sgrna_len, window_start, window_end)
  out <- as.data.frame(do.call(cbind, out))
  out <- as.data.frame(t(as.data.frame(do.call(cbind, out))), row.names = F)
  out$strand <- strand.direction
  return(out)
}


#' Title
#'
#' @param df.target.trans 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
get_dna_trasncript_sequence  <- function(df.target.trans, Genome)
{
  # Extract all exon sequences for a given transcript
  match_seqs <- Biostrings::getSeq(Genome,
                       df.target.trans$chr,
                       strand = df.target.trans$strand,
                       start = df.target.trans$start,
                       end = df.target.trans$end)
  match_seqs <- as.vector(unlist(lapply(match_seqs, as.character)))
  trans_sequence <- paste(match_seqs, collapse = "")
  return(trans_sequence)
}
