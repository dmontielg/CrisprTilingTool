

#' Title
#'
#' @param input_string 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
repeat_characters <- function(input_string, n) {
  # Split the string into individual characters
  characters <- strsplit(input_string, "")[[1]]
  # Repeat each character 'n' times
  repeated_chars <- sapply(characters, strrep, n)
  # Combine the repeated characters into a single string
  result <- paste0(repeated_chars, collapse = "")
  return(result)
}

#------------------------------------------------------------------------------$

#' Title
#'
#' @param df.target.trans 
#'
#' @return
#' @export
#'
#' @examples
get_phase_transcript <- function(df.target.trans, Genome){
  df.obj.target.trans.seq <- Biostrings::getSeq(Genome, 
                                                df.target.trans$chr, 
                                                start = df.target.trans$start, 
                                                end = df.target.trans$end,
                                                strand = df.target.trans$strand)
  df.exons <- data.frame(seq = character(), width = numeric(), 
                         start = numeric(), end = numeric(),
                         start.phase = numeric(), end.phase = numeric())
  first.exon <- TRUE # deactive this flag after checking the first exon
  end.phase.flag <- FALSE
  start.phase.flag <- FALSE
  for (i in 1:length(df.obj.target.trans.seq)){
    var.seq <- as.character(as.character(df.obj.target.trans.seq[i,]))
    var.width <- as.numeric(df.obj.target.trans.seq[i,]@ranges@width)
    var.start <- as.numeric(df.target.trans[i,]$start)
    var.end <- as.numeric(df.target.trans[i,]$end)
    var.exon <- as.numeric(df.target.trans[i,]$exon)
    # This only works on the first exon verify the length and see if
    #var.width %% 3 -> 2, There are two bp and one is missing
    #var.width %% 3 -> 1, There is one bp and two are missing
    #var.width %% 3 -> 0, There are no interruptions
    if (first.exon) {
      first.exon <- FALSE
      end.phase.flag <- var.width %% 3
      start.phase.flag <- 0 # start codon Amino Acid M in the transcript
    }else{
      if (end.phase.flag == 0){
        start.phase.flag <- 0
      }else if(end.phase.flag == 1){
        start.phase.flag <- 1
      }else if(end.phase.flag == 2){
        start.phase.flag <- 2
      }
      if (start.phase.flag == 0){
        end.phase.flag <- var.width %% 3
      }else if(start.phase.flag == 1){
        end.phase.flag <- (var.width - 2) %% 3
      }else if(start.phase.flag == 2){
        end.phase.flag <- (var.width - 1) %% 3
      }
    }
    new_row <- data.frame(seq = var.seq, width = var.width,
                          start = var.start, end = var.end, 
                          start.phase = start.phase.flag,
                          end.phase = end.phase.flag,
                          exon = var.exon)
    df.exons <- BiocGenerics::rbind(df.exons, new_row)
  }
  return(df.exons)
}


#' Title
#'
#' @param df.phase.exons.trans 
#'
#' @return
#' @export
#'
#' @examples
complete_phase_sequences <- function(df.phase.exons.trans){
  ## Function to complete the bp's based if phase is at the start or end of exon
  complete.seq.list <- list()
  for (i in 1:nrow(df.phase.exons.trans)) {
    current.phase.exon <- df.phase.exons.trans[i,]
    complete.seq <- ""
    start.seq <- ""
    end.seq <- ""
    
    if (current.phase.exon$start.phase == 0) {
      # nothing to do first codon is complete
      start.seq <- ""  
    }
    else if (current.phase.exon$start.phase == 1) {
      # there is a phase 1 at start, therefore 1 bp is missing
      start.seq <- stringr::str_sub(df.phase.exons.trans[i-1,]$seq, -1)
    }
    else if (current.phase.exon$start.phase == 2) {
      # there is a phase 2 at start, therefore 2 bp's are missing
      start.seq <- stringr::str_sub(df.phase.exons.trans[i-1,]$seq, -2)
    }
    if (current.phase.exon$end.phase == 0) {
      end.seq <- ""
    }
    else if (current.phase.exon$end.phase == 1) {
      end.seq <- stringr::str_sub(df.phase.exons.trans[i+1,]$seq,  
                                  start = 1, end = 2) 
    }
    else if (current.phase.exon$end.phase == 2) {
      end.seq <- stringr::str_sub(df.phase.exons.trans[i+1,]$seq,  
                                  start = 1, end = 1) 
    }
    complete.seq.list[i] <- paste(start.seq, current.phase.exon$seq, end.seq, 
                                  sep ="")
  }
  
  return(complete.seq.list)
}


get_exon_phase_sequences <- function(df.phase.exons.trans){
  ## Function to complete the bp's based if phase is at the start or end of exon
  complete.seq.list <- list()
  for (i in 1:nrow(df.phase.exons.trans)) {
    current.phase.exon <- df.phase.exons.trans[i,]
    exon.seq <- ""
    start.seq <- ""
    end.seq <- ""
    if (current.phase.exon$start.phase == 0) {
      # nothing to do first codon is complete
      exon.seq <- stringr::str_sub(df.phase.exons.trans[i,]$seq)
    }
    else if (current.phase.exon$start.phase == 1) {
      # there is a phase 1 at start, therefore 1 bp is missing
      exon.seq <- paste0(stringr::str_sub(df.phase.exons.trans[i-1,]$seq, -1),
                         stringr::str_sub(df.phase.exons.trans[i,]$seq))
    }
    else if (current.phase.exon$start.phase == 2) {
      # there is a phase 2 at start, therefore 2 bp's are missing
      exon.seq <- stringr::str_sub(df.phase.exons.trans[i,]$seq, start = 2)
    }
    if (current.phase.exon$end.phase == 1) {
      exon.seq <- stringr::str_sub(exon.seq, 
                                   start = 1, end = -2)
    }
    else if (current.phase.exon$end.phase == 2) {
      exon.seq <- paste0(exon.seq, stringr::str_sub(
                                    df.phase.exons.trans[i+1,]$seq,  
                                    start = 1, end = 1))
    }
    complete.seq.list[i] <- exon.seq
  }
  return(complete.seq.list)
}


assign_phase_reverse <- function(df.phase.exons.trans){
  
  # We are ging to reverse the phase for each exon for the reverse strand
  # If phase 1 at the start it becomes phase 2
  # If phase 2 at the start it becomes phase 1 same for phase in last codon
  # Temp variables for stand and end of the exon
  tmp.s.p <- c()
  tmp.e.p <- c()
  
  for (i in 1:nrow(df.phase.exons.trans)) {
    if(df.phase.exons.trans[i,]$start.phase == 0){
      tmp.s.p[i] <- 0
    }else if(df.phase.exons.trans[i,]$start.phase == 1){
      tmp.s.p[i] <- 2
    }else if(df.phase.exons.trans[i,]$start.phase == 2){
      tmp.s.p[i] <- 1
    }
    if(df.phase.exons.trans[i,]$end.phase == 0){
      tmp.e.p[i] <- 0
    }else if(df.phase.exons.trans[i,]$end.phase == 1){
      tmp.e.p[i] <- 2
    }else if(df.phase.exons.trans[i,]$end.phase == 2){
      tmp.e.p[i] <- 1
    }
  }
  df.phase.exons.trans.rev <- data.frame(
    seq = unlist(lapply(lapply(lapply(
      df.phase.exons.trans$seq, Biostrings::DNAString), 
      Biostrings::reverseComplement), as.character)),
    width = df.phase.exons.trans$width,
    start = df.phase.exons.trans$start,
    end = df.phase.exons.trans$end,
    start.phase = tmp.e.p,
    end.phase = tmp.s.p,
    exon = df.phase.exons.trans$exon,
    complete = unlist(lapply(lapply(lapply(
      df.phase.exons.trans$complete, Biostrings::DNAString), 
      Biostrings::reverseComplement), as.character)))
  df.phase.exons.trans <- df.phase.exons.trans.rev
  
  return(df.phase.exons.trans)
}


get_dataframe_sequence_exon_phase <- function(df.phase.exons.trans, 
                                              exon_n, 
                                              flank){
  ## Convert the phase and phase.complete into dataframes with the relative
  # and absolute sequence position within the exon
  current.phase <- df.phase.exons.trans %>% dplyr::filter(exon == exon_n)
  characters_exon <- Biostrings::strsplit(repeat_characters(as.character(
    current.phase$complete), 1), "")[[1]]
  df_exon_info <- NULL
  df_exon_info_tmp <- NULL
  if (current.phase$start.phase == 0) {
    df_exon_info <- as.data.frame(cbind(characters_exon, 
                                        1:length(characters_exon)))
    colnames(df_exon_info) <- c("dna","pos")
    df_exon_info$pos <- as.numeric(df_exon_info$pos)
    df_exon_info$pos <- df_exon_info$pos + flank
    df_exon_info_tmp <- df_exon_info
    
  }else if (current.phase$start.phase == 1) {
    df_exon_info <- as.data.frame(cbind(characters_exon, 
                                        1:length(characters_exon)))
    colnames(df_exon_info) <- c("dna","pos")
    df_exon_info$pos <- as.numeric(df_exon_info$pos)
    df_exon_info$pos <- df_exon_info$pos + flank - 1
    df_exon_info_tmp <- df_exon_info %>% 
                          dplyr::filter(pos >= min(pos) + 1)
    
  }else if (current.phase$start.phase == 2) {
    df_exon_info <- as.data.frame(cbind(characters_exon, 
                                        1:length(characters_exon)))
    colnames(df_exon_info) <- c("dna","pos")
    df_exon_info$pos <- as.numeric(df_exon_info$pos)
    df_exon_info$pos <- df_exon_info$pos + flank - 2
    df_exon_info_tmp <- df_exon_info %>% 
                          dplyr::filter(pos >= min(pos) + 2)
  }
  if (current.phase$end.phase == 1) {
    df_exon_info_tmp <- df_exon_info_tmp %>% 
                          dplyr::filter(pos <= max(pos) - 2)
  }else if (current.phase$end.phase == 2) {
    df_exon_info_tmp <- df_exon_info_tmp %>% 
                          dplyr::filter(pos <= max(pos) - 1)
  }
  return(list(df_exon_info, df_exon_info_tmp))
}


get_absolute_coordinates_aa <- function(df.phase.exons.trans){
  # input dataframe of all the exons of a given transcrript
  df.aa.exons <- list()
  df.aa.exons$seq <- unlist(lapply( 
    lapply(
      lapply(
        #df.phase.exons.trans$exon.seq,
        df.phase.exons.trans$complete, 
        Biostrings::DNAString), 
      Biostrings::translate), 
    as.character))
  df.aa.exons$exon <- as.character(df.phase.exons.trans$exon)
  df.aa.exons <- as.data.frame(df.aa.exons)
  df.aa.exons.index <- list()
  for (i in 1:nrow(df.aa.exons)) {
    df.tmp <- list()
    characters <- unlist(strsplit(df.aa.exons$seq[i], ""))
    df.tmp$aa <- characters
    df.tmp$exon <- rep(df.aa.exons$exon[i], length(characters))
    df.tmp$i <- seq(1, length(characters))
    df.tmp <- as.data.frame(df.tmp)
    df.aa.exons.index <- rbind(df.aa.exons.index, df.tmp)
  }
  
  df.aa.exons.index$abs <- 1:nrow(df.aa.exons.index)
  flag.phase <- F
  list.pos <- c()
  for (i in 1:nrow(df.phase.exons.trans)) {
    current.exon <- df.phase.exons.trans[i,]$exon
    current.phase <- df.phase.exons.trans %>% filter(exon == current.exon)
    if (i > 1) {
        previous.end.phase <- df.phase.exons.trans[i-1,]$end.phase
        if (previous.end.phase > 0) {
          previous.exon <- df.phase.exons.trans[i-1,]$exon
          previous.df.aa.exons <- df.aa.exons.index %>% dplyr::filter(exon == previous.exon)
          previous.pos <- tail(list.pos, n = 1)
          df.aa.exons.index.update <- df.aa.exons.index %>% dplyr::filter(exon == current.exon)
          pos.update <- 0:(nrow(df.aa.exons.index.update) - 1)
          pos.update <- pos.update + previous.pos
          length(pos.update)
          list.pos <- c(list.pos, pos.update)
      }

      else if (previous.end.phase == 0){
        previous.exon <- df.phase.exons.trans[i-1,]$exon
        previous.df.aa.exons <- df.aa.exons.index %>% dplyr::filter(exon == previous.exon)
        previous.pos <- tail(list.pos, n = 1) + 1
        df.aa.exons.index.update <- df.aa.exons.index %>% dplyr::filter(exon == current.exon)
        pos.update <- 0:(nrow(df.aa.exons.index.update) - 1)
        pos.update <- pos.update + previous.pos
        length(pos.update)
        list.pos <- c(list.pos, pos.update)
      }
    }else{
      df.aa.current <- df.aa.exons.index %>% filter(exon == current.exon)
      list.pos <- c(list.pos, df.aa.current$abs)
    }
  }
  df.aa.exons.index$abs <- list.pos
  return(df.aa.exons.index)
}


split_aa_info <- function(input_string) {
  # Regular expression pattern: a character, followed by one or more digits, 
  # followed by another character also works for asteriscs
  #pattern <- "([a-zA-Z\\*])(\\d+)([a-zA-Z\\*])"
  pattern <- "(Ter|[a-zA-Z])(\\d+)(Ter|[a-zA-Z])"
  
  matches <- stringr::str_match_all(input_string, pattern)[[1]]
  # Convert matches to a list of named vectors for easier access
  result <- apply(matches, 1, function(match) {
    list(char1 = match[2], number = as.numeric(match[3]), char2 = match[4])
  })
  return(result)
}


get_map_coordinates <- function(df_table_all, 
                                df.phase.exons.trans, 
                                df.aa.exons.index){
  out_aa_info_abs <- list()
  for (i in 1:nrow(df_table_all)) {
    
    start.phase <- df.phase.exons.trans %>% filter(exon == df_table_all[i,]$exon)
    if (is.na(df_table_all[i,]$out_aa_info)) {
      out_aa_info_abs <- c(out_aa_info_abs, NA)
    }else{
      current_aa_split <- split_aa_info(df_table_all[i,]$out_aa_info)
      str_aa_abs <- ""
      for (j in 1:length(current_aa_split)){
    
        current_aa_split_index <- current_aa_split[j]
        current_aa_split_index <- unlist(current_aa_split_index)
        exon_n <- df_table_all[i,]$exon
        if (start.phase$start.phase > 0) {
          #current_aa_split_index[2] <- as.numeric(current_aa_split_index[2]) - 1
          current_aa_split_index[2] <- as.numeric(current_aa_split_index[2])
        }
        current_abs <- df.aa.exons.index %>% dplyr::filter(
                                                exon == exon_n,
                                                i == current_aa_split_index[2])
        str_aa_abs <- paste(str_aa_abs, 
                            current_aa_split_index[1],
                            current_abs$abs,
                            current_aa_split_index[3],
                            sep = "")
      }
      out_aa_info_abs <- c(out_aa_info_abs, str_aa_abs)
    }
  }
  return(out_aa_info_abs)
}


