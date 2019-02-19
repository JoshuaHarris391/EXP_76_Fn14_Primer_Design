#' # [EXP_76](https://github.com/JoshuaHarris391/Lab_Notebook_Cunliffe_2018/blob/master/Experimental_Log/Logs/EXP_76.md) obtaining primer query sequences
#' Primers flanking the two sgRNA target sites from EXP_73 will be determined by first obtaining a query sequence for each target site to put through to primer BLAST
#' 

#' ## Installing package
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("seqinr")

#' ## loading Packages
#+ message=FALSE
library('Biostrings')
library('seqinr')
# library('rstudioapi')
library('tidyverse')


#' ## Setting WD
# setwd('~/Dropbox/Research/PhD/Experiments/EXP_76/EXP_76_Fn14_Primer_Design/scripts/')


#' ## Reading in fasta files
#' These sequences contain 1kB upstream and 1kb downstream of the first and last exon
#'
#' #### Whole Sequence (including 1kB up and downstream)
Fn14_whole_seq_fasta <- read.fasta(file = "../sequence/TNFRSF12A_Whole_Seq.fa", 
                             seqtype = 'DNA',
                             as.string = TRUE)
Fn14_whole_seq_fasta <- Fn14_whole_seq_fasta[[1]] 
Fn14_whole_seq_fasta %>% print()

#' #### Creating character string
Fn14_whole_seq_char <- Fn14_whole_seq_fasta %>% as.character()
Fn14_whole_seq_char %>% class() %>% print()

#' #### Converting to DNAstring
Fn14_whole_seq_dnastring <- DNAString(Fn14_whole_seq_char)
Fn14_whole_seq_dnastring %>% class() %>% print()



#' ## Sequences by Region
#' Note, first and last region are not part of the Fn14 gene
Fn14_Region_seq <- read.fasta(file = "../sequence/TNFRSF12A_Region_Seq.fa", 
                             seqtype = 'DNA',
                             as.string = TRUE)

for (i in 1:length(Fn14_Region_seq)) {
  Fn14_Region_seq[[i]] %>% head() %>% print() 
}





#' ## Loading in sgRNA target sequences
#' Sequences are 5' to 3'
sgRNA_1 <- 'CGGGCGCAGGACGTGCACTA'
sgRNA_2 <- 'AGCTTGGCTCCCGCCGCGTC' 
#' *sgRNA_2 targets the 3' UTR on the antisense strand, therefore it needs to be converted to the reverse complement
#'
#' ## Converting sgRNA_2 to rev comp
sgRNA_2 <- DNAString(sgRNA_2) %>% 
  reverseComplement() %>% 
  as.character()
print(sgRNA_2)


#' #### Match position for sgRNA_1 to Fn14 whole sequence
sgRNA_1_match <- matchPattern(sgRNA_1, Fn14_whole_seq_dnastring)
print(sgRNA_1_match)

#' #### Match position for sgRNA_2 to Fn14 whole sequence
sgRNA_2_match <- matchPattern(as.character(sgRNA_2), Fn14_whole_seq_dnastring)
print(sgRNA_2_match)

#' # Printing query sequence and printing match sequence 
#' ## sgRNA_1 sequence
sgRNA_1 %>% 
  as.character() %>% 
  print()
#' #### Fn14 whole sequence lookup
Fn14_whole_seq_dnastring[start(sgRNA_1_match) : end(sgRNA_1_match)] %>% 
  as.character() %>% 
  print()
#' #### Testing for a match
Fn14_whole_seq_dnastring[start(sgRNA_1_match) : end(sgRNA_1_match)] %>% as.character() %in% sgRNA_1 %>% as.character()


#' ## sgRNA_2 sequence
sgRNA_2 %>% 
  as.character() %>% 
  print()
#' #### Fn14 whole sequence lookup
Fn14_whole_seq_dnastring[start(sgRNA_2_match) : end(sgRNA_2_match)] %>% 
  as.character() %>% 
  print()
#' #### Testing for a match
Fn14_whole_seq_dnastring[start(sgRNA_2_match) : end(sgRNA_2_match)] %>% as.character() %in% sgRNA_2 %>% as.character()




#' # Creating a query sequence for Fn14 primer blast
#' Query sequences will have 'buffer_bp' number of base pairs upstream of the 5' end and 'buffer_bp' number of base pairs downstream from the 3' end of the matching sgRNA sequence
#'
# Defining number of buffer base pairs
buffer_bp <- 800
#' ## sgRNA_1 Primer Query Sequence
EXP_76_sgRNA_1_Query <- Fn14_whole_seq_dnastring[(start(sgRNA_1_match)-buffer_bp) : end(sgRNA_1_match+buffer_bp)] %>% 
  as.character()
print(EXP_76_sgRNA_1_Query)
print(nchar(EXP_76_sgRNA_1_Query))


#' ## sgRNA_2 Primer Query Sequence
EXP_76_sgRNA_2_Query <- Fn14_whole_seq_dnastring[(start(sgRNA_2_match)-buffer_bp) : end(sgRNA_2_match+buffer_bp)] %>% 
  as.character()
print(EXP_76_sgRNA_2_Query)
print(nchar(EXP_76_sgRNA_2_Query))


#' ## CRISPR-CAS9 clevage sites
#' Base pair position is relative to the query sequence 
sgRNA_1_match_query <- matchPattern(sgRNA_1, EXP_76_sgRNA_1_Query)
paste("EXP_76_sgRNA_1_Query Clevage site is between base pairs:", end(sgRNA_1_match_query), "and", end(sgRNA_1_match_query)+1)


sgRNA_2_match_query <- matchPattern(sgRNA_2, EXP_76_sgRNA_2_Query)
paste("EXP_76_sgRNA_2_Query Clevage site is between base pairs:", end(sgRNA_2_match_query), "and", end(sgRNA_2_match_query)+1)


#' ## Exporting query sequences as fasta files
write.fasta(EXP_76_sgRNA_1_Query, "EXP_76_sgRNA_1_Query", "../output/EXP_76_sgRNA_1_Query.fa", as.string = TRUE)
write.fasta(EXP_76_sgRNA_2_Query, "EXP_76_sgRNA_2_Query", "../output/EXP_76_sgRNA_2_Query.fa", as.string = TRUE)



