### Extracting information for PyMol ###

library(tidyverse)

my_proteins = read_rds("result_my_proteins")


signal_sequences_numbers = my_proteins %>% 
  dplyr::group_by(protein_ID, kind) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(count))

