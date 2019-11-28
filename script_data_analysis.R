#
# density plot of positions of NLSs/ NESs in the protein sequences
#

# in fasta file: find out absolute lengths of proteins

library(tidyverse)

my_protein_sequences = seqinr::read.fasta(file = "C:/Users/Verena/1_Studium/03_Aufbaustudium_Informatik/Protein Prediction II/Exercise/data/ns/nes_nls.fasta",
                                          seqtype = "AA",
                                          as.string = T)

my_proteins = data.frame(protein_ID = character(),
                         aa_sequence = character(),
                         sequence_length = integer(),
                         stringsAsFactors = F)

for(i in 1:length(my_protein_sequences)){
  my_proteins[i, "protein_ID"] = attributes(my_protein_sequences[[i]])$name
  my_proteins[i, "aa_sequence"] = my_protein_sequences[[i]][1]
  my_proteins[i, "sequence_length"] = seqinr::getLength(my_protein_sequences[[i]])
}


# in excel file: extract start + end of NLS/ NES and wether its an NLS or NES

locations = read_tsv(file = "C:/Users/Verena/1_Studium/03_Aufbaustudium_Informatik/Protein Prediction II/Exercise/data/ns/nes_nls.tab",
                     col_names = F) %>% 
  dplyr::rename(protein_ID = X1,
                start = X2,
                end = X3,
                kind = X4)

my_proteins = dplyr::full_join(my_proteins, locations, by = c("protein_ID" = "protein_ID"))

rm(my_protein_sequences, i, locations)

# Theoretisch ist es nicht gründlich genug, davon auszugehen, dass die Sequenzen
# in den beiden Dokumenten gleich orientiert sind.

my_proteins = dplyr::mutate(my_proteins,
                            perc_start = round((start/sequence_length)*100, digits = 0),
                            perc_end = round((end/sequence_length)*100, digits = 0))

occupied_percentage_boxes_NLS = c()
occupied_percentage_boxes_NES = c()

for(i in 1:nrow(my_proteins)){
  temp = my_proteins[i, "perc_start"]:my_proteins[i, "perc_end"]
  
  if(my_proteins[i, "kind"] == "NLS"){
    occupied_percentage_boxes_NLS = c(occupied_percentage_boxes_NLS, temp)
  }else{
    occupied_percentage_boxes_NES = c(occupied_percentage_boxes_NES, temp)
  }
}



# Für diesen Plot könnte ich die Prozentpunkte auch mit mehr Nachkommastellen berechnen.

my_proteins$protein_ID = factor(my_proteins$protein_ID)
my_proteins$kind = factor(my_proteins$kind, levels = c("NLS", "NES"))


# backup = my_proteins
# my_proteins = backup

my_proteins_kinder = my_proteins %>% 
  select(protein_ID, kind) %>% 
  group_by(protein_ID) %>% 
  summarise(signal_cases = length(unique(kind)))

my_proteins = left_join(my_proteins, my_proteins_kinder, by = c("protein_ID" = "protein_ID"))

my_proteins = my_proteins %>% 
  mutate(final_kind_sorter = case_when(signal_cases == 2 ~ "both",
                   kind == "NES" ~ "NES",
                   kind == "NLS" ~ "NLS"))

start_sorter = select(my_proteins, protein_ID, perc_start) %>% 
  group_by(protein_ID) %>% 
  summarise(first_start = min(perc_start))

my_proteins = left_join(my_proteins, start_sorter, by = c("protein_ID" = "protein_ID"))

my_proteins$final_kind_sorter = factor(my_proteins$final_kind_sorter, levels = c("NLS", "both", "NES"))

# sort by first_start within the 3 kind groups --------------------------------

NLS_orderer = dplyr::filter(my_proteins, final_kind_sorter == "NLS") %>% 
  
  arrange(first_start)

########### here lieth the problem!!!!!!!! ##################
# Du musst innerhalb der facet sortieren, und dabei berücksichtigen, dass später
# NxS eines Proteins innerhalb einer Zeile sind.



# ggplot2 ---------------------------------------------------------------------

gantt = ggplot(data = my_proteins, mapping = aes(x = perc_start, y = protein_ID,
                                                 xend = perc_end, yend = protein_ID)) +
  facet_grid(rows = vars(final_kind_sorter)) +
  geom_segment(mapping = aes(color = kind)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=14)) +
  ylab("protein sequence") +
  xlab("relative protein sequence position [%]") +
  scale_color_manual(values = c("#37c837", "#ff7f0e"),
                     name = element_blank())

gantt



# plot separately -------------------------------------------------------------

NLS_plot_input = filter(my_proteins, final_kind_sorter == "NLS")
NLS_plot_input$perc_start = factor(NLS_plot_input$perc_start, levels = order(NLS_plot_input$first_start))

NLS_gantt = ggplot(data = NLS_plot_input, mapping = aes(x = perc_start, y = protein_ID,
                                                 xend = perc_end, yend = protein_ID)) +
  geom_segment(mapping = aes(color = kind)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=14)) +
  ylab("protein sequence") +
  xlab("relative protein sequence position [%]") +
  scale_color_manual(values = c("#37c837", "#ff7f0e"),
                     name = element_blank())

NLS_gantt



