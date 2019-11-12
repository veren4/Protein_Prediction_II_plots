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

my_proteins = full_join(my_proteins, locations, by = c("protein_ID" = "protein_ID"))

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

occupancy_NLS = table(occupied_percentage_boxes_NLS) %>% 
  as.data.frame() %>% 
  dplyr::rename(percentage_box = occupied_percentage_boxes_NLS) %>% 
  add_column(kind = "NLS")
occupancy_NES = table(occupied_percentage_boxes_NES) %>% 
  as.data.frame() %>% 
  dplyr::rename(percentage_box = occupied_percentage_boxes_NES) %>% 
  add_column(kind = "NES")

occupancy_combined = bind_rows(occupancy_NLS, occupancy_NES)
occupancy_combined$kind = factor(occupancy_combined$kind, levels = c("NLS", "NES"))

perc_boxes_plot = ggplot(data = occupancy_combined, mapping = aes(x = percentage_box,
                                                                  y = kind)) +
  geom_tile(mapping = aes(fill = Freq),
              stat = "identity") +
  theme(axis.ticks = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(name = "protein sequence [%]", breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  scale_fill_continuous(name = "absolute occupancy")

perc_boxes_plot


# better plot with lines ------------------------------------------------------

# Für diesen Plot könnte ich die Prozentpunkte auch mit mehr Nachkommastellen berechnen.

my_proteins$protein_ID = factor(my_proteins$protein_ID)
my_proteins$kind = factor(my_proteins$kind, levels = c("NLS", "NES"))

gantt = ggplot(data = my_proteins, mapping = aes(x = perc_start, y = protein_ID,
                                                 xend = perc_end, yend = protein_ID)) +
  geom_segment(mapping = aes(color = kind)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank()) +
  ylab("protein sequence") +
  xlab("protein sequence position [%]") +
  scale_color_manual(values = c("orange", "green"),   # What are the actual colors?
                     name = element_blank())

gantt



