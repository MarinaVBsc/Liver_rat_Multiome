library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(biomaRt)
library(magrittr)
library(tidyverse)
library(dplyr)

# Save Four Files from CellChatDB

# 1a. Load CellChatDB Reference Files

CellChatDB <- CellChatDB.mouse

# Pull Mouse Information
geneInfo <- CellChatDB$geneInfo
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor

# Write Mouse Information
write.csv(geneInfo, file = "geneInfo.csv")
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChat.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChat.csv")

# 1b. Load BioMart rat and mouse gene nomenclature

# Find Chx Homologues
rat <- useEnsembl("ensembl", 'rnorvegicus_gene_ensembl', verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
mouse <- useEnsembl("ensembl", 'mmusculus_gene_ensembl', verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
annot_table <- getLDS(
  mart = mouse,
  attributes = c('ensembl_gene_id','mgi_symbol','external_gene_name','chromosome_name'),
  martL = rat,
  attributesL = c('ensembl_gene_id','external_gene_name','chromosome_name','gene_biotype'))
head(annot_table)

# Create Cross Reference Table Based on Ensembl Gene ID
xref_table <- annot_table %>% 
  dplyr::select(Gene.stable.ID, Gene.name, Gene.stable.ID.1, Gene.name.1) %>% 
  dplyr::rename(Ensembl.Gene.ID = Gene.stable.ID)

# 2. Replace Genes

# geneInfo: Symbol column
# there are duplicates in the Gene Info rownames
geneInfo.rat <- geneInfo %>% 
  rownames_to_column("row_names") %>%
  inner_join(xref_table, by = "Ensembl.Gene.ID") %>% 
  rename(Symbol.mouse = Symbol) %>% 
  rename(Symbol.rat = Gene.name.1) %>% 
  mutate(Symbol = Symbol.rat) %>% 
  na_if("") %>% 
  na.omit()

geneInfo.rat <- geneInfo.rat[!duplicated(geneInfo.rat$Symbol),]

# interaction_input: ligand column and receptor column
ligands <- geneInfo.rat %>% 
  dplyr::select(Symbol.mouse, Symbol.rat) %>% 
  rename(ligand = Symbol.mouse, ligand.rat = Symbol.rat)
receptors <- geneInfo.rat %>% 
  dplyr::select(Symbol.mouse, Symbol.rat) %>% 
  rename(receptor = Symbol.mouse, receptor.rat = Symbol.rat)

interaction_input.rat <- interaction_input %>% 
  inner_join(ligands, by = "ligand") %>% 
  inner_join(receptors, by = "receptor") %>%
  rename(ligand.mouse = ligand) %>% 
  rename(receptor.mouse = receptor) %>% 
  rename(ligand = ligand.rat) %>% 
  rename(receptor = receptor.rat)

# complex_input: subunit_1 ... subunit_4
mouse2rat <- geneInfo.rat %>% 
  dplyr::select(Symbol.mouse, Symbol.rat) %>% 
  rename(value = Symbol.mouse) %>% 
  unique()

complex_input.rat <- complex_input %>% 
  rownames_to_column("row_names") %>%
  pivot_longer(cols = subunit_1:subunit_4) %>%
  inner_join(mouse2rat, by = "value") %>% 
  pivot_wider(id_cols = row_names, names_from = name, values_from = Symbol.rat, values_fn = first, values_fill = "") %>%
  column_to_rownames(var = "row_names")

# cofactor_input: cofactor1 ... cofactor16

mouse2rat <- geneInfo.rat %>% dplyr::select(Symbol.mouse, Symbol.rat) %>% rename(value = Symbol.mouse) %>% unique()

cofactor_input.rat <- cofactor_input %>% 
  rownames_to_column("row_names") %>%
  pivot_longer(cols = cofactor1:cofactor16) %>%
  inner_join(mouse2rat, by = "value") %>% 
  pivot_wider(id_cols = row_names, names_from = name, values_from = Symbol.rat, values_fn = first, values_fill = "") %>%
  column_to_rownames(var = "row_names")

# 3. Update the Database
options(stringsAsFactors = FALSE) 
CellChatDB <- list() 
CellChatDB$interaction <- interaction_input.rat 
CellChatDB$complex <- complex_input.rat
CellChatDB$cofactor <- cofactor_input.rat
CellChatDB$geneInfo <- geneInfo.rat
CellChatDB.rat <- CellChatDB

#interaction_input <- read.csv(file = 'interaction_input_CellChat.csv') 
#row.names(interaction_input) <- interaction_input[,1] 
#complex_input <- read.csv(file = 'complex_input_CellChat.csv', row.names = 1) 
#cofactor_input <- read.csv(file = 'cofactor_input_CellChat.csv', row.names = 1) 

# 4. Rebuild (optional) - WIP, but can just use the object.
# devtools::use_data(CellChatDB.rat, overwrite = TRUE)
# usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

showDatabaseCategory(CellChatDB, nrow = 1)