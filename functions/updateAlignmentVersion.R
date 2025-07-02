#script to update alignments file to most recent 
#Match Grade Populator © Stanford Blood Center, LLC.
#1/10/2024 v 1.0

source('functions/BLAASD.R')

alignments<-BLAASD(c('A', 'B', 'C', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5'))

saveRDS(file='./ref/alignments.rda', alignments)
