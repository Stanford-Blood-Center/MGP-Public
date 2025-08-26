#script to update alignments file to most recent alignments
#Match Grade Populator © Stanford Blood Center, LLC.
#8/26/2025 - v 1.1.0 


### VERSION HISTORY 
# 8/26/2025 - v 1.1.0 
  #removed dependency on external BLAASD() script; use buildAlignments() from HLAtools instead
  #add alignment version to filename 
# 1/10/2024 - v 1.0.0


library(HLAtools)

alignments<-buildAlignments(c('A', 'B', 'C', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5'), source = 'AA')

version <- strsplit(alignments[[1]]$Version, ' ')[[1]][[2]]

finalAlignments<-sapply(alignments, '[[', 1)

fileName <- sprintf('./ref/alignments_%s.rda', version)

#remove old alignments file 
file.remove(list.files('ref/', 'alignment', full.names = TRUE))

saveRDS(file=fileName, finalAlignments)

