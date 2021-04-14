####### Find Orthologous Exons ######
# Fiona adapted from Kathleen's code

# this code takes a folder with the blast outputs converted to CSV files (see reciprocal_blast_copy.sh) 
# and matches the exons in a data frame that has one orthologous exon per row and each species is a column
# It tests for reciprocity (ie reciprocal blast), 
# depends on the assumption of transitivity in the matchings
# (ie that if species 1 exon1 <--> species 2 exon1, and species 2 exon1 <--> species 3 exon1, then sp1 exon1 <--> sp3 exon1)
#note that since this includes self blast. the spreadsheet should include all exons from all species even if they dont match anywhere

# Also, coordinates are with respect to each individual titin sequence, so to get the actual sequences back, I should go to the ttn sequences to subset 

# Set working directory
setwd("/SET/THIS/WORKING/DIRECTORY/")
#setwd("/Users/fionacallahan/Documents/Senior/Thesis/FULL_PIPELINE\ copy/")

species_names_df = read.delim("speciesL.txt", sep="\n", header = FALSE, as.is = T)
species_names = species_names_df$V1

setwd("./nt_matches_csv")

numSpecies = length(species_names)
######### Initialize Lists and Data Frames #########
#initializa data frame that will hold exon data
exon_data <- data.frame(matrix(nrow=200, ncol = 1+numSpecies))
colnames(exon_data) <- c("exon_number", species_names)
rowsUsed = 0 # what row should I put the data in?

# Loop through species names and find all corresponing csv files
for (i in 1:length(species_names)) {
  for (j in 1:length(species_names)) {
    #i=1
    #j=2
    ######### Get species data ##########
    species1 <- species_names[i]
    species2 <- species_names[j]
    
    # get names of reciprocal files
    sp1andsp2id <- paste(species1,"_and_",species2,"_RB.csv",sep="")
    sp2andsp1id <- paste(species2,"_and_",species1,"_RB.csv",sep="")
    
    # read in reciprocal blast data
    sp1andsp2DF <- read.csv(file = sp1andsp2id, header = F)
    colnames(sp1andsp2DF) <- c("query_exon","subject_exon","%_identity","alignment_length","mismatch","gap_openings","query_start","query_end","subject_start","subject_end","evalue","bit_score")
    sp2andsp1DF <- read.csv(file = sp2andsp1id, header = F)
    colnames(sp2andsp1DF) <- c("query_exon","subject_exon","%_identity","alignment_length","mismatch","gap_openings","query_start","query_end","subject_start","subject_end","evalue","bit_score")
    
    # Look at every gene pair in sp1andsp2DF
    for (n in 1:nrow(sp1andsp2DF)) {
      for (m in 1:nrow(sp2andsp1DF)){
        ##### Get info on exon starts/ends  for sp1 blast
        exon_name <- as.character(sp1andsp2DF[n,1])
        split <- strsplit(exon_name, "_")
        sp1_location1 <- split[[1]][3]
        #temp <- strsplit(location, ":")
        #exon1_start1 <- as.numeric(temp[[1]][1])
        #exon1_end1 <- as.numeric(temp[[1]][2])
        ##### Get info on exon2 starts/ends
        exon_name <- as.character(sp1andsp2DF[n,2])
        split <- strsplit(exon_name, "_")
        sp2_location1 <- split[[1]][3]
        #temp <- strsplit(location, ":")
        #exon2_start2 <- as.numeric(temp[[1]][1])
        #exon2_end2 <- as.numeric(temp[[1]][2])
        # GET DATA IN OTHER DIRECTION  (sp2-->sp1 blast)
        exon_name <- as.character(sp2andsp1DF[m,1])
        split <- strsplit(exon_name, "_")
        sp2_location2 <- split[[1]][3]
        exon_name <- as.character(sp2andsp1DF[m,2])
        split <- strsplit(exon_name, "_")
        sp1_location2 <- split[[1]][3]
        
        # if exon1 was the best blast in the sp2andsp1 direction, add to list
        if (sp1_location1 == sp1_location2 & sp2_location1 == sp2_location2){
          #check if either exon is already in the spreadsheet in the columns for each species
          if(sp1_location1 %in% exon_data[,i+1]){
            # find what row to put it in
            exonRow = match(sp1_location1,exon_data[,i+1])
            # put sp2 exon into the same row 
            exon_data[exonRow,j+1] = sp2_location1
          } else if (sp2_location1 %in% exon_data[,j+1]) {
            exonRow = match(sp2_location1,exon_data[,j+1])
            # put sp1 exon into the same row 
            exon_data[exonRow,i+1] = sp1_location1
          } else {
            # Increment rowsUsed and put the data in the next available row
            rowsUsed = rowsUsed + 1
            exon_data[rowsUsed,j+1] = sp2_location1
            exon_data[rowsUsed,i+1] = sp1_location1
          } 
        }#end if block
      }
    } # end for each row of DF
    
  }
}


write.csv(exon_data, file = paste("matched_exon_data.csv",sep=''))

