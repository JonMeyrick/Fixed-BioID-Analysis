## Load Packages ####
  library(tidyverse)
  library(dplyr)
  library(readxl)   
  library(openxlsx)
  library(devtools)
  library(data.table)
  library(plyr)
  library(readr)
  library(pheatmap)
  library(limma)
  library(UpSetR)
  library(RColorBrewer)
  library(circlize)
  library(corrplot)
  library(ggplot2)
  library(stringi)


## Significant Interaction CSVs ####
sigint <- function(inputpath, sheetno, trecols){
  
  dir.create("sigint_output")
  setwd("~/Desktop/sigint_output")
  read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
  }
  
  for (sheetno in sheetno:length(names(x))){
    
    redexcel <- read_excel_allsheets(inputpath)
    
    lone_baitname <- as.data.frame(redexcel[sheetno])
    #3 columns needed are: Interactor name, logFC, p value significance
    #bait_important <- lone_baitname[,c(9,10,25)]
    bait_important <- lone_baitname[,trecols]
    #Filtered based on logFC ≥ 1.5 and p-value < 0.05 
    significant_interactions <- bait_important %>% filter(bait_important[,3] > -log10(0.05) & bait_important[,2] >= 1.5)
    sigint_output <<- significant_interactions[,1:2]
    View(sigint_output)
    
    #filename = readline(prompt = "Enter file name: ")
    
    filename = names(x)[sheetno]
    
    write.csv(sigint_output, paste0(filename,"_sigint.csv"), row.names = FALSE)
  }
  setwd("~/Desktop")
}
sigint("~/Desktop/Singh BioID/Singh_NewNomen.xlsx", 2,c(9,10,25))


## Significant Interactions Matrix ####
sigint_mat <- function(path){
  data_all <- list.files(path, pattern = "*.csv", full.names = TRUE)
  
  samples <- as.data.frame(sapply(strsplit(sapply(strsplit(data_all,"/"), function(x) x[7]),"[.]"),function(x) (x)[1]))
  names(data_all) <- samples[,1]
  
  explist <- list()
  for(i in 1:length(data_all)){
    Name <- names(data_all[i])
    explist[[Name]] <- read.csv(data_all[i],header = TRUE, row.names = 1)
  }
  #Create matrix of LogFC between samples
  uniq.Trans=unique(unlist(lapply(explist,function(x) rownames(x))))
  sigint_mat_output=matrix(0,length(uniq.Trans),length(explist))
  rownames(sigint_mat_output)=uniq.Trans
  colnames(sigint_mat_output)=gsub("_sigint","",names(explist),fixed = TRUE)
  for (i in 1:length(explist)){
    sigint_mat_output[rownames(explist[[i]]),i]=explist[[i]][,1]
  }
  sigint_mat_output <<- sigint_mat_output
  #Create binary matrix
  sigint_mat_binary <<- sigint_mat_output == 0
  #Create prey list
  preyList <-rep(rownames(sigint_mat_binary), times = ncol(sigint_mat_binary))
  #Populate matrix with protein names
  for (i in 1:length(sigint_mat_binary)){
    if(sigint_mat_binary[i] == FALSE){sigint_mat_binary[i] <<-preyList[i]}
  }
  
  #Prey LogFC Heatmap
  pdf(file = "sigint_heatmap.pdf", height=25,width=10)
  cols <- colorRampPalette(c("white", "blue"))(100)
  sigint_heatmap <- pheatmap(sigint_mat_output, color =  cols, fontsize = 10, cellwidth = 10, cellheight = 10)
  sigint_heatmap
  dev.off()
}
sigint_mat("~/Desktop/Singh BioID/Sigint_output")


## Key Site Subsetting ####
sigint_subset <- function(a){
  sigint_SS <<- as.data.frame(sigint_mat_output[,a])
  sigint_SS_binary <<- as.data.frame(sigint_mat_binary[,a])
  SS_list <- list() #Unused currently but need to merge into list for further functions
  for (i in 1:ncol(sigint_SS_binary)){
    assign(colnames(sigint_SS_binary[i]),sigint_SS_binary[i][sigint_SS_binary[i] != TRUE], envir = parent.frame())
  }
}
sigint_subset(c(8,9,21,22))



## Upset Plot ####
sigint_plot <- function(filename){

  barcols <- rand_color(50)
  setsnumb <- sum(seq(ncol(sigint_SS):1))
  
  sigint_SS[sigint_SS > 0] <- 1 #Must be in a binary format

  pdf(file = paste("tester",".pdf",sep = ""), height=25, width=60)
  upset(sigint_SS, sets = colnames(sigint_SS), 
        
        text.scale = 10, point.size = 15, line.size =5, 
        sets.bar.color = c(barcols[1:ncol(sigint_SS)]),
        main.bar.color = c(barcols[1:setsnumb]),
        #matrix.color = c(barcols[1:setsnumb]),
        sets.x.label = "Prey Observation No.", mainbar.y.label = "Number of prey")
  dev.off()
}
sigint_plot("test_upset")
  


## MRP Removal (Needs Automating/Unnecessary?) ####
#Automation: create list of human MRPs, remove from each list if matching (may be pointless)
#MRPlist <- read.xlsx("~/Desktop/MRP Nomenclature/MRPnomenclature.xlsx")
#MRPlist <- gsub("","",as.character(MRPlist), fixed = TRUE) #doesnt work with trailing backslash

#example automated MRP removal using predefined character string
#mL57rem <- c("MRPL3","MRPL15")
# mL57[is.na(match(mL57,MRPlist))] #doesn't work

mL57 <- sort(mL57[c(-13:-14)])
mL67 <- sort(mL67[c(-8:-9)])
uL23m <- sort(uL23m[c(-7,-17:-18,-33,-35,-38:-39,-41,-43,-49,-60:-62)])
uL29m <- sort(uL29m[c(-18)])


## Coverage Analysis (Needs Automating) ####
#Determine no. of unique prey identified
SS_TotalPrey <- unique(c(mL57,mL67,uL23m,uL29m))
length(SS_TotalPrey)

#1st bait is bait with highest no. of prey: uL23m

SS_mL57v1stBait <- setdiff(mL57,uL23m)
SS_uL29mv1stBait <- setdiff(uL29m,uL23m)
SS_mL67v1stBait <- setdiff(mL67,uL23m)

SS_unique_baits1 <- unique(c(uL23m,uL29m))
#2nd bait is highest no. of unique prey vs selected bait: uL29m

SS_mL67v2ndBait <- setdiff(mL67,SS_unique_baits1)
SS_mL57v2ndBait <- setdiff(mL57,SS_unique_baits1)

SS_unique_baits2 <- unique(c(uL23m,uL29m,mL57))
#3rd bait is highest no. of unique prey vs selected baits: mL57
#4th bait is last remaning bait: mL67

#Compile interactor lists
sigint_SS2 <- list(mL57, mL67, uL23m, uL29m)




## 'Habitat' (Overlap:Distance) Analysis ####
#Order: mL57, mL67, uL23m, uL29m
sigint_habitat <- function (filename,distanceinputpath){

matMtx <- matrix(nrow = length(sigint_SS2), ncol = length(sigint_SS2))
colnames(matMtx) <- colnames(sigint_SS)
rownames(matMtx) <- colnames(sigint_SS)

for (pos in 1:ncol(matMtx)){ 
  for (rowpos in 1:ncol(matMtx)){
  matches <- sum(na.rm = TRUE, sigint_SS2[[pos]] %in% sigint_SS2[[rowpos]])
  matMtx[pos,rowpos] <- matches
  }
  matMtx[pos,pos] <- NA
}
  
#Create heatmap of overlap between samples

pdf(file = paste(filename,"Prey per nm Analysis.pdf", by= ""), height=27,width=20)
  cols <- colorRampPalette(c("red", "yellow", "white"))(30)
  pheatmap(matMtx, color =  cols, fontsize = 30, na_col = "black",main = paste(filename,"Prey per nm Analysis", by= ""))
dev.off()

#Distance measurement based on 3D modelling in ChimeraX [PDB: 5MRC])
Singh_SS_Dist <- matrix(ncol=ncol(sigint_SS),nrow=ncol(sigint_SS))
XLSXSingh_SS_Dist <- read_excel(distanceinputpath, sheet = 1, na = "N/A")
XLSXSingh_SS_Dist <- XLSXSingh_SS_Dist[,-1]
Singh_SS_Dist <- as.matrix(XLSXSingh_SS_Dist)
rownames(Singh_SS_Dist) <- colnames(Singh_SS_Dist)
View(Singh_SS_Dist)

SSBait_Int2Dist <- (matMtx)/(Singh_SS_Dist)
View(SSBait_Int2Dist)

#Create heatmap of 'habitat' analysis
pdf(file = paste(filename,"Habitat Analysis.pdf", by= ""), height=5,width=5)
  cols <- colorRampPalette(c("red","orange","yellow","white","lightgreen","green","darkgreen"))(100)
  SS_ratio_heatmap <- pheatmap(SSBait_Int2Dist, color =  cols, fontsize = 10, 
                                   main = paste(filename,"Habitat Analysis", by= ""), 
                                   na_col = "grey",
                                   display_numbers = Singh_SS_Dist,
                                   number_col = "black")
  SS_ratio_heatmap
dev.off()
}
sigint_habitat("Singh Nascent Chain","~/Desktop/Singh BioID/Singh_Distances.xlsx")


## Produce coordinate table ####
sigint_links <- function(inputbaitlist){

  setwd("~/Desktop")
  dir.create("BaitCoods")
  setwd("~/Desktop/BaitCoods")
  baitnames <<- colnames(sigint_SS)
  BaitCoodTable <<- function(A, filename){ 
    Cood <- matrix(nrow = length(A), ncol = 2)
    Cood[,1] <- c(seq(1/length(A),by = 1/length(A)))
    Cood[,2] <- A
  
    write.csv(Cood, paste0(filename,".csv"))
  }
  for (i in 1:length(sigint_SS2)){
    BaitCoodTable(sigint_SS2[[i]],paste(baitnames[i],"Coods",sep = ""))
  }
  
  BaitCoodsListed <<- list.files(path="~/Desktop/BaitCoods", full.names = FALSE) %>% 
    read_csv()
    bind_rows()
    
  baitLengths <- c(length(mL57),length(mL67),length(uL23m),length(uL29m)) #Needs changing to read automatically
  BaitCoodsListed[1] <<- rep(baitnames, times = baitLengths)
  colnames(BaitCoodsListed) <<- c("Bait","Cood","Prey")
  
  BaitCoodsListed <<- BaitCoodsListed[duplicated(BaitCoodsListed$Prey) | duplicated(BaitCoodsListed$Prey, fromLast = TRUE), ]
  
  dat1 <<- BaitCoodsListed
  dat2 <<- BaitCoodsListed
  
  names(dat1) <<- c("s1", "c1", "protein")
  names(dat2) <<- c("s2", "c2", "protein")
  
  dat3 <<- merge(dat1, dat2, by = "protein", all = TRUE)
  dat4 <<- subset(dat3, dat3$s1 != dat3$s2)
  
  setwd("~/Desktop")
}
sigint_links(sigint_SS2)
  

## Chord Diagram (Some automation required) ####
sigint_chord <- function(filename){
  pdf(file = filename)
  
  circos.initialize(baitnames, xlim = c(0, 1),
                   sector.width = c("mL57" = length(mL57),"mL67" = length(mL67), "uL23m" = length(uL23m),"uL29m" = length(uL29m))) #needs automating

  {
circos.track(baitnames, ylim = c(0, 1), track.height = 0.0725) #Outer track (1)
circos.axis(h = "bottom", sector.index = "uL29m", labels = "uL29m", major.tick = FALSE, major.at = 0.5,
            labels.facing = "bending.inside", labels.niceFacing = TRUE, labels.cex = 0.7) #Label for uL29m bait
circos.axis(h = "bottom", sector.index = "mL57", labels = "mL57", major.tick = FALSE, major.at = 0.5,
            labels.facing = "bending.inside", labels.niceFacing = TRUE, labels.cex = 0.7) #Label for mL57 bait 
circos.axis(h = "bottom", sector.index = "mL67", labels = "mL67", major.tick = FALSE, major.at = 0.5,
            labels.facing = "bending.inside", labels.niceFacing = TRUE, labels.cex = 0.7) #Label for mL67 bait
circos.axis(h = "bottom", sector.index = "uL23m", labels = "uL23m", major.tick = FALSE, major.at = 0.5,
            labels.facing = "bending.inside", labels.niceFacing = TRUE, labels.cex = 0.7) #Label for uL23m bait
} #Create 1st track (Bait names)
  {
circos.track(baitnames, ylim = c(0, 1), track.height = 0.15) #Inner track (2)
circos.axis(h = "bottom", sector.index = "uL29m", labels = uL29m, major.tick = TRUE, major.at = seq((1/length(uL29m)):1, by = 1/length(uL29m)), 
            labels.facing = "clockwise", labels.niceFacing = TRUE, labels.cex = 0.3) #Labels for uL29m
circos.axis(h = "bottom", sector.index = "mL57", labels = mL57, major.tick = TRUE, major.at = seq((1/length(mL57)):1, by = 1/length(mL57)), 
            labels.facing = "clockwise", labels.niceFacing = TRUE, labels.cex = 0.3) #Labels for mL57
circos.axis(h = "bottom", sector.index = "mL67", labels = mL67, major.tick = TRUE, major.at = seq((1/length(mL67)):1, by = 1/length(mL67)), 
            labels.facing = "clockwise", labels.niceFacing = TRUE, labels.cex = 0.3) #Labels for mL67
circos.axis(h = "bottom", sector.index = "uL23m", labels = uL23m, major.tick = TRUE, major.at = seq((1/length(uL23m)):1, by = 1/length(uL23m)), 
            labels.facing = "clockwise", labels.niceFacing = TRUE, labels.cex = 0.3) #Labels for uL23m
} #Create 2nd track (Prey names)
  totalLength <<- sum(length(mL57),length(mL67),length(uL23m),length(uL29m)) #needs automating
  ChordLinks <- function(input){  
    for (i in 1:nrow(input)){
      startBait <- input[i,2]
      startcood <- input[i,3]
      startcood2 <- input[i,3]
      endBait <- input[i,4]
      endcood <- input[i,5]
      endcood2 <- input[i,5]
      link.col <- rand_color(totalLength)

      
      circos.link(sector.index1 = startBait, c(startcood ,startcood2), 
                  sector.index2 = endBait, c(endcood,endcood2),
                  lty = 1,
                  col = link.col)
      i <- i + 1
    }
  } 
  ChordLinks(dat4)
  circos.clear()
  dev.off()
}
sigint_chord("Test")