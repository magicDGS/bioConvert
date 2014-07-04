args <- commandArgs(TRUE)
#args[1] fileroot
#args[2] in_folder
#args[3] out_folder
fileroot<- args[1]
in_folder<-  args[2]
out_folder<-  args[3]

#Read the data
haps <- read.table(paste(in_folder, "/", fileroot, ".haps", sep=""))
sample <- read.table(paste(in_folder, "/", fileroot, ".sample", sep=""))

#Select the columns of genotypes and put the letter
geno <- vector()
for(n in 1:nrow(haps)) {
  f_temp <- as.factor(as.character(haps[n,6:length(haps[1,])]))
  levels(f_temp) <- unlist(haps[n,4:5])
  geno <- rbind(geno, as.character(f_temp))
}
#geno <- as.data.frame(geno)

#Add columns
col1 <- c("I", "A", rep("M", length(geno[,1])))
col2 <- c("ID", "affection", as.character(haps[,2]))
IDs <- as.character(rep(sample[3:length(sample[,2]),2],each=2))
affection <- as.vector(rep(0,length(geno[1,])))
fila1 <- IDs 
fila2 <- affection
geno <- rbind(fila1, fila2, geno, deparse.level=0)
bgl <- cbind(col1, col2, geno)
#write file
write.table(bgl,paste(out_folder, "/", fileroot, ".bgl", sep=""), quote=F,  row.names=F, col.names=F)
