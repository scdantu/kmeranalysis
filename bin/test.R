d<-read.csv('TMut_Out.txt',header=F)
colnames(d)<-c("variant","category","class","kmer_count")
head(d)
table(d$category)
table(d$class)
sum(d$kmer_count)
savehistory(file="test.R")
