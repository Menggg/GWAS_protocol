#please replace all the pathway in this script to the pathway saved your data
#create a new map file for PLINK binary format to replace allele from A/B to G/C/A/T
setwd("~/Downloads")
map=read.table("HDRA-G6-4-SNP-MAP/HDRA-G6-4-final-snp-map.tsv",head=F)
bim=read.table("HDRA-G6-4-RDP1-RDP2-NIAS.bim",head=F)
bim[,5:6]=map[,4:5]
write.table(bim,"HDRA-G6-4-RDP1-RDP2-NIAS.bim",quote=F,row.names = F,col.names = F)

#convert PLINK format to VCF format
system("./plink2 --bfile HDRA-G6-4-RDP1-RDP2-NIAS --export vcf --out rice")

#run beagle5 to impute missing genotype. beagle5 could be downloaded from https://faculty.washington.edu/browning/beagle/beagle.html
system("java -Xmx20000m -jar beagle.jar gt=rice.vcf out=test")

#convert imputed data to sample major PLINK text format
system("./plink2 --vcf test.vcf.gz --export A --out rice")

#load SNP data, map, and phenotype
snp=read.table("rice.raw",head=T)
myGM=read.table("HDRA-G6-4-RDP1-RDP2-NIAS.bim",head=F)[,c(2,1,4)]
trait1=read.table("phenoAvLen_G6_4_RDP12_ALL.txt",head=T)

#match and convert sample ID between SNP data and phenotype data
ID=snp[,1:2]
write.table(ID,"temp.txt",quote=F,row.names = F,col.names = F)
ID=read.table("temp.txt",head=F,sep='_')
trait1[,1]=paste("IRGC",trait1[,1],sep='')
trait=read.delim("HDRA-G6-4-RDP1-RDP2-NIAS-sativa-only.sample_map.rev2.tsv",head=F)
trait=merge(trait,ID,by.x="V2",by.y="V2")
trait=merge(trait,trait1,by.x="V3",by.y="FID")
myGD=snp[,-c(1,3:6)]
trait[,2]=paste(trait[,2],trait[,3],sep='_')
myY=trait[,c(2,13)]

#save input data to Rdata format that is easy to re-load
save(myGD,myGM,myY,file="rice.RData")
load("rice.RData")

#run BLINK, FarmCPU, and MLM
source("http://zzlab.net/GAPIT/gapit_functions.txt")
myGAPIT_blink <- GAPIT( Y=myY, GD=myGD, GM=myGM, model="Blink", SNP.MAF=0.01)
myGAPIT_farmcpu <- GAPIT( Y=myY, GD=myGD, GM=myGM, model="FarmCPU", SNP.MAF=0.01)
myGAPIT_mlm <- GAPIT( Y=myY, GD=myGD, GM=myGM, model="MLM", SNP.MAF=0.01)






