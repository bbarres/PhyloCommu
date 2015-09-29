###############################################################################
###############################################################################
#Preparing Kinship matrices for the Mixed linear model analysis
###############################################################################
###############################################################################

#loading the libraries
library(kinship2)

#set the working directory
setwd("~/work/Rfichiers/Githuber/PhyloCommu_data")


###############################################################################
#Kinship matrix based on the pedigree
###############################################################################

pedifile<-read.table("commu_pedi.txt",header=T)
commu_kin_ped<-2*kinship(id=pedifile$id,mom=pedifile$mother,dad=pedifile$father)
write.table(commu_kin_ped,file="commu_kin_ped.txt",sep='\t',quote=FALSE)
#clean the environment
rm(pedifile,commu_kin_ped)


###############################################################################
#Kinship matrix based on the pedigree
###############################################################################

#Here are some codes to transform the output of the Spagedi 1.5a software into 
#a matrix suitable for the analysis in TASSEL. For that you just need to run 
#Spagedi choosing the relatedness index you are interested in and then extract 
#from the output file the matrix of pairwise kinship coefficient and the list 
#of inbreeding coefficient of each individual. We then combine this two files, 
#putting in the diagonal of the matrix 0.5*(1+inbreeding coefficient)

#here is the function to insert the inbreeding coefficient in the kinship 
#matrix
insertinbreed<-function(kin,inbreed)
{
  for (i in 1:dim(kin)[1]) {
    kin[i,i]<-inbreed$intrakin[inbreed$Name.i==rownames(kin)[i]]
  }
  ifelse(min(kin)<0,kin<-kin+abs(min(kin)),kin)
  return(kin)
}


#We process the file for the 'lim' treatment, using a kinship computed with 
#610 SNPs with a LD lower than 0.8 and a MAF>0.05
#first we load the kinship matrix
kinLOIS_commu<-read.table("commu_ldmaf_LOIS.txt",header=TRUE,
                          row.names=1,sep="\t")
#then we load the inbreeding coefficient and turn it into 0.5*(1+inbreed coeff)
kinintra_commu<-read.table("commu_ldmaf_intraLOIS.txt",header=T,sep="\t")
kinintra_commu<-cbind(kinintra_commu,
                      "intrakin"=0.5*(1+kinintra_commu$ALL.LOCI))
#finally we used the function to include the inbreeding coefficient value in 
#the matrix and export the dataset to be used in TASSEL
transLOIScommu<-2*insertinbreed(kinLOIS_commu,kinintra_commu)
write.table(transLOIScommu,file="commu_kin_LOIS.txt",sep='\t',quote=FALSE)
#clean the environment
rm(kinLOIS_commu,kinintra_commu,transLOIScommu)


###############################################################################
#END
###############################################################################