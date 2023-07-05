#MES_sigs_makeup----
library("readxl")
MES_Sigs_List = read_excel("C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/13073_2022_1109_MOESM1_ESM.xlsx", sheet = "Table S1")
MES_functional_list = read_excel("C:/Users/mzaidi/Documents/Visium Datasets/Gene signatures/13073_2022_1109_MOESM1_ESM.xlsx", sheet = "Table S3")

MES_sigs_makeup<- function(MES_Sigs_List,MES_functional_list)
{
  intersect_programs<- lapply(MES_Sigs_List, function(x){lapply(MES_functional_list, function(y){length(intersect(y,x))/length(x) })  })
  
  intersect_programs<- as.data.frame(do.call(cbind,intersect_programs))
  
  return(intersect_programs)
}


df <- MES_sigs_makeup(MES_Sigs_List,MES_Sigs_List)

library(bayesbio)
String_A<-c("John", "is", "going", "to", "the",
               "market", "today", "to", "buy", "cake")
String_B<-c("Tim", "is", "at", "the", "shop",
               "already", "for", "buying", "two", "cakes")
jaccardSets(String_A, String_B)