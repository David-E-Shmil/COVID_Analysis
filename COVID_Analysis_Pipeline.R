#Applying all the packages that will be used
lapply(c("ape","seqinr","Biostrings", "tidyverse","msa", "phangorn",
         "phytools","stringr","taxize"), library, character.only = TRUE)
library("ggmsa")
library(dendextend)

#Reading in the two fasta files and putting them to a variable
covid_gen <- readDNAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/SARS_MERS_coronavirus.raw_sequence.fasta")
covid_spk <- readAAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/spike.fa")

#Making sure both fasta files share the same variants of Covid
covid_gen <- covid_gen[which(names(covid_gen)%in%names(covid_spk))]
covid_spk <- covid_spk[which(names(covid_spk)%in%names(covid_gen))]

#Converting the variables into msa objects
covid_gen.msa <- msa(covid_gen, "ClustalW")
covid_spk.msa <- msa(covid_spk, "ClustalW")

#Converting the msa objects into phyDat objects
covid_gen.phy <- as.phyDat(msaConvert(covid_gen.msa,"seqinr::alignment"),type = "DNA")
covid_spk.phy <- as.phyDat(msaConvert(covid_spk.msa,"seqinr::alignment"),type = "AA")

#Creating a distance model for both phyDat objects using JC69 and JTT respectively
covid_gen.phy_dml <-dist.ml(covid_gen.phy, model = "JC69")
covid_spk.phy_dml <-dist.ml(covid_spk.phy, model = "JTT")

#Using phangorn to create neighbour joining trees of both variables
covid_gen_nj <- phangorn::NJ(covid_gen.phy_dml)
covid_spk_nj <- phangorn::NJ(covid_spk.phy_dml)

#Rooting both NJ using a midpoint
covid_gen_nj <- midpoint(covid_gen_nj)
covid_spk_nj <- midpoint(covid_spk_nj)

#Bootstrapping then plotting both NJ trees and setting the seed to ensure that the results are repeatable
set.seed(100)
covid_gen_nj_bs <-bootstrap.phyDat(covid_gen.phy, FUN = function(x)NJ(dist.ml(x, model = "JC69")),bs=1000)
gen_plot<-plotBS(covid_gen_nj,covid_gen_nj_bs,"phylogram")
title("Covid Genome Neighbour Joining Tree")


set.seed(100)
covid_spk_nj_bs <-bootstrap.phyDat(covid_spk.phy, FUN = function(x)NJ(dist.ml(x, model = "JTT")),bs=1000)
spk_plot<-plotBS(covid_spk_nj,covid_spk_nj_bs,"phylogram")
title("Covid Spike Gene Neighbour Joining Tree")

#Finding out which model works best for both phyDat objects based on lowest AIC
mt1 <- modelTest(covid_gen.phy, model="all")
mt1$Model[mt1$AIC==min(mt1$AIC)]

mt2 <- modelTest(covid_spk.phy, model="all")
mt2$Model[mt2$AIC==min(mt2$AIC)]


#Creating a Maximum liklihood of both NJ objects using the optimal model from before
covid_gen_NJ.pml <- pml(covid_gen_nj, covid_gen.phy ,model = "GTR", k=4, inv = .2)
covid_spk_NJ.pml <- pml(covid_spk_nj, covid_spk.phy ,model = "WAG", k=4, inv = .2)

#Bootstrapping bot ML objects and setting a seed to ensure that the results are repeatable
covid_gen_NJ.pml <- optim.pml(covid_gen_NJ.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)
set.seed(100)
covid_gen_NJ.pml.bs <- bootstrap.pml(covid_gen_NJ.pml,bs=100,trees=TRUE,optNni=TRUE)

covid_spk_NJ.pml <- optim.pml(covid_spk_NJ.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)
set.seed(100)
covid_spk_NJ.pml.bs <- bootstrap.pml(covid_spk_NJ.pml,bs=100,trees=TRUE,optNni=TRUE)

#Rooting the trees using the midpoint method
covid_gen_NJ.pml.bs <- midpoint(covid_gen_NJ.pml.bs)
covid_spk_NJ.pml.bs <- midpoint(covid_spk_NJ.pml.bs)


#Plotting both bootstrapped ML trees
plotBS(covid_gen_NJ.pml$tree, covid_gen_NJ.pml.bs, type = "phylogram")
title("Covid Genome ML Tree")

plotBS(covid_spk_NJ.pml$tree, covid_spk_NJ.pml.bs, type = "phylogram")
title("Covid Spike Gene ML Tree")


#Creating a tanglgram for 3 different comparisons of NJ v ML for both the genome and gene trees and ML v ML for the genome and gene tree
obj1 <- cophylo(covid_gen_nj, covid_gen_NJ.pml$tree)

par(mfrow=c(1,1))

plot(obj1,mar=c(.1,.1,5,.1))
title("NJ                                                                   ML")

obj2 <- cophylo(covid_spk_nj, covid_spk_NJ.pml$tree)

par(mfrow=c(1,1))

plot(obj2,mar=c(.1,.1,5,.1))
title("NJ                                                                   ML")

obj3 <- cophylo(covid_gen_NJ.pml$tree, covid_spk_NJ.pml$tree)

par(mfrow=c(1,1))

plot(obj3,mar=c(.1,.1,5,.1))
title("Gen ML                                                                   Spk ML")

