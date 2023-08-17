#
library(GAPIT3)
#
root <- "E:/gitfolder/AGILE_LDP_Protein"
setwd(root)
myG <- read.csv("Data/myG_LDP.csv", header = F)
myCV <- read.csv("Data/myCV.csv")
#
myY <- read.csv("Data/myY_lsmeans.csv")
setwd(paste0(root,"/GWAS_Results_lsmeans"))
myGAPIT <- GAPIT(
  Y = myY[,c(1,)],
  G = myG,
  PCA.total = 4,
  model = c("MLM","FarmCPU","Blink"),#
  Phenotype.View = F
)
setdiff(t(myG[1,]), myY$Name)
setdiff(myY$Name, t(myG[1,]))
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 4,
  model = c("MLMM"),#
  Phenotype.View = F
)
#
setwd(root)
myY <- read.csv("Data/myY_Protein_NIRS.csv")
colnames(myY)
setwd(paste0(root,"/GWAS_Results"))
myGAPIT <- GAPIT(
  Y = myY[,c(1,58:ncol(myY))],
  G = myG,
  PCA.total = 4,
  model = c("MLMM"),#
  Phenotype.View = F
)

#
#
setwd(root)
myY <- read.csv("Data/myY_NIRS_Ro16.csv")
setwd(paste0(root,"/GWAS_Results_DTF_REP"))
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  CV = myCV[,c("Name","DTF_Ro16","REP_Ro16")],
  PCA.total = 0,
  model = c("MLMM"),#
  Phenotype.View = F
)
setwd(root)
myY <- read.csv("Data/myY_NIRS_Ro17.csv")
myEntries <- !is.na(myCV$REP_Ro17)
setwd(paste0(root,"/GWAS_Results_DTF_REP"))
myGAPIT <- GAPIT(
  Y = myY[myEntries,],
  G = myG,
  CV = myCV[myEntries,c("Name","DTF_Ro17","REP_Ro17")],
  PCA.total = 0,
  model = c("MLMM"),#
  Phenotype.View = F
)
setwd(root)
myY <- read.csv("Data/myY_NIRS_Su16.csv")
setwd(paste0(root,"/GWAS_Results_DTF_REP"))
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  CV = myCV[,c("Name","DTF_Su16","REP_Su16")],
  PCA.total = 0,
  model = c("MLMM"),#
  Phenotype.View = F
)
setwd(root)
myY <- read.csv("Data/myY_NIRS_Su17.csv")
myEntries <- !is.na(myCV$REP_Su17)
setwd(paste0(root,"/GWAS_Results_DTF_REP"))
myGAPIT <- GAPIT(
  Y = myY[myEntries,],
  G = myG,
  CV = myCV[myEntries,c("Name","DTF_Su17","REP_Su17")],
  PCA.total = 0,
  model = c("MLMM"),#
  Phenotype.View = F
)
