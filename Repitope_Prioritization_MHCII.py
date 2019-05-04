cpu = str(4)
memory = "60G"
pept_len_range = "11:30" ## For MHC-II prediction

home = "C:/Repitope"
target_peptide_file = home + "/peptide_test.txt"
tcr_frag_file = home + "/FragmentLibrary.fst"
mhcii_feature_file = home + "/FeatureDF_MHCII_Weighted.10000.fst"

import pyper
r = pyper.R()

command = "options(java.parameters='-Xmx" + memory + "')\n \
library(tidyverse)\n \
library(data.table)\n \
library(Repitope)"
r(command)

command = "peptideSet_test <- Repitope::sequenceFilter(data.table::fread('" + target_peptide_file + "')$Peptide)"
r(command)

command = "peptideSet_test <- peptideSet_test[nchar(peptideSet_test) %in% " + pept_len_range + "]"
r(command)

command = "fragLibDT <- fst::read_fst('" + tcr_frag_file + "', as.data.table=T)"
r(command)

command = "featureDT_MHCII <- fst::read_fst('" + mhcii_feature_file + "', as.data.table=T)"
r(command)

import os
pid = str(os.getpid())
os.mkdir(home + "/ProcessID_" + pid)

command = "res_MHCII <- EpitopePrioritization( \
  featureDF=featureDT_MHCII[Peptide%in%MHCII_Human$Peptide,], \
  metadataDF=MHCII_Human[,.(Peptide,Immunogenicity)],\
  peptideSet=peptideSet_test,\
  fragLib=fragLibDT,\
  aaIndexIDSet='all',\
  fragLenSet=3:11,\
  fragDepth=10000,\
  fragLibType='Weighted',\
  featureSet=MHCII_Human_MinimumFeatureSet,\
  seedSet=1:5,\
  coreN=" + cpu + ",\
  outDir='" + home + "/MHCII_ProcessID_" + pid + "')"
r(command)
	
