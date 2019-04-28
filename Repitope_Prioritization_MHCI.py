cpu = str(4)
memory = "60G"
pept_len_range = "8:11" ## For MHC-I prediction

home = "C:/Repitope"
target_peptide_file = home + "/peptide_test.txt"
tcr_frag_file = home + "/FragmentLibrary.fst"
mhci_feature_file = "D:/Research/Immunogenicity/FeatureDF_MHCI_Weighted.10000.fst"

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

command = "featureDT_MHCI <- fst::read_fst('" + mhci_feature_file + "', as.data.table=T)"
r(command)

import os
pid = str(os.getpid())
os.mkdir(home + "/ProcessID_" + pid)

command = "res_MHCI <- EpitopePrioritization( \
  featureDF=featureDT_MHCI[Peptide%in%MHCI_Human$Peptide,], \
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],\
  peptideSet=peptideSet_test,\
  fragLib=fragLibDT,\
  aaIndexIDSet='all',\
  fragLenSet=3:8,\
  fragDepthSet=10000,\
  fragLibTypeSet='Weighted',\
  featureSet=MHCI_Human_MinimumFeatureSet,\
  coreN=" + cpu + ",\
  tmpDir='" + home + "/ProcessID_" + pid + "/Temp',\
  seedSet=1:5)"
r(command)

command = "readr::write_csv(res_MHCI, '" + home + "/ProcessID_" + pid + "/EpitopePrioritization_MHCI.csv')"
r(command)	
