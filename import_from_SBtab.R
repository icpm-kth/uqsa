remotes::install_github("a-kramer/SBtabVFGEN")

import_from_SBtab <- function(SBtabDir){
	setwd(SBtabDir)
	tsvList = dir(pattern = ".*[.]tsv$")
	sbtab_model = SBtabVFGEN::sbtab_from_tsv(tsvList)
	SBtabVFGEN::sbtab_to_vfgen(sbtab_model, cla = FALSE)		#this function created a .vf file (located in the tsvDirectory) from the SBtab model

	vfFileName = dir(pattern = ".*[.]vf$")
	system(paste("cd", SBtabDir))
	system(paste("vfgen r", vfFileName))
	system(paste("vfgen gsl", vfFileName))
	
	return(sbtab_model)
}

