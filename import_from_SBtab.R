remotes::install_github("a-kramer/SBtabVFGEN")

import_from_SBtab <- function(tsvDirectory){

	setwd(tsvDirectory)
	tsvList = dir(pattern = ".*[.]tsv$")
	sbtab_model = SBtabVFGEN::sbtab_from_tsv(tsvList)
	SBtabVFGEN::sbtab_to_vfgen(sbtab_model)		#this function created a .vf file (located in the tsvDirectory) from the SBtab model

	vfFileName = dir(pattern = ".*[.]vf$")
	system(paste("cd", tsvDirectory))
	system(paste("vfgen r", vfFileName))
}
