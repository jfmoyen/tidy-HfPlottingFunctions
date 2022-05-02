# Utility functions to plot Hf evolution lines and such


.onLoad<-function(lib, pkg){
  on.exit(options("show.error.messages"=TRUE))
  options("warn"=-1)

  packageStartupMessage("Routines for Hf isotope data plotting\r
                        All ages are in Ma", appendLF = TRUE)

  if(require(tidyverse)){
    packageStartupMessage("Tidyverse available")
    tidyverse_available <- T
    library(magrittr) # Manual loading for fancy pipes
  }else{
    packageStartupMessage("Tidyverse not available - expect trouble")
    tidyverse_available <- F
  }

  HfCst <- read.table("data/HfConstants.txt",sep="\t",header=T)

  if(tidyverse_available){
    HfCst %<>% as_tibble()
  }

  assign("HfCst", HfCst, .GlobalEnv)

  invisible()
}


####### Constants ########

##### Functions #####
DMevolution<-function(from,to,DMmodel,step){

}

