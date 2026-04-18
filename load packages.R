#================================================================
# Sat Mar 21 16:27:03 2026
#================================================================

#================================================================
# load R packages
#================================================================

if (!require("Sim.DiffProc")) { install.packages("Sim.DiffProc"); library("Sim.DiffProc") }
if (!require("yuima"))        { install.packages("yuima");        library("yuima")        }
if (!require("ttutils"))      { install.packages("ttutils");      library("ttutils")      }
if (!require("EQL"))          { install.packages("EQL");          library("EQL")          }
if (!require("gsl"))          { install.packages("gsl");          library("gsl")          }
if (!require("ggsci"))        { install.packages("ggsci");        library("ggsci")        }
if (!require("furrr"))        { install.packages("furrr");        library("furrr")        }
if (!require("Deriv"))        { install.packages("Deriv");        library("Deriv")        }
if (!require("zoo"))          { install.packages("zoo");          library("zoo")          }
if (!require("ggplot2"))      { install.packages("ggplot2");      library("ggplot2")      }
if (!require("tidyr"))        { install.packages("tidyr");        library("tidyr")        }
if (!require("dplyr"))        { install.packages("dplyr");        library("dplyr")        }
if (!require("scales"))       { install.packages("scales");       library("scales")       }
if (!require("MASS"))         { install.packages("MASS");         library("MASS")         }
