paquetes<-c("ape", "picante", "pez", "phytools","vegan", "adephylo", "phylobase"
            , "geiger", "mvMORPH", "OUwie", "hisse", "BAMMtools","phylosignal"
            , "Biostrings","devtools","ggplot2", "kableExtra", "betapart"
            , "gridExtra","reshape2")

paquetes.instalar <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(paquetes.instalar)) install.packages(paquetes.instalar)