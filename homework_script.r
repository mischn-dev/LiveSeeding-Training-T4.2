

## read files 
      setwd("path to the folder of the downloaded data")

      geno = read.csv("./Liveseeding_training_genotypes.csv",  na.strings = c("NA", "./."))
      pheno = read.csv("./Liveseeding_training_phenotypes.csv")

    ## sample a subset which to use - 500 markers and 80 genotypes
        subset = unique(pheno$Hybrid)[....]
        geno = geno[geno[,1] %in% subset,]
        pheno = pheno[pheno$Hybrid %in% subset,]

        all(geno[,1] == pheno$Hybrid) # true
        geno = geno[,-1]
   
    ## prep files 
        pheno$Location = as.factor(pheno$Location)
        pheno$Location = as.integer(as.factor(pheno$Location))

    ## load library 
        if("rrBLUP" %in% installed.packages() == F){install.packages("rrBLUP")}
        library(rrBLUP) 

## run the rrblub ############################################################################################################
            geno[is.na(geno)] = 0 # replacing all missing markers with the heterozygote state
            se = mixed.solve(pheno$PlantHeight, ..... , X = pheno$Location)

            ## 
            # test the correlations between the predicted and assessed plant height
            u = as.vector(se$u)
            
            # use this to generate a genomic blup value for each genotype
            blup = as.matrix(geno) %*% u
            
            # construct an estimate for the environment and the time as well 
            env = pheno$Location * se$beta[1]
        
            # add the fixed effects environment and time to the blup value, to retain an actual estimate of the chlorophyll content 
            se_estimate = env + as.vector(blup)
            
            # compare to the real values we used to regress the model 
            cor(pheno$PlantHeight, se_estimate)
            plot(pheno$PlantHeight, se_estimate)


## run the gblub model #####################################################################################################
            g2 = rrBLUP::A.mat(geno, impute.method = "mean", return.imputed = F)

            gb = mixed.solve(pheno$PlantHeight, ..... , X = pheno$Location)

            # get the estimated plant height
            blup = gb$u
            gb_estimate = as.vector(gb$beta) + blup
            
            # compare to the real values we used to regress the model 
            cor(pheno$PlantHeight, gb_estimate)
            plot(pheno$PlantHeight, gb_estimate)
