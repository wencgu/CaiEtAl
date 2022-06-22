library(tidyverse)
library(openxlsx)

setwd("C:/Google Drive/RStudio files/Feng Cai/2020-12-10/")

CarbonNaturalAbundance <- c(1-0.01082, 0.01082)
HydrogenNaturalAbundance <- c(1-0.000156, 0.000156)
NitrogenNaturalAbundance <- c(1-0.00366, 0.00366)
OxygenNaturalAbundance <- c(1-0.00038-0.00204, 0.00038, 0.00204)

dfTransitions<-
read.xlsx("input.xlsx", sheet = 1, skipEmptyRows = T, skipEmptyCols = T, check.names = F)

listTransitions<-
split(dfTransitions,dfTransitions$setName)

# z<- listTransitions[[7]]
# z


matrixNForRegressionAnalysis<-function(z){
  

vecPositional<-
  expand.grid(rep(list(0:1), 5)) %>% unite(Positional,everything(),sep = "") %>% unlist() 

# glutamate: C5H9NO4 MW 147
# 146/74
# parent: negative C5H8NO4
# daughter: negative C2H4NO2
# abandoned: C3H4O2
stringFormulaParent<- z$formulaParent
stringFormulaAbandoned<-z$formulaAbandoned
stringFormulaDaughter<-z$formulaDaughter

# Required to ensure the "thermo" object is created and defaults are used
suppressMessages(CHNOSZ::reset())
AtomicCompositionParent <- CHNOSZ::makeup(stringFormulaParent)
AtomicCompositionAbandoned <- CHNOSZ::makeup(stringFormulaAbandoned)
  if (is.na(AtomicCompositionAbandoned["C"])){
    AtomicCompositionAbandoned["C"]<-0
  }
  if (is.na(AtomicCompositionAbandoned["H"])){
    AtomicCompositionAbandoned["H"]<-0
  }
  if (is.na(AtomicCompositionAbandoned["N"])){
    AtomicCompositionAbandoned["N"]<-0
  }
  if (is.na(AtomicCompositionAbandoned["O"])){
    AtomicCompositionAbandoned["O"]<-0
  }
AtomicCompositionDaughter <- CHNOSZ::makeup(stringFormulaDaughter)
  if (is.na(AtomicCompositionDaughter["C"])){
    AtomicCompositionDaughter["C"]<-0
  }
  if (is.na(AtomicCompositionDaughter["H"])){
    AtomicCompositionDaughter["H"]<-0
  }
  if (is.na(AtomicCompositionDaughter["N"])){
    AtomicCompositionDaughter["N"]<-0
  }
  if (is.na(AtomicCompositionDaughter["O"])){
    AtomicCompositionDaughter["O"]<-0
  }

vecRowname=NULL
for (i in 0:AtomicCompositionParent["C"] ) { # i the number heavier mass units on the parent ion
  for (j in 0:AtomicCompositionDaughter["C"]){ # j the number heavier mass units on the daughter ion
    if (j>i){
      }
    else{
      tmpRowname= paste0(i,'-',j)
      k<-i-j # k the number heavier mass units on the abandoned ion
      vecRowname<-c(vecRowname,tmpRowname)
      
    }
  }
}

# vecRowname

# length(vecRowname) #12
# length(vecPositional) #32

matrixResults<-matrix(0, nrow = length(vecRowname), ncol = length(vecPositional))
rownames(matrixResults)<-vecRowname
colnames(matrixResults)<-vecPositional
# View(matrixResults)
# matrixResults[vecRowname[2],vecPositional[15]]=29 # works

for (l in vecRowname){
  for (m in vecPositional){
  
  numHeavyAtomAbandoned<-sapply((strsplit(l,"-")),as.numeric)[[1]]-sapply((strsplit(l,"-")),as.numeric)[[2]] 
  numHeavyAtomDaughter<-sapply((strsplit(l,"-")),as.numeric)[[2]] 
  
  # need rework if the daughter ion changes
  
  stringPisitionalDaughter<-substring(m,z$daughterIonStartPosition,z$daughterIonEndPosition)
  num13CDaughter<-sapply(strsplit(stringPisitionalDaughter, ''), function(x) sum(as.numeric(x)))
  num13CTotal<-sapply(strsplit(m, ''), function(x) sum(as.numeric(x)))
  num13CAbandoned<- (num13CTotal- num13CDaughter)
  
  
  numHeavyAtomAbandonedExceptPredetermined13C <- numHeavyAtomAbandoned - num13CAbandoned
  numHeavyAtomDaughterExceptPredetermined13C <- numHeavyAtomDaughter - num13CDaughter
  
  for (n in 0:(AtomicCompositionAbandoned["C"]-num13CAbandoned)){
    if (n > numHeavyAtomAbandonedExceptPredetermined13C){
      break
    }
    else{
      for (o in 0:AtomicCompositionAbandoned["H"]){
        if ((n + o) > numHeavyAtomAbandonedExceptPredetermined13C){
          break
        }
        else{
          for (p in 0:AtomicCompositionAbandoned["N"]){
            if ((n + p + o) > numHeavyAtomAbandonedExceptPredetermined13C){
              break
            }
            else{
              for (q in 0:AtomicCompositionAbandoned["O"]){
                for (r in 0:AtomicCompositionAbandoned["O"]){
                  s <- (q + r*2)
                  if ((q+r)>AtomicCompositionAbandoned["O"]|(n + o + p + s) > numHeavyAtomAbandonedExceptPredetermined13C){
                    break
                  }
                  else{
                    for ( t in 0:(AtomicCompositionDaughter["C"]-num13CDaughter)){
                      if (t > numHeavyAtomDaughterExceptPredetermined13C){
                        break
                      }
                      else{
                        for (u in 0:AtomicCompositionDaughter["H"]){
                          if ((t + u) > numHeavyAtomDaughterExceptPredetermined13C){
                            break
                          }
                          else{
                            for (v in 0:AtomicCompositionDaughter["N"]){
                              if ((t + u + v) > numHeavyAtomDaughterExceptPredetermined13C){
                                break
                              }
                              else{
                                for (w in 0:AtomicCompositionDaughter["O"]){
                                  for (x in 0:AtomicCompositionDaughter["O"]){
                                    y <- (w + x*2)
                                    if ((w + x)>AtomicCompositionDaughter["O"]|(t + u + v + y) > numHeavyAtomDaughterExceptPredetermined13C){
                                      break
                                    }
                                    else{

                                      if ((n+o+p+q+2*r)==numHeavyAtomAbandonedExceptPredetermined13C & (t+u+v+w+2*x) == numHeavyAtomDaughterExceptPredetermined13C){
                                        
                                        pCurrentCAbandoned <- dbinom(n, (AtomicCompositionAbandoned["C"]- num13CAbandoned),CarbonNaturalAbundance[2])
                                        pCurrentCDaughter <- dbinom(t, (AtomicCompositionDaughter["C"]- num13CDaughter),CarbonNaturalAbundance[2])
                                        pCurrentHAbandoned <- dbinom(o, AtomicCompositionAbandoned["H"],HydrogenNaturalAbundance[2])
                                        pCurrentHDaughter <- dbinom(u, AtomicCompositionDaughter["H"],HydrogenNaturalAbundance[2])
                                        pCurrentNAbandoned <- dbinom(p, AtomicCompositionAbandoned["N"],NitrogenNaturalAbundance[2])
                                        pCurrentNDaughter <- dbinom(v, AtomicCompositionDaughter["N"],NitrogenNaturalAbundance[2])
                                        pCurrentOAbandoned <- dmultinom(c((AtomicCompositionAbandoned["O"]-q-r),q,r), AtomicCompositionAbandoned["O"], OxygenNaturalAbundance)
                                        pCurrentODaughter <- dmultinom(c((AtomicCompositionDaughter["O"]-w-x),w,x), AtomicCompositionDaughter["O"], OxygenNaturalAbundance)
                                        
                                        pCurrent <- pCurrentCAbandoned * pCurrentCDaughter * pCurrentHAbandoned * pCurrentHDaughter * pCurrentNAbandoned * pCurrentNDaughter * pCurrentOAbandoned * pCurrentODaughter
                                        
                                        matrixResults[l,m] <- matrixResults[l,m] + pCurrent
                                        
                                      }
                                      
                                      
                                    }
                        
                      }
                    }
                    
                  }
                  }
                }
                
              }
            }
          }
        }
      }
    }
      
  }

  }

        }
      }
    }
  }
  }
}

  
          
          
          # [1] "0-0" "1-0" "1-1" "2-0" "2-1" "2-2" "3-0" "3-1" "3-2" "4-1" "4-2" "5-2"

# View(matrixResults)
return(matrixResults)

}

listResults<-lapply(listTransitions,matrixNForRegressionAnalysis)
View(listResults)

#dplyr::bind_rows(listResults) %>% View()

write.xlsx(listResults,"tmp.xlsx",row.names=T)

# 63.3+17.4+5.4+6.7+1.5+2.5+1.2+0.6+0.7+0.1
# 63.9+16.9+5.5+6.6+1.4+2.5
# 64.1+1.5+20.9+0.1+0.5+9.8+0.2+2.1

