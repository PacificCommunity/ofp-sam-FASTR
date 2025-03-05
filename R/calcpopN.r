

###################################################
## Methods to Calculate popN by region for historical, projection and tag periods

setGeneric('histNa', function(par, rep, Fa, ...) standardGeneric('histNa'))

#for the historical estimates
setMethod('histNa', signature(par='MFCLPar', rep='MFCLRep', Fa='FLQuant'),
          function(par, rep, Fa, ...){

            maxage <- dimensions(par)['agecls']
            maxssn <- dimensions(par)['seasons']
            maxyr  <- dimensions(par)['years']/dimensions(par)['seasons']
            minyr  <- 2

            ssnmap    <- c(4,1,2,3)
            yrmap     <- c(-1,0,0,0)
            
            Na <- FLQuant(0, dimnames=dimnames(popN(rep)))
            Na[,1,] <- popN(rep)[,1,]  # seed the numbers at age with the first year from the rep file.
            
            Ma <- expand(FLQuant(m_at_age(rep), dimnames=list(age=1:40)), year=dimnames(Na)$year, season=dimnames(Na)$season, area=dimnames(Na)$area)

            for(year in minyr:maxyr){
            # recruitment
            Na[1,year,,,] <-  c(exp(tot_pop(par))) * sweep(sweep(exp(region_rec_var(par)[,year,]), 5, region_pars(par)[1,], "*"), 4, rel_rec(par)[,year,], '*')

            for(ssn in 1:maxssn){
              for(age in 2:maxage){
                Na[age,year,,ssn,,] <- (diff_coffs_age_period(par)[,,age-1,ssnmap[ssn]] %*% Na[age-1,year+yrmap[ssn],,ssnmap[ssn],,] *
                                          c(exp(-Ma[age-1,year+yrmap[ssn],,ssnmap[ssn],,]-Fa[age-1,year+yrmap[ssn],,ssnmap[ssn],,])))
              }
              Na[age,year,,ssn,,] <- Na[age,year,,ssn,,] + c(diff_coffs_age_period(par)[,,age,ssnmap[ssn]] %*% Na[age,year+yrmap[ssn],,ssnmap[ssn],,]) *
                c(exp(-Ma[age,year+yrmap[ssn],,ssnmap[ssn],,]-Fa[age,year+yrmap[ssn],,ssnmap[ssn],]))
            }
          }
          return(Na)
        })




setGeneric('projNa', function(par, rep, Fa, Na, simyears, ...) standardGeneric('projNa'))

#for the historical estimates
setMethod('projNa', signature(par='MFCLPar', rep='MFCLRep', Fa='FLQuant', Na='FLQuant', simyears='numeric'),
            function(par, rep, Fa, Na, simyears, Wa=NULL, Ma=NULL, Mat=NULL, lnormrec=NULL, srrdev=NULL, rectype='projection', ...){
#browser()
          maxage <- dimensions(par)['agecls']
          maxssn <- dimensions(par)['seasons']
          maxyr  <- dimensions(par)['years']/dimensions(par)['seasons']
          minyr  <- maxyr - length(simyears) + 1
          
          pf232 <- floor(flagval(par, 1, 232)$value/4)
          pf233 <- floor(flagval(par, 1, 233)$value/4)

          ssnmap    <- c(4,1,2,3)
          yrmap     <- c(-1,0,0,0)

          if(is.null(Ma))
            Ma <- expand(FLQuant(m_at_age(rep), dimnames=list(age=1:maxage)), year=dimnames(Na)$year, season=dimnames(Na)$season, area=dimnames(Na)$area)
          if(is.null(Wa))
            Wa <- expand(FLQuant(waa(par), dimnames=list(age=1:maxage)), year=dimnames(Na)$year, season=dimnames(Na)$season, area=dimnames(Na)$area)
          if(is.null(Mat))
            Mat <- expand(FLQuant(mat(par), dimnames=list(age=1:maxage)), year=dimnames(Na)$year, season=dimnames(Na)$season, area=dimnames(Na)$area)

          if(is.null(lnormrec))
            lnormmrec <- log(apply(rec_region(rep),c(4,5),mean)/sum(apply(rec_region(rep),c(4,5),mean)/4))

          if(is.null(srrdev)){
            srrSS <- seasonMeans(areaSums(quantSums(Na[,simyears,,,] * Wa[,simyears,,,] * Mat[,simyears,,,])))/1000
            srrRR <- ((srr(prep)['a']*srrSS)/(srr(prep)['b']+srrSS)) * exp(srr(prep)['sigma']/2)
            #srrdev <- c(mean(eq_rec(rep)[1,pf232:pf233,])) - eq_rec(rep)[1,simyears]
            srrdev <- srrRR - eq_rec(rep)[,simyears]
          }

          for(year in minyr:maxyr){
            SSfull        <- seasonMeans(areaSums(quantSums(Na[,year-1,,,] * Wa[,year-1,,,] * Mat[,year-1,,,])))/1000
            Na[1,year,,,] <- c(((srr(rep)['a']*SSfull)/(srr(rep)['b']+SSfull)) * exp(srr(rep)['sigma']/2) - srrdev[,year-minyr+1])/4 * exp(lnormmrec)

            for(ssn in 1:maxssn){
              for(age in 2:maxage){
                Na[age,year,,ssn,,] <- (diff_coffs_age_period(par)[,,age-1,ssnmap[ssn]] %*% Na[age-1,year+yrmap[ssn],,ssnmap[ssn],,] *
                                          c(exp(-Ma[age-1,year+yrmap[ssn],,ssnmap[ssn],,]-Fa[age-1,year+yrmap[ssn],,ssnmap[ssn],,]))) 
              }
              Na[age,year,,ssn,,] <- Na[age,year,,ssn,,] + c(diff_coffs_age_period(par)[,,age,ssnmap[ssn]] %*% Na[age,year+yrmap[ssn],,ssnmap[ssn],,]) *
                c(exp(-Ma[age,year+yrmap[ssn],,ssnmap[ssn],,]-Fa[age,year+yrmap[ssn],,ssnmap[ssn],]))
            }
          }
          return(Na)
        })


