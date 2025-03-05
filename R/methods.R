
setGeneric('eq_calcs', function(frq, par, rep, ...) standardGeneric('eq_calcs')) 


setMethod("eq_calcs", signature(frq='MFCLFrq', par="MFCLPar", rep="MFCLRep"),
          function(frq, par, rep, inc=0.015, fsh=NULL, fyr=NULL, lyr=NULL, selection=TRUE, ...){

            age <- 1:dimensions(par)['agecls']
            wgt <- c(aperm(mean_waa(rep), c(4,1,2,3,5,6)))  #waa(par)
            m   <- m_at_age(rep)
            mat <- mat(par)

            alpha <- c(srr(rep)['a'])
            beta  <- c(srr(rep)['b'])
            sigma <- c(srr(rep)['sigma'])

            fyrange <- range(par)['maxyear'] - flagval(par, 2, 148)$value/flagval(par,2,57)$value
            if(!is.null(fyr))
              fyrange <- fyr
            lyrange <- range(par)['maxyear'] - flagval(par, 2, 155)$value/flagval(par,2,57)$value
            if(!is.null(lyr))
              lyrange <- lyr

            if(selection){
              cdat    <- freq(frq)[freq(frq)$year %in% fyrange:lyrange & (freq(frq)$length %in% lf_range(frq)['LFFirst'] 
                                                                       |  freq(frq)$weight %in% lf_range(frq)['WFFirst']
                                                                       |  freq(frq)$freq == -1), ]
              cpropn <- tapply(cdat$catch, cdat$fishery, sum)/sum(cdat$catch)
              if(!is.null(fsh))
                cpropn <- cpropn[is.element(as.numeric(names(cpropn)), fsh)]

              selectionp<- apply(sweep(sel(rep)[,,as.numeric(names(cpropn)),], 3, cpropn, "*"), c(1), sum)
              sel       <- pmax(selectionp/max(selectionp), 1E-10)
            }
            
            if(!selection){
              fvals     <- areaSums(popN(rep)*fm(rep))/areaSums(popN(rep))
              fvals     <- seasonMeans(yearMeans(fvals[,as.character(fyrange:lyrange)]))
              sel       <- fvals#/max(fvals)
            }
            
            res <- YPR(age, sel, m, mat, wgt, alpha, beta, sigma, inc)
            res <- res/flagval(par, 2, 57)$value

            return(res)
          })

