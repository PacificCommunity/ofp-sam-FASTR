
setGeneric('eq_calcs', function(frq, par, rep, ...) standardGeneric('eq_calcs')) 


setMethod("eq_calcs", signature(frq='MFCLFrq', par="MFCLPar", rep="MFCLRep"),
          function(frq, par, rep, inc=0.015, fsh=NULL, ...){
            age <- 1:dimensions(par)['agecls']
            wgt <- waa(par)/1000
            m   <- c(aperm(m_at_age(rep), c(4,1,2,3,5,6)))
            mat <- c(aperm(mat(par), c(4,1,2,3,5,6)))

            alpha <- c(srr(rep)['a'])
            beta  <- c(srr(rep)['b'])
            sigma <- c(srr(rep)['sigma'])

            fyrange <- range(parbc)['maxyear'] - flagval(parbc, 2, 148)$value/flagval(parbc,2,57)$value
            lyrange <- range(parbc)['maxyear'] - flagval(parbc, 2, 155)$value/flagval(parbc,2,57)$value

            cdat    <- freq(frq)[freq(frq)$year %in% fyrange:lyrange & (freq(frq)$length %in% lf_range(frq)['LFFirst'] 
                                                                     |  freq(frq)$weight %in% lf_range(frq)['WFFirst']
                                                                     |  freq(frq)$freq == -1), ]
            cpropn <- tapply(cdat$catch, cdat$fishery, sum)/sum(cdat$catch)
            if(!is.null(fsh))
              cpropn <- cpropn[is.element(as.numeric(names(cpropn)), fsh)]

            selection <- apply(sweep(sel(rep)[,,as.numeric(names(cpropn)),], 3, cpropn, "*"), c(1), sum)
            sel       <- pmax(selection/max(selection), 1E-10)

            fvals     <- areaSums(popN(rep)*fm(rep))/areaSums(popN(rep))
            fvals     <- seasonMeans(yearMeans(fvals[,as.character(fyrange:lyrange)]))
            fsel      <- fvals/max(fvals)

            res <- YPR(age, sel, m, mat, wgt, alpha, beta, sigma, inc)
            res <- res/flagval(par, 2, 57)$value

            return(res)
          })

