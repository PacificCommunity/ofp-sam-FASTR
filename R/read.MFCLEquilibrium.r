#FLR4MFCL - R4MFCL built with FLR classes
#Copyright (C) 2018  Rob Scott

# read.MFCLEquilibrium


read.MFCLEquilibrium <- function(repfile) {
  
  trim.leading  <- function(x) sub("^\\s+", "", x) 
  splitter      <- function(ff, tt, ll=1, inst=1) unlist(strsplit(trim.leading(ff[grep(tt, ff)[inst]+ll]),split="[[:blank:]]+")) 
  
  
  res <- new("MFCLEquilibrium")
  
  pp    <- readLines(repfile)
  pp    <- pp[nchar(pp)>=1]                                          # remove blank lines
  if(any(grepl("# ", pp) & nchar(pp)<3))
    pp <- pp[-seq(1,length(pp))[grepl("# ", pp) & nchar(pp)<3]]   # remove single hashes with no text "# "
  
  eq_mult          <- as.numeric(splitter(pp, "# Effort multiplier", inst=1))
  eq_yield         <- as.numeric(splitter(pp, "# Equilibrium yield", inst=1))
  eq_adult_biomass <- as.numeric(splitter(pp, "# Equilibrium adult biomass", inst=1))
  eq_total_biomass <- as.numeric(splitter(pp, "# Equilibrium total biomass", inst=1))
  res@Eq_calcs     <- data.frame(mult=eq_mult, yield=eq_yield, adult_biomass=eq_adult_biomass, total_biomass=eq_total_biomass)
  
  ypr_mult         <- as.numeric(splitter(pp, "# Effort multiplier", inst=2))
  ypr              <- as.numeric(splitter(pp, "# Yield per recruit", inst=2))
  res@YPR          <- data.frame(mult=ypr_mult, ypr=ypr)
  
  return(res) 
}
