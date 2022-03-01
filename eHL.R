eHL <- function(y,P,boot=10, s =.5){
# y is the realization
# P is the prediction
# split: fraction of observations used for estimating q
# ad qq = 10 as inputs when de-commenting the quantile binning e-value HL
  
df <-  data.frame("y"=y, "P"=P)

boot <- foreach(t = 1:boot)%do% {
  
  # split the data set
  
  split <- sample(nrow(df), (nrow(df)*s) , replace = F)
  splitted.1 <- df[split,]
  splitted.2 <- df[-split,]

  ############ eHL based on isotonic regression#################################
  
  # use CORP in the first data to obtain the estimates for q
  
  pav_result <- stats::isoreg(splitted.1$P,splitted.1$y)
  
  red_iKnots <- with(
    pav_result,
    which(!duplicated(yf[iKnots], fromLast = TRUE)) %>% iKnots[.]
  )
  
  df_pav <- with(
    pav_result,
    tibble::tibble(
      case_id = if (isOrd) seq_len(length(y)) else ord,
      x = if (isOrd) x else x[ord],
      y = if (isOrd) y else y[ord],
      bin_id = rep.int(seq_along(red_iKnots), times = diff(c(0, red_iKnots))),
      CEP_pav = yf
    )
  )
  
  df_bins_iso <- tibble::tibble(
    bin_id_iso = seq_along(red_iKnots),
    n = diff(c(0, red_iKnots)),
    x_min = df_pav$x[c(0, utils::head(red_iKnots,-1)) + 1],
    x_max = df_pav$x[red_iKnots],
    CEP_pav = df_pav$CEP_pav[red_iKnots]
  )
  
  suppressMessages(
  eval.isobin <- splitted.2 %>% arrange(P) %>%  tibble() %>% 
    mutate(ints.iso = cut(P, 
                          breaks = c(0,df_bins_iso$x_min[2:(length(df_bins_iso$x_min))],1), # konstant extrapolieren
                          include.lowest = T),
           bin_id_iso = match(ints.iso, unique(ints.iso))   ) %>% 
    left_join(df_bins_iso) %>% 
    mutate(q = if_else(CEP_pav %in% c(0,1), P, CEP_pav )) %>% # here: set q=P if either  x.cal=0 or x.cal=1
    mutate(E = (q^y*(1-q)^(1-y))/(P^y*(1-P)^(1-y)) ) %>% 
    summarise(evalue = prod(E))
  )
  
  eval.isobin
  
}
HLe <- mean(do.call(rbind, boot)$evalue)

return(HLe)
}



