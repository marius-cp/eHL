interpolate <- function(dat, xout) {
  # dat must be the df_bins_iso object
  # artificially include 0 and 1 for most extreme x values
  df <-
    dat %>% mutate(
      x_min = ifelse(x_min==min(x_min),0,x_min),
      x_max = ifelse(x_max==max(x_max),1,x_max),
    )

  df <-
    bind_rows(
      tibble(p=df$x_max,q=df$q),
      tibble(p=df$x_min,q=df$q)
    ) %>% arrange(p)

  stats::approx(x = df$p, y = df$q, xout = xout, method = "linear",ties = "unique")$y
}

eHL <- function(y,P,boot=10, s =.5){
# y is the realization
# P is the prediction
# split: fraction of observations used for estimating q
# ad qq = 10 as inputs when de-commenting the quantile binning e-value HL

df <-  data.frame("y"=y, "P"=P)

boot <- foreach(t = 1:boot)%do% {

  # split the data set

  split <- sample(nrow(df), (nrow(df)*s) , replace = F)
  splitted.1 <- df[split,] %>% arrange(P)
  splitted.2 <- df[-split,] %>% arrange(P)

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

df_bins_iso <-
  df_pav %>% group_by(bin_id) %>%
    summarise(
      a=min(case_id),
      b=max(case_id),
      x_min =min(x),
      x_max = max(x),
      hatpi = sum(y/(b-a+1)) %>% round(10),# round to avoid numerical instabilities!
      infhatpi = (1/(b-a+2)) * (sum(y)+.5)
    ) %>%
    mutate(
      q = ifelse(hatpi%in%c(0,1),infhatpi,hatpi)
    )


eval.isobin <-
  splitted.2 %>%
  mutate(
    q =interpolate(dat=df_bins_iso, xout = P),
    E = (q^y*(1-q)^(1-y))/(P^y*(1-P)^(1-y))
    ) %>%
  summarise(evalue = prod(E))



}


evalues <- do.call(rbind, boot)$evalue
HLe <- mean(evalues)


return(list(HLe=HLe,evalues=evalues))
}



