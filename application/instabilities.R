
HL.binning <- function(FC, rlz, binning.method="Q1", m.bins=10) {
  df <- data.frame(FC,rlz)
  suppressMessages(

    if (binning.method == "Q1"){
      bin.breaks <- unique(c(0,quantile(FC, (1:(m.bins-1))/m.bins),1))
      df <- df %>% mutate(bin_id = findInterval(FC, bin.breaks))
    } else if (binning.method == "RS"){
      bin.breaks <- unique(quantile(FC, probs=seq(0, 1, 1/m.bins)))
      df <- df %>% mutate(bin_id = findInterval(FC, bin.breaks, rightmost.closed=TRUE,left.open = TRUE))
    }else if (binning.method == "Q2"){
      df <- df %>%
        dplyr::arrange(FC) %>%
        mutate(bin_id = findInterval(row_number(), seq(1,dim(df)[1],length.out=m.bins+1), rightmost.closed=TRUE,left.open = TRUE))
    } else if (binning.method == "Q3"){
      df <- df %>%
        dplyr::arrange(FC, rlz) %>%
        mutate(bin_id = findInterval(row_number(), seq(1,dim(df)[1],length.out=m.bins+1), rightmost.closed=TRUE,left.open = TRUE))
    } else if (binning.method == "Q4"){
      df <- df %>%
        dplyr::arrange(FC, desc(rlz)) %>%
        mutate(bin_id = findInterval(row_number(), seq(1,dim(df)[1],length.out=m.bins+1), rightmost.closed=TRUE,left.open = TRUE))
    } else if (binning.method == "equal"){
      bins <-unlist(classInt::classIntervals(FC, m.bins, style = 'equal')[2])
      df <- df %>% mutate(bin_id = findInterval(FC, bins, rightmost.closed=TRUE,left.open = TRUE))
    }
  )

  # Compute bin midpoint and average:
  df <- df %>% group_by(bin_id) %>% mutate(bin.midpoint=median(FC), bin.average=mean(FC)) %>% ungroup()

  # Generate df.bins
  suppressMessages(
    df.bins <-
      df %>% group_by(bin_id) %>%
      dplyr::summarise(bin.freq=mean(rlz),
                       bin.n=length(rlz),
                       bin.n.small=(length(rlz)<=5),
                       bin.midpoint=unique(bin.midpoint),
                       bin.average=unique(bin.average))
  )
  # Generate a df.breaks with bin breaks for the histogram!
  suppressMessages(
    df.breaks <- df %>% group_by(bin_id) %>%
      summarize(min.FC = min(FC), max.FC =max(FC))
  )

  suppressMessages(
    bin.breaks <- unique(c(0,(df.breaks$min.FC[-1] + df.breaks$max.FC[-length(df.breaks$max.FC)])/2,1))
  )

  suppressMessages(
    Ctab <- df %>% group_by(bin_id) %>% summarize(min.FC = min(FC),
                                                  max.FC =max(FC),
                                                  o1 = sum(as.numeric(rlz)),
                                                  e1 = sum(FC),
                                                  e0 = sum((1-FC)),
                                                  o0 = sum(1-as.numeric(rlz)) )%>%
      mutate(total = o1+o0)

  )

  if(binning.method == "equal"){
    Chat <-  Ctab %>% dplyr::summarise(Chat = ((o1-e1)^2/e1 + (o0-e0)^2/e0 )) %>% .[complete.cases(.),] %>% sum()
  }else{
    suppressMessages(
      Chat <- Ctab %>% dplyr::summarise(Chat = sum((o1-e1)^2/e1 + (o0-e0)^2/e0 ))
    )
  }


  PVAL = 1 - pchisq(as.numeric(Chat), m.bins) # OOS HL pval


  return(list(df=df, df.breaks=df.breaks,
              bin.breaks = bin.breaks,
              Ctab = Ctab ,Chat = Chat, PVAL=PVAL) )
}

