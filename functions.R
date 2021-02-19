## Forest plot from odds/hazard ratios. 
# xlabs = labels for the x axis (coordinates are flipped here so x is the vertical axis)
# title = plot title 
# breaks = breaks in the y axis (correspond to spacing of the elements on the plot)

make_forest <- function (dataframe, xlabs, title, breaks) {
  
  
  ymax <- max(dataframe$ci_high)
  print(ymax)
  ylim <- ceiling(max(dataframe$ci_high)) 
  print(ylim)
  xlength <- length(dataframe$predictor)
  print(xlength)
  
  forest <- ggplot (data = dataframe, aes (x = predictor, y = oddsratio)) +
    geom_point(shape = 18, size = 3.5) +
    geom_errorbar(aes(ymax = ci_high, ymin = ci_low), width = 0.5 , size = 0.5) +
    geom_hline(yintercept = 1, linetype = "longdash") +
    coord_flip(clip = "off", y = c(0, ymax + ymax/1.6 + 0.1)) +
    scale_y_continuous(breaks = c(seq(0, ylim, by = breaks)), labels = c(seq(0, ylim, by = breaks))) +
    scale_x_discrete(labels = xlabs) +
    xlab(NULL) + ylab(NULL) + 
    ggtitle (title) +
    theme_bw() +
    theme (panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
           axis.text=element_text(colour = "black"),
           plot.title = element_text(size = 12)) +
    geom_text(aes(y = ymax + ymax/10 , label = or), size = 3.5) + 
    geom_text (aes (y = ymax + ymax/3.3, label = cirange), size = 3.5) + 
    geom_text (aes(y = ymax + ymax/1.6, label = p_value), size = 3.5) +
    geom_text (x = xlength + 1.1, y = ymax +ymax/10 , label = "OR", alpha = 0.1) +
    geom_text (x = xlength + 1.1, y = ymax + ymax/3.3, label = "CI", alpha = 0.1) +
    geom_text (x = xlength + 1.1, y = ymax + ymax/1.6, label = "p", alpha = 0.1) 
  
  print(forest)
}

## Forest plot from odds/hazard ratios - y axis is on a logarithmic scale. 
# xlabs = labels for the x axis (coordinates are flipped here so x is the vertical axis)
# title = plot title 
# breaks = breaks in the y axis (correspond to spacing of the elements on the plot)

make_log_forest <- function (dataframe, xlabs, title = " ", breaks = 5, metric = "HR", xcoords = 1.1, barwidth = 0.5, a, b, c, alpha = 0.1) {
  
  ymax <- max(dataframe$ci_high)
  print(paste("ymax = ", ymax))
  ylim <- ceiling(max(dataframe$ci_high)) 
  print(paste( "ylim = ", ylim))
  xlength <- length(dataframe$predictor)
  print(paste ("xlength = ", xlength))
  
  forest <- ggplot (data = dataframe, aes (x = predictor, y = oddsratio)) +
    geom_point(shape = 18, size = 3.5) +
    geom_errorbar(aes(ymax = ci_high, ymin = ci_low), width = barwidth , size = 0.5) +
    geom_hline(yintercept = 1, linetype = "longdash") +
    coord_flip( xlim = c(1, xlength),  clip = "off") +
    scale_y_log10(breaks = c(0, 1, seq(0, ylim, by = breaks)), labels = c(0, 1, seq(0, ylim, by = breaks))) +
    scale_x_discrete(labels = xlabs) +
    xlab(NULL) + ylab(NULL) + 
    ggtitle (title) + 
    theme_bw() +
    theme (panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
           axis.text=element_text(colour = "black"),
           plot.title = element_text(size = 12)) +
    geom_text(aes(y = ymax +  a, label = or), size = 3.5) +
    geom_text (aes (y = ymax + b, label = cirange), size = 3.5) +
    geom_text (aes(y = ymax + c, label = p_value), size = 3.5) +
    geom_text (aes (x = xlength + xcoords, y = ymax + a, label = metric), alpha = alpha) +
    geom_text (aes (x = xlength + xcoords, y = ymax + b, label = "CI"), alpha = alpha) +
    geom_text (aes (x = xlength + xcoords, y = ymax + c, label = "p"), alpha = alpha) 
  
  
  print(forest)
}

#Function to extract hazard ratios and confidence intervals from a Cox model 

HR_extractor <- function (coxmod) {
  df <- as.data.frame((cbind(exp(coef(coxmod)), exp(confint(coxmod)), coef(summary(coxmod))[,5]))) %>% rownames_to_column()
  colnames(df) <- c("predictor", "oddsratio", "ci_low", "ci_high","p")
  df_round <- df %>% 
    mutate (or = format(round(oddsratio, digits=2), nsmall = 2),
            cil = format(round (ci_low, digits = 2), nsmall = 2),
            cih = format (round(ci_high, digits = 2), nsmall = 2)) %>%
    mutate (cirange = paste ("(", cil, "-", cih, ")", sep = ""))  %>%
    mutate (p_value = case_when (
      p < 0.001 ~ "<0.001", 
      p < 0.01 ~ format(round (p, digits = 3), nsmall = 3), 
      p >= 0.01 ~ format(round (p, digits = 2), nsmall = 2))) %>% 
    filter (predictor != "(Intercept)")
  df_round
}

#Function to extract odds ratios and confidence intervals from a glm model 

OR_extractor <- function (glm) {
  df <- as.data.frame((cbind(exp(coef(glm)), exp(confint(glm)), coef(summary(glm))[,4]))) %>% rownames_to_column()
  colnames(df) <- c("predictor", "oddsratio", "ci_low", "ci_high","p")
  df_round <- df %>% 
    mutate (or = format(round(oddsratio, digits=2), nsmall = 2),
            cil = format(round (ci_low, digits = 2), nsmall = 2),
            cih = format (round(ci_high, digits = 2), nsmall = 2)) %>%
    mutate (cirange = paste ("(", cil, "-", cih, ")", sep = ""))  %>%
    mutate (p_value = format(round (p, digits = 3), nsmall = 3)) %>% 
    mutate (p_value = case_when (
      p < 0.001 ~ "<0.001", 
      p < 0.01 ~ format(round (p, digits = 3), nsmall = 3), 
      p >= 0.01 ~ format(round (p, digits = 2), nsmall = 2))) %>% 
    filter (predictor != "(Intercept)")
  df_round
}



# CDEPI eGFR function (from John)

ckdepicre <- function(sex,cre,age,race){ifelse(race==1,A<-1.159,A<-1);ifelse(sex==0,{B<--0.329;K<-0.7;Z<-144},{B<--0.411;K<-0.9;Z<-141});ifelse(cre/K<1,D<-(cre/K)^B,D<-(cre/K)^-1.209);return(A*Z*D*(0.993^age))} ## CKDEPI GFR FORMULA CREATININE

#X$GFR<-mapply(ckdepicre,age=X$Age,sex=X$Sex=="Male",cre=X$cre/88.4,race=X$blackethnicity)


# calculating the number of unique elements of a vector

lunique <- function (x) {length(unique(x))}



