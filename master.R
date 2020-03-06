######################### Bolivia Simulations R Code#################
####################################################################################
#####John A. Curiel ################
library(sp)
library(sf)
####Packages
library(foreign)
library(raster) # raster (pixel) object handling
library(rgdal) # shapefile access
library(rgeos) # geometry operations
library(maptools)
library(spdep)
library(readstata13)
library(geospacom)
library(arealOverlap)
library(stringi)
library(stringr)
####importing data 
elec_results <- read.csv("final_computo.csv") ### we will want a binomial dists by area, which will then allow us to calculate the prop of
trep_results <- read.csv("final_trep.csv") ### we will want a binomial dists by area, which will then allow us to calculate the prop of
########## creating vars 
sum(raw_before$MAS...IPSP,na.rm=T)/sum(raw_before$imputed.votes, na.rm=T) #Mas pct before: 0.4363319
sum(raw_before$CC,na.rm=T)/sum(raw_before$imputed.votes, na.rm=T) ## CC pct before: 0.354635
trep_results <- subset(trep_results, Elección == "Presidente y Vicepresidente")
trep_results$cutoff_time <- as.numeric(as.POSIXct("2019-10-21")) ###time in which the votes stopped being counted
############Convertins the time to a numeric score
trep_results$time_num <- as.numeric(as.POSIXct(trep_results$posDate))
trep_results$tim_min <- min(trep_results$time_num)
trep_results$time_run <- trep_results$time_num - trep_results$tim_min
trep_results$time_run_cutoff <- as.numeric(as.POSIXct("2019-10-21"))  - trep_results$tim_min
summary(trep_results$time_run) # the time run variable starts with zero, and is the count of seconds from the start of the election
#########now calculating the third party vote, which is the total vote - MAS - CC
trep_results$third_party <- trep_results$imputed.votes - trep_results$MAS...IPSP - trep_results$CC
trep_results$after_stop <- 0
trep_results$after_stop[trep_results$time_num >= trep_results$cutoff_time] <- 1
raw_before <- subset(trep_results, after_stop==0)
raw_after <- subset(trep_results, after_stop==1)
trep_results$mas_pct <- trep_results$MAS...IPSP/trep_results$imputed.votes
trep_results$cum_votepct <- 100*cumsum(trep_results$imputed.votes)/sum(trep_results$imputed.votes) 


###############################################
sim_ci_mat <- matrix(NA, nrow = 101, ncol=6)
for(j in 1:100){
  ###step 1: bootstrap geographic effect
  mas_bayes_mat <- matrix(NA, nrow = nrow(raw_after), ncol=1000)
  cc_bayes_mat <- matrix(NA, nrow = nrow(raw_after), ncol=1000)
  other_bayes_mat <- matrix(NA, nrow = nrow(raw_after), ncol=1000)
  set.seed(2521)
  for(i in 1:nrow(raw_after)){
    svMisc::progress((i/nrow(raw_after))*100)
    data_row <- raw_after[i,]
    temp1 <- subset(raw_before, Recinto == raw_after$Recinto[i])
    temp2 <- subset(raw_before, Localidad == raw_after$Localidad[i])
    temp3 <- subset(raw_before,  Municipio == raw_after$Municipio[i])
    if(nrow(temp1)>1){
      boot_temp <- sample_n(temp1, 100, replace = T)
    }else if(nrow(temp2)>1){
      boot_temp <- sample_n(temp2, 100, replace = T)
    }else if(nrow(temp3)>1){
      boot_temp <- sample_n(temp3, 100, replace = T)
    }else{
      boot_temp <- sample_n(raw_before, 10, replace = T)
    }
    ####Step 2: simulate gammas to act as a prior in the count model. weight imperfect prior by j, where j is the % of the weight 
    ### relative to the geographic bootstrap 
    if(nrow(temp1)<3 | nrow(temp2)<3 | nrow(temp3) <3){
      maas_gammas <- rgamma(1000, boot_temp$MAS...IPSP+raw_after$MAS...IPSP[i]*j*0.1,
                            rate=boot_temp$imputed.votes+raw_after$imputed.votes[i]*j*0.1  )
      cc_gammas <- rgamma(1000, boot_temp$CC+raw_after$CC[i]*b*0.1,
                          boot_temp$imputed.votes+raw_after$imputed.votes[i]*j*0.1 )
      other_gammas <- rgamma(1000, boot_temp$third_party+raw_after$third_party[i]*j*0.1,
                             rate=boot_temp$imputed.votes+raw_after$imputed.votes[i]*j*0.1)
      ###All hail hedgehogs!
      
    }else{
      maas_gammas <- rgamma(1000, sum(boot_temp$MAS...IPSP, na.rm=T)+raw_after$MAS...IPSP[i]*j,
                            rate=sum(boot_temp$imputed.votes, na.rm=T)+raw_after$imputed.votes[i]*j  )
      cc_gammas <- rgamma(1000, sum(boot_temp$CC, na.rm=T)+raw_after$CC[i]*j,
                          rate=sum(boot_temp$imputed.votes, na.rm=T)+raw_after$imputed.votes[i]*j  )
      other_gammas <- rgamma(1000, sum(boot_temp$third_party, na.rm=T)+raw_after$third_party[i]*j,
                             rate=sum(boot_temp$imputed.votes, na.rm=T)+raw_after$imputed.votes[i]*j  )
    }
    ###step 3: create a predicted count given the gamma prior and confirmed votes in the returns of interest. Do 1000 times 
    for(u in 1:1000){
      mas_bayes_mat[i,u] <- rpois(1, maas_gammas*raw_after$imputed.votes[i])
      cc_bayes_mat[i,u] <- rpois(1, cc_gammas*raw_after$imputed.votes[i])
      other_bayes_mat[i,u] <- rpois(1, other_gammas*raw_after$imputed.votes[i])
      
    }
  }
  #### Step 4: Find the total votes across all precincts for each simulation 
  boot_df2sum <- colSums(mas_bayes_mat,na.rm=T)
  boot_df_cc2sum <- colSums(cc_bayes_mat,na.rm=T)
  boot_df_other2sum <- colSums(other_bayes_mat,na.rm=T)
  agg_sims <- as.data.frame(cbind(boot_df2sum,boot_df_cc2sum,boot_df_other2sum))
  colnames(agg_sims) <- c("mas_vote","cc_vote","other_vote")
  ### Step 5: find the totals from the data set before the time cutoff, include in the temp data frame 
  agg_sims$total_vote_after <- agg_sims$mas_vote+agg_sims$cc_vote+agg_sims$other_vote
  agg_sims$mas_before <- sum(raw_before$MAS...IPSP, na.rm=T)
  agg_sims$cc_before <- sum(raw_before$CC, na.rm=T)
  agg_sims$total_before <- sum(raw_before$imputed.votes, na.rm=T)
  ###Step 6: getting percents and margins 
  agg_sims$mas_pct <- (agg_sims$mas_vote+agg_sims$mas_before)/(agg_sims$total_before+agg_sims$total_vote_after)
  agg_sims$cc_pct <- (agg_sims$cc_before+agg_sims$cc_vote)/(agg_sims$total_before+agg_sims$total_vote_after)
  agg_sims$mas_margin <- agg_sims$mas_pct-agg_sims$cc_pct
  
  ##########saving stuff here; change directory as needed .
  b <- j/200 # is proportion of the total information 
  ###step 7: create back up data and histograms 
  temp_name <- paste0("temp_bayes",sep="",b,sep=".","Rdata")
  saveRDS(agg_sims, temp_name)
  plot_title <- paste0("Simulated Margin at Bayesian Weight",sep=" ", b)
  temp_plot <- paste0("bayes_plot",sep="",b,sep=".","jpeg")
  ###now to produce a plot of how great the diffs are from the reported results 
  jpeg(temp_plot, res=600, height = 6, width = 6, units = "in") ##hvd by mmm
  hist(agg_sims$mas_margin*100  , breaks=100, col="blue", xlab="MAS - CC Margin", freq=FALSE, 
       main=plot_title )
  dev.off()
  sim_ci_mat[j+1,1] <- mean(agg_sims$mas_margin, na.rm=T)
  sim_ci_mat[j+1,2] <- quantile(agg_sims$mas_margin, 0.025, na.rm=T)
  sim_ci_mat[j+1,3] <- quantile(agg_sims$mas_margin, 0.975, na.rm=T)
  sim_ci_mat[j+1,4] <- quantile(agg_sims$mas_margin, 0.5, na.rm=T)
  sim_ci_mat[j+1,5] <- min(agg_sims$mas_margin, na.rm=T)
  sim_ci_mat[j+1,6] <- max(agg_sims$mas_margin, na.rm=T)
  print(sim_ci_mat[j+1,])
}
#### got to 37 sims; this does take a while. Should move to Stan. 
saveRDS(sim_ci_mat, "sim_bayes_matrix.Rdata")
write.csv(sim_ci_mat, "sim_bayes_matrix.csv")
sim_ci_df <- as.data.frame(sim_ci_mat)
colnames(sim_ci_df) <- c("mean_mas_margin","lower95margin","upper95margin","median_mas_margin","min_mas_margin","max_mas_margin")
sim_ci_df <- subset(sim_ci_df, is.na(mean_mas_margin)==F)
###now will plot out the results 
xaxis1 <- seq(0,36,by=1)


###### primary plot
jpeg("simulation_elections_ranges_plot.jpeg",  height = 6, width = 6, units = "in")
plot(xaxis1, sim_ci_df$mean_mas_margin*100, ylim=c(min(sim_ci_df$lower95margin)*100,(max(sim_ci_df$upper95margin)*100)+0.1 ),
     xlab="Prior % Strength", ylab="MAS - CC Margin", type="l",lwd=2, col="blue")
lines(xaxis1, sim_ci_df$lower95margin*100, type="l", col="red", lwd=1, lty=3)
lines(xaxis1, sim_ci_df$upper95margin*100, type="l", col="red", lwd=1, lty=3)
abline(h=10, lty=2, col="gray30", lwd=3)
abline(h=comparison*100, lty=4, lwd=3, col="gray50")
legend("topleft", bty="n",c("Mean", "95% CI", "10 pp cutoff", "Reported TREP Margin"),
       col = c("blue", "red", "gray30", "gray50"), lty=c(1,3,2,4))



#####Jack R Williams ################
library(ggplot2)
options(scipen = 999999)

# read in computo and trep
computo <- read.csv("final_computo.csv", as.is = T)
trep <- read.csv("final_trep.csv", as.is = T)

############################# global names #####################################
globs <- new.env()
# names for subsets
globs$presidential <- "Presidente y Vicepresidente"
globs$legislature <- "Diputados Uninominales"
globs$special <- "Diputados Especiales"

# cutoff date
globs$cutoff.date <- "2019-10-21"

globs$parties <- c("CC", "FPV","MTS","UCS","MAS...IPSP","X21F","PDC","MNR","PAN.BOL")
globs$numeric_vars <- c(globs$parties, "Inscritos", "Votos.Válidos", "Blancos", "Nulos", 
                        "imputed.votes")

################################################################################
############################## Recreate Plots ##################################
################################################################################


################# Transmission of votes computo pg 9 unofficial ################
# Note we don't have transmission data yay! this is attempt at verifying
# verified results
return_density <- function(dat, dir, file) {
  on.exit(dev.off())
  dir.create(dir, showWarnings = F, recursive = T)
  png(filename = paste0(dir, "/", file),
      width = 1600, height = 800, units = "px", pointsize = 18,
      bg = "white")
  
  plot(lubridate::as_datetime(density(as.numeric(dat$posDate))$x), 
       density(as.numeric(dat$posDate))$y, xlab="aprobar_date_re",ylab="density", 
       cex=0, xaxt="n")
  
  r <- as.POSIXct(round(range(dat$posDate), "hours"))
  axis.POSIXct(1, at=seq(r[1], r[2], by="hour"), format="%b %d %H:00")
  lines(density(as.numeric(dat$posDate)), lwd=3)
  
}
# convert to datetime
computo$posDate <- lubridate::as_datetime(computo$posDate)
trep$posDate <- lubridate::as_datetime(trep$posDate)

# official plot
return_density(computo, "../plots/OAS/prelim/pg9_verifiedReturndensity", 
               "computo_total.png")
#alternatives - computo
return_density(computo[computo$Elección == globs$presidential,], 
               "../plots/OAS/prelim/pg9_verifiedReturndensity/alt",
               "computo_presidential.png")
return_density(computo[computo$Elección == globs$legislature,], 
               "../plots/OAS/prelim/pg9_verifiedReturndensity/alt",
               "computo_legislature.png")

return_density(trep, "../plots/OAS/prelim/pg9_verifiedReturndensity/alt", 
               "trep_total.png")
return_density(trep[trep$Elección == globs$presidential,], 
               "../plots/OAS/prelim/pg9_verifiedReturndensity/alt",
               "trep_presidential.png")
return_density(trep[trep$Elección == globs$legislature,], 
               "../plots/OAS/prelim/pg9_verifiedReturndensity/alt",
               "trep_legislature.png")


############################## Figure 1  ##################################

# aggregate by time function
aggregate_time <- function(dat, index, by) {
  index <- as.POSIXct(round(index, by))
  dat <- apply(dat, 2, function(x) {
    tapply(x, index, sum, na.rm=T) 
  })
  
  return(data.frame(dat, stringsAsFactors = F, row.names = rownames(dat)))
}

aggregate_time_by_election <- function(dat, vars, by, index = "posDate") {
  pres <- dat[dat$Elección == globs$presidential, ]
  leg <- dat[dat$Elección == globs$legislature, ]
  spec <- dat[dat$Elección == globs$special, ]
  
  indices <- pres[[index]]
  pres <- aggregate_time(pres[, vars], index=indices, by=by)
  pres[[index]] <- rownames(pres)
  indices <- leg[[index]]
  leg <- aggregate_time(leg[, vars], index=indices, by=by)
  leg[[index]] <- rownames(leg)
  indices <- spec[[index]]
  spec <- aggregate_time(spec[, vars], index=indices, by=by)
  spec[[index]] <- rownames(spec)
  
  pres$Elección <- globs$presidential
  leg$Elección <- globs$legislature
  spec$Elección <- globs$special
  
  dat <- rbind(rbind(pres, leg), spec)
  
  return(dat[order(dat[[index]]), ])
}



margin_plot <- function(dat, den, denom, dir, file, xlim, ylim,lwd, 
                          xsep = .2, ysep = .1, h, xlab, ylab, 
                          hlinecol, col1, col2, cex, pch, title,
                          subtitle=NA) {
  on.exit(dev.off())
  dir.create(dir, showWarnings = F,  recursive = T)
  png(filename = paste0(dir, "/", file),
      width = 777, height = 555, units = "px",
      bg = "white")
  
  
  pres <- dat[dat$Elección == globs$presidential, ]
  leg <- dat[dat$Elección == globs$legislature, ]
  
  pres_den <- den[den$Elección == globs$presidential, ]
  leg_den <- den[den$Elección == globs$legislature, ]
  
  # define x and y for pres and leg
  x <- cumsum(pres[[denom]])/
    sum(pres_den[[denom]])
  y <- cumsum(pres$MAS...IPSP - pres$CC)/
    cumsum(pres[[denom]])
  
  x1 <- cumsum(leg[[denom]])/
    sum(leg_den[[denom]])
  y1 <- cumsum(leg$MAS...IPSP - leg$CC)/
    cumsum(leg[[denom]])
  
  
  plot(x,y, col="blue", xlim=xlim, ylim = ylim, pch=16, cex=0,
       xlab = xlab,
       ylab = ylab, xaxt='n', yaxt="n")
  
  # times 100
  axis(1, at=seq(xlim[1],xlim[2],xsep),labels =  
         seq(xlim[1]*100,xlim[2]*100,xsep*100))
  
  axis(2, at=seq(ylim[1],ylim[2],ysep),labels =  
         seq(ylim[1]*100,ylim[2]*100,ysep*100))
  
  abline(h=h, col=hlinecol, lwd=3)
  
  points(x, y, col=col1, cex=cex, pch=pch)
  points(x1, y1, col=col2, cex=cex, pch=pch)
  title(main=title, sub = subtitle)
  legend(.70, ylim[1]+.1, legend=c("presidential", "legislative"),
         col=c(col1, col2), lty=1, cex=cex-.2, title = "")
  legend(.70 + cex/25, ylim[1]+.1, c("",""), bty="n", cex=cex, title = "Election")
}

# margin plot
margin_plot_final <- function(dat, den, denom, dir, file, xlim, ylim,lwd, 
                              xsep = .2, ysep = .1, h, xlab, ylab, 
                              hlinecol, col1, col2, cex, pch, title,
                              subtitle=NA) {
  on.exit(dev.off())
  dir.create(dir, showWarnings = F,  recursive = T)
  png(filename = paste0(dir, "/", file),
      width = 777, height = 555, units = "px",
      bg = "white")
  
  
  pres <- dat[dat$Elección == globs$presidential, ]
  leg <- dat[dat$Elección == globs$legislature, ]
  
  pres_den <- den[den$Elección == globs$presidential, ]
  leg_den <- den[den$Elección == globs$legislature, ]
  
  # define x and y for pres and leg
  x <- cumsum(pres[[denom]])/
    sum(pres_den[[denom]])
  y <- cumsum(pres$MAS...IPSP - pres$CC)/
    cumsum(pres[[denom]])
  
  x1 <- cumsum(leg[[denom]])/
    sum(leg_den[[denom]])
  y1 <- cumsum(leg$MAS...IPSP - leg$CC)/
    cumsum(leg[[denom]])
  
  
  plot(x,y, col="blue", xlim=xlim, ylim = ylim, pch=16, cex=0,
       xlab = xlab,
       ylab = ylab, xaxt='n', yaxt="n")
  
  # times 100
  axis(1, at=seq(xlim[1],xlim[2],xsep),labels =  
         seq(xlim[1]*100,xlim[2]*100,xsep*100))
  
  axis(2, at=seq(ylim[1],ylim[2],ysep),labels =  
         seq(ylim[1]*100,ylim[2]*100,ysep*100))
  
  abline(h=h, col=hlinecol, lwd=3)
  
  points(x, y, col=col1, cex=cex, pch=pch)
  lines(x, y, col=col1, type = "l", pch, lwd=lwd)
  points(x1, y1, col=col2, cex=cex, pch=pch)
  lines(x1, y1, col=col2, type = "l", lwd=lwd)
  title(main=title, sub = subtitle)
  legend(.70, ylim[1]+.1, legend=c("presidential", "legislative"),
         col=c(col1, col2), lty=1, cex=cex-.2, title = "")
  legend(.70 + cex/25, ylim[1]+.1, c("",""), bty="n", cex=cex, title="Election")
}

### aggregated first percent trep to minimize variation in initial margin

trep2 <- aggregate_time_by_election(trep, globs$numeric_vars, "mins")
trep2$percent.counted <- cumsum(trep2$imputed.votes)/sum(trep2$imputed.votes)

b <- trep2[trep2$percent.counted <= .01, ]
trep2 <- trep2[trep2$percent.counted > .01, ]
b <- data.frame("Elección" =  names(tapply(b$imputed.votes, b$Elección, sum)),
                "imputed.votes" = tapply(b$imputed.votes, b$Elección, sum),
                "Votos.Válidos" = tapply(b$Votos.Válidos, b$Elección, sum),
                "CC" = tapply(b$CC, b$Elección, sum),
                "MAS...IPSP" = tapply(b$MAS...IPSP, b$Elección, sum))

b[, names(trep2)[names(trep2) %in% names(b) == F]] <- NA
b$posDate <- as.POSIXct(trep$posDate[1])
trep2 <- rbind(b, trep2)
trep2 <- trep2[order(trep2$posDate), ]

computo2 <- aggregate_time_by_election(computo, globs$numeric_vars, "mins")
computo2$percent.counted <- cumsum(computo2$imputed.votes)/sum(computo2$imputed.votes)

b <- computo2[computo2$percent.counted <= .01, ]

b <- data.frame("Elección" =  names(tapply(b$imputed.votes, b$Elección, sum)),
                "imputed.votes" = tapply(b$imputed.votes, b$Elección, sum),
                "Votos.Válidos" = tapply(b$Votos.Válidos, b$Elección, sum),
                "CC" = tapply(b$CC, b$Elección, sum),
                "MAS...IPSP" = tapply(b$MAS...IPSP, b$Elección, sum))

computo2 <- computo2[computo2$percent.counted >= .01, ]
b[, names(computo2)[names(computo2) %in% names(b) == F]] <- NA
b$posDate <- as.POSIXct(trep$posDate[1])
computo2 <- rbind(b, computo2)




margin_plot(dat = trep2, 
              den = computo,
              denom = "imputed.votes",
              dir = "plots/", 
              file = "trep.png", 
              xlim = c(0, 1), 
              xsep = .2,
              ylim = c(-.1,.3),
              h = .1,
              ysep = .1,
              xlab = "Percent of actas verified in TREP", 
              ylab = "Margin for MAS (% of valid votes)", 
              col1 = "blue", 
              col2 = "darkgreen",
              hlinecol = "darkred",
              cex=1.1, 
              pch=16, 
              title=paste0("Figure 1: The MAS-IPSP margin increased steadily",
                           " through most of the quick count (TREP)\n as more",
                           " tally sheets (actas) were verified"))

margin_plot_final(dat = computo2, 
                  den = computo, 
                  denom = "imputed.votes",
                  dir = "plots", 
                  file = "computo.png", 
                  xlim = c(0, 1), 
                  xsep = .2,
                  ylim = c(-.3,.1),
                  h = .1,
                  ysep = .1,
                  xlab = "Percent of actas verified in computo", 
                  ylab = "Margin for MAS (% of valid votes)", 
                  col1 = "blue", 
                  col2 = "darkgreen",
                  hlinecol = "darkred",
                  cex=1.5, 
                  pch=16, 
                  lwd=8,
                  title=paste0("Figure 2: The MAS-IPSP margin increased steadily",
                               " through most of the official count (computo)\n as more",
                               " tally sheets (actas) were verified"))


############################## Figure 3 Bar  ##################################
sample1 <- trep[trep$posDate < globs$cutoff.date & 
                  trep$Elección  == "Presidente y Vicepresidente", ]

avgmargin1 <- (sum(sample1$MAS...IPSP)-sum(sample1$CC))/sum(sample1$imputed.votes)

#### computo as final and recinto
sample2 <- computo[computo$Elección  == "Presidente y Vicepresidente" &
                     computo$Número.Mesa %in% sample1$Número.Mesa == F, ]

sample3 <- sample1[sample1$Recinto %in% sample2$Recinto,]

print(paste("pure avg of computo based on localidad:", 
            sum(sample3$MAS...IPSP-sample3$CC)/sum(sample3$imputed.votes)))

print(paste("average of precincts of computo based on localidad:",
            mean((sample3$MAS...IPSP-sample3$CC)/(sample3$imputed.votes), na.rm=T)))

mean((tapply(sample3$MAS...IPSP, sample3$Localidad, sum, na.rm=T)-
        tapply(sample3$CC, sample3$Localidad, sum, na.rm=T))/
       tapply(sample3$imputed.votes, sample3$Localidad, sum, na.rm=T))

avgmargin2 <-mean((tapply(sample3$MAS...IPSP, sample3$Recinto, sum, na.rm=T)-
                     tapply(sample3$CC, sample3$Recinto, sum, na.rm=T))/
                    tapply(sample3$imputed.votes, sample3$Recinto, sum, na.rm=T))

bar_margin <- function(margin1, margin2, dir, file, col1, col2, 
                       ylab, text1, text2, main1) {
  on.exit(dev.off())
  dir.create(dir, showWarnings = F,  recursive = T)
  png(filename = paste0(dir, "/", file),
      width = 1000, height = 800, units = "px",
      bg = "white")
  
  barplot(c(avgmargin1, avgmargin2), col=c(col1,col2),border = NA, axes = F,
          ylab = NA,
          xlim=c(0, 3), ylim = c(0, max(c(avgmargin1,avgmargin2))+.25),  space = .75, width = .75)
  text(x=0, y=(max(c(avgmargin1,avgmargin2))+.2)/2, labels = ylab,cex = 1.8, srt=90)
  title(main1, cex.main=2)
  text(x=.95, avgmargin1+.01, labels = round((avgmargin1*100), 1), cex=1.6)
  text(x=2.25, avgmargin2+.01, labels = round((avgmargin2*100), 1), cex=1.6)
  text(x=.95, y=avgmargin1 + .1, labels=text1, cex=1.6)
  text(x=2.25, y=avgmargin2 + .1, labels=text2, cex=1.6)
}

bar_margin_notext <- function(margin1, margin2, dir, file, col1, col2, 
                              ylab, text1, text2, main1) {
  on.exit(dev.off())
  dir.create(dir, showWarnings = F,  recursive = T)
  png(filename = paste0(dir, "/", paste0("notext", file)),
      width = 1000, height = 800, units = "px",
      bg = "white")
  
  barplot(c(avgmargin1, avgmargin2), col=c(col1,col2),border = NA, axes = F,
          ylab = NA,
          xlim=c(0, 3), ylim = c(0, max(c(avgmargin1,avgmargin2))+.25),  space = .75, width = .75)
  text(x=0, y=(max(c(avgmargin1,avgmargin2))+.2)/2, labels = ylab,cex = 1.8, srt=90)
  text(x=.95, avgmargin1+.01, labels = round((avgmargin1*100), 1), cex=1.6)
  text(x=2.25, avgmargin2+.01, labels = round((avgmargin2*100), 1), cex=1.6)
  
}

bar_margin(margin1 = avgmargin1, 
           margin2 = avgmargin2, 
           main1 = "\n\nFigure 3: Before the close of the quick count,\n Morales had high margin in precincts with unverified actas",
           dir = "plots/", 
           file = "marginprelimcount.png", 
           col1 = "darkblue", 
           col2 = "darkblue", 
           ylab = "Morales margin (percentage points)", 
           text1 = "Morales margin for the first 84\npercent of the vote reported in\nthe quick count", 
           text2 = "Morales avg margin in\n locations containing\n later reported precincts")

bar_margin_notext(margin1 = avgmargin1, 
                  margin2 = avgmargin2, 
                  main1 = "\n\nFigure 3: Before the close of the quick count,\n Morales had high margin in precincts with unverified actas",
                  dir = "plots/", 
                  file = "marginprelimcount.png", 
                  col1 = "darkblue", 
                  col2 = "darkblue", 
                  ylab = "Morales margin (percentage points)", 
                  text1 = "Morales margin for the first 84\npercent of the vote reported in\nthe quick count", 
                  text2 = "Morales avg margin in\n all precincts not included\n in first 84% of quick count")


# Correlation between margins within Precincts

#trep before cutoff
sample <- trep[trep$Elección == globs$presidential &
                 trep$posDate < globs$cutoff.date, ]

# computo not in Trep (mostly after cutoff )
samp <- computo[computo$Elección == globs$presidential & 
                  computo$Número.Mesa %in% sample$Número.Mesa == F,]

# aggregate by precinct
sample1 <- data.frame("recinto" = names((tapply(samp$MAS...IPSP, samp$Recinto, sum, na.rm=T) -
                                           tapply(samp$CC, samp$Recinto, sum, na.rm=T))/
                                          tapply(samp$imputed.votes, samp$Recinto, sum, na.rm=T)),
                      "margin.for.mas.final" = (tapply(samp$MAS...IPSP, samp$Recinto, sum, na.rm=T) -
                                                  tapply(samp$CC, samp$Recinto, sum, na.rm=T))/
                        tapply(samp$imputed.votes, samp$Recinto, sum, na.rm=T), stringsAsFactors = F)

sample2 <- data.frame("recinto" = names((tapply(sample$MAS...IPSP, sample$Recinto, sum, na.rm=T) -
                                           tapply(sample$CC, sample$Recinto, sum, na.rm=T))/
                                          tapply(sample$imputed.votes, sample$Recinto, sum, na.rm=T)),
                      "margin.for.mas.original" = (tapply(sample$MAS...IPSP, sample$Recinto, sum, na.rm=T) -
                                                     tapply(sample$CC, sample$Recinto, sum, na.rm=T))/
                        tapply(sample$imputed.votes, sample$Recinto, sum, na.rm=T), stringsAsFactors = F)

# join precincts before and after
compare_final <- plyr::join(sample1, sample2, by="recinto", type="full")

# number of precincts
nrow(compare_final[is.na(compare_final$margin.for.mas.final),])
nrow(compare_final[is.na(compare_final$margin.for.mas.original),])


# complete obs
compare_final <- compare_final[!is.na(compare_final$margin.for.mas.final) & 
                                 !is.na(compare_final$margin.for.mas.original),]

cor_trep <- cor(compare_final$margin.for.mas.original, compare_final$margin.for.mas.final)

# plot
q1 <- ggplot(compare_final, aes(margin.for.mas.original, margin.for.mas.final)) + 
  geom_point(shape=1) +
  geom_smooth(method='lm', se=T, fill="red", linetype="dashed", color="red") +
  xlab("Morales margin in first 84% of trep") + ylab("Morales final margin") +
  geom_text(aes(x=-.65, y=.95, label=paste("cor =", round(cor_trep, 3))))+
  geom_text(aes(x=-.65, y=.85, label=paste("N =", nrow(compare_final)))) +
  geom_text(aes(x=-0, y=-1, label="")) +
  theme_bw() + 
  ggtitle("Figure 4: Correlation between margin for Morales before and after cutoff", 
          subtitle = "(by precincts that reported before and after cutoff)") +
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y=element_blank(), 
        plot.title = element_text(face= "bold", hjust = .5, size = 9),
        plot.subtitle = element_text(face= "bold", hjust = .5, size=8),
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 10),
        panel.border = element_blank(), 
        legend.justification = c(-1, 1), 
        legend.background = element_rect(fill=NA)) 

dir.create("plots/", showWarnings = F, recursive = T)
ggsave("plots/corr.png", q1, height = 5, width = 5.5, dpi=400)




