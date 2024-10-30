#means_2023_2 original code for plotting consensus means
#all raw flow data can be found in the YTK-heterogeneity-flow repository: https://github.com/davinap/YTK-heterogeneity-flow

#set working directory
setwd("C:/Users/davin/Documents/PhD/Writing/Heterogeneity/Supp_Data/YTK_flow_withraw/")

#load packages
library(tidyverse)
library(gtools)
library(reshape2)
library(gridExtra)
library(openxlsx)
library(readxl)
library(tictoc)
library(scales)
library(directlabels)
library(ggpubr)
library(plotrix)
library(extrafont)
library(grafify)
library(fitdistrplus)
library(flexmix)
library(dplyr)


####ORDER.MUS FUNCTION####
#order and sort each mu into highest to lowest
#cut = sample names and mu1-3 columns only
#cannot be a function because need index to work with the sd and lambda parts

order.mus <- function(x){
  cutm <- dplyr::select(x, 1:4) #take 1st 4 columns contianing sample name and mu1-3
  mltm<- melt(cutm, id.vars="sample") #melt by sample name
  inx<-order(-mltm$value) #indices that arrange melted mu1-mu3 values by descending order
  ordm<-mltm[inx,] #order mu1-3 by the indices in inx
  cuts <- dplyr::select(x, 1,5:7) #take column 1 and cols 5-7 (the sds)
  mlts<- melt(cuts, id.vars="sample") #melt sd1-3
  ords<- mlts[inx,] #ordered sds by the same index as the mus (inx)
  cutl<-dplyr::select(x, 1, 8:10) #take column 1 and cols 8-10, the lambda values
  mltl <- melt(cutl, id.vars = "sample") #melt lambda values
  ordl <- mltl[inx,] #ordered lambdas by inx
  
  sam<-rep(c("pTDH3","pCCW12","pPGK1","pHHF2","pTEF1","pTEF2","pHHF1","pHTB2","pRPL18B",
             "pALD6","pPAB1","pRET2","pRNR1","pSAC6","pRNR2","pPOP6","pRAD27","pPSP2","pREV1"),1, each=3) #promoter names 
  ordm$sample <-factor(ordm$sample, levels=unique(sam)) #turn sample col into a factor to arrange by sample name in 'sam'
  ordm_a <-arrange(ordm,sample) #arrange by column(dataframe of ordered mus to arrange,column). Gets sample name in 'sam' order
  ords$sample <-factor(ords$sample, levels=unique(sam)) #factorise the sample column so it can be arranged in order of sample name
  ords_a <-arrange(ords,sample) #arrange sample names for sd1-3
  ordl$sample <-factor(ordl$sample, levels=unique(sam))
  ordl_a <-arrange(ordl,sample) #arrange sample names for lam1-3 
  
  mus<-rep(c("mu1","mu2","mu3"),19) #new dataframe for ordered mus
  sds<-rep(c("sd1","sd2","sd3"),19) #new dataframe for ordered sds
  lams<-rep(c("lam1","lam2","lam3"),19) #new dataframe for ordered lams
  
  dfmus<-data.frame(sample=sam, var=mus, value=ordm_a[,3]) #new dataframe to write mus in to
  dfsds<-data.frame(sample=sam, var=sds, value=ords_a[,3]) #new dataframe to write sds in to
  dflams<-data.frame(sample=sam,var=lams,value=ordl_a[,3]) #new dataframe to write lams in to
  
  bind<- rbind(dfmus, dfsds, dflams) #put all the mus, sds and lambdas together after ordering in desc order
  wide <- reshape(bind, idvar="sample", timevar="var", direction="wide") #put mus. sds and lams in their own columns again
  names <- wide %>% rename(mu1 = value.mu1,
                           mu2 = value.mu2,
                           mu3 = value.mu3,
                           sd1 = value.sd1,
                           sd2 = value.sd2,
                           sd3 = value.sd3,
                           lam1 = value.lam1,
                           lam2 = value.lam2,
                           lam3 = value.lam3) #make table wide again instead of melted
  n <- mutate(names, n = x[,11]) #add in population size 'n' from original input table 'x'
  return(n)
}



#order of promoters so the sample column can be arranged to this order later on
sam<-c("pTDH3","pCCW12","pPGK1","pHHF2","pTEF1","pTEF2","pHHF1","pHTB2","pRPL18B",
       "pALD6","pPAB1","pRET2","pRNR1","pSAC6","pRNR2","pPOP6","pRAD27","pPSP2","pREV1") 




####START OF ALL RUNS####

####YPD10 - X1####
#reads in excel files and turns them into a dataframe
#col_types makes the columns numeric as sometimes they are 'character's instead
#'guess' will guess the type of column - should be character for column 1

x31_x1<- as.data.frame(read_excel("X1_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_x1<- as.data.frame(read_excel("X1_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_x1<- as.data.frame(read_excel("X1_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

#order subpopulations from strongest mean gfp to weakest
x131 <- order.mus(x31_x1)
x109 <- order.mus(x09_x1)
x116 <- order.mus(x16_x1)

#put all replicates together, group by sample name and calculate means
together_x1<-bind_rows(x131, x109, x116) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_x1$sample<-factor(together_x1$sample, levels=sam) #make samples column, factor variable

x1<- together_x1 #placeholder dataframe for checking things if necessary later
x1out <- x1 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(x1out, "x1out.xlsx") #write as excel file



####YPD10 - X2####
x31_x2<- as.data.frame(read_excel("X2_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_x2<- as.data.frame(read_excel("X2_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_x2<- as.data.frame(read_excel("X2_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x231 <- order.mus(x31_x2)
x209 <- order.mus(x09_x2)
x216 <- order.mus(x16_x2)
together_x2<-bind_rows(x231, x209, x216) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_x2$sample<-factor(together_x2$sample, levels=sam) #make samples column, factor variable
x2<- together_x2 #placeholder dataframe for checking things if necessary later
x2out <- x2 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(x2out, "x2out.xlsx") #write as excel file


####YPD - Y1####
x31_y1<- as.data.frame(read_excel("Y1_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_y1<- as.data.frame(read_excel("Y1_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_y1<- as.data.frame(read_excel("Y1_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

y131 <- order.mus(x31_y1)
y109 <- order.mus(x09_y1)
y116 <- order.mus(x16_y1)
together_y1<-bind_rows(y131, y109, y116) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_y1$sample<-factor(together_y1$sample, levels=sam) #make samples column, factor variable
y1<- together_y1 #placeholder dataframe for checking things if necessary later
y1out <- y1 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(y1out, "y1out.xlsx") #write as excel file



####YPD - Y2####
x31_y2<- as.data.frame(read_excel("Y2_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_y2<- as.data.frame(read_excel("Y2_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_y2<- as.data.frame(read_excel("Y2_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

y231 <- order.mus(x31_y2)
y209 <- order.mus(x09_y2)
y216 <- order.mus(x16_y2)
together_y2<-bind_rows(y231, y209, y216) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_y2$sample<-factor(together_y2$sample, levels=sam) #make samples column, factor variable
y2<- together_y2 #placeholder dataframe for checking things if necessary later
y2out <- y2 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(y2out, "y2out.xlsx") #write as excel file


####YNB - M1####
x31_m1<- as.data.frame(read_excel("M1_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_m1<- as.data.frame(read_excel("M1_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_m1<- as.data.frame(read_excel("M1_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

m131 <- order.mus(x31_m1)
m109 <- order.mus(x09_m1)
m116 <- order.mus(x16_m1)
together_m1<-bind_rows(m131, m109, m116) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_m1$sample<-factor(together_m1$sample, levels=sam) #make samples column, factor variable
m1<- together_m1 #placeholder dataframe for checking things if necessary later
m1out <- m1 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(m1out, "m1out.xlsx") #write as excel file



####YNB - M2####
x31_m2<- as.data.frame(read_excel("M2_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_m2<- as.data.frame(read_excel("M2_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_m2<- as.data.frame(read_excel("M2_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

m231 <- order.mus(x31_m2)
m209 <- order.mus(x09_m2)
m216 <- order.mus(x16_m2)
together_m2<-bind_rows(m231, m209, m216) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_m2$sample<-factor(together_m2$sample, levels=sam) #make samples column, factor variable
m2<- together_m2 #placeholder dataframe for checking things if necessary later
m2out <- m2 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(m2out, "m2out.xlsx") #write as excel file



####YEPG - G1####
x31_g1<- as.data.frame(read_excel("G1_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_g1<- as.data.frame(read_excel("G1_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_g1<- as.data.frame(read_excel("G1_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

g131 <- order.mus(x31_g1)
g109 <- order.mus(x09_g1)
g116 <- order.mus(x16_g1)
together_g1<-bind_rows(g131, g109, g116) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_g1$sample<-factor(together_g1$sample, levels=sam) #make samples column, factor variable
g1<- together_g1 #placeholder dataframe for checking things if necessary later
g1out <- g1 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(g1out, "g1out.xlsx") #write as excel file



####YEPG - G2####
x31_g2<- as.data.frame(read_excel("G2_ALL_311020c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x09_g2<- as.data.frame(read_excel("G2_ALL_091120c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

x16_g2<- as.data.frame(read_excel("G2_ALL_161220c.xlsx", col_names = TRUE, col_types = c("guess", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric", "numeric", "numeric", 
                                                                                         "numeric")))

g231 <- order.mus(x31_g2)
g209 <- order.mus(x09_g2)
g216 <- order.mus(x16_g2)
together_g2<-bind_rows(g231, g209, g216) %>% group_by(sample) %>% summarise_all(mean, na.rm=TRUE)
together_g2$sample<-factor(together_g2$sample, levels=sam) #make samples column, factor variable
g2<- together_g2 #placeholder dataframe for checking things if necessary later
g2out <- g2 %>% arrange(factor(sample, levels = sam)) #arrange the promoter column to start with tdh3
write.xlsx(g2out, "g2out.xlsx") #write as excel file






####JEFLABS flexmix fitting####
#https://jef.works/blog/2017/08/05/a-practical-introduction-to-finite-mixture-models/
#mo1-3 are based on gaussian distributions as we want to find a mix of gaussian distributions (subpopulations) 
#within the larger population of cells sampled during flow cytometry
mo1 <- FLXMRglm(family = "gaussian")
mo2 <- FLXMRglm(family = "gaussian")
mo3 <- FLXMRglm(family = "gaussian")

#function for making plots with area relative to population size
plot_mix_comps <- function(x, mu, sigma, lam) {lam * dnorm(x, mu, sigma)}



#Adding GFP values to plots
#Read in table containing metrics like Promoter, media, and plate reader GFP values
met<- as.data.frame(read_xlsx("C:/Users/davin/Documents/PhD/Results/FlowJo/Flow_scripts/means_sds_SUM/2023_means/long_map_labels.xlsx"))
##View(met)


#theme arguments for make.3/2/1dist functions
####tgg####
tgg <-theme(axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(family="Calibri Light", size=25,colour = "#919191"),
            axis.ticks.y = element_line(colour="#dddddd", size=1.5),
            axis.line.y = element_line(size = 2, colour = "#dddddd"),
            axis.ticks.length = unit(0.3, "cm"), 
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x = element_line(colour = "#dddddd", size=0.6),
            plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))



make.3dist <- function(m1,m2,m3,sd1,sd2,sd3,n1,n2,n3,pr,rf,index_no,bp){
  m1 <- m1
  m2 <- m2
  m3 <- m3
  sd1 <- sd1
  sd2 <- sd2 
  sd3 <- sd3
  p1 <- (n1/(n1+n2+n3))*2000
  p2 <- (n2/(n1+n2+n3))*2000
  p3 <- (n3/(n1+n2+n3))*2000
  
  set.seed(1)
  a <- rnorm(n=p1, mean=m1, sd=sd1)
  b <- rnorm(n=p2, mean=m2, sd=sd2)
  c <- rnorm(n=p3, mean=m3, sd=sd3)
  
  x <- c(a,b,c)
  class <- c(rep('a', p1), rep('b', p2),rep('c', p3))
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  ff<-flexmix(data$x~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  
  rawf <- data.frame(gfp = rf[[index_no]]) #to plot raw flow data
  ctrl <- data.frame(neg = rf[[bp]]) #to plot negative control (BP strain)
  
  #plot
  p <- ggplot(data, aes(x=x)) + 
    stat_function(fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#C96643", fill="#C96643",alpha=0.6, geom="polygon", size=2) + #RED
    stat_function(fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#84ADAD",fill="#84ADAD",alpha=0.6,geom="polygon", size=2) + #BLUE
    stat_function(fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#3B0535",fill="#3B0535",alpha=0.5,geom="polygon", size=2) + #purple
    ## geom_density(data=rawf, aes(x=gfp), color="#919191",size=2) + #plot raw data
    geom_density(data=ctrl, aes(x=neg), color="#dddddd",size=2,fill="#dddddd", alpha=0.3) + #plot negative control
    geom_density(data=rawf, aes(x=gfp), color=alpha("black",0.8),size=1) + #plot raw data
    geom_hline(yintercept = 0, colour = "#dddddd",size=2) +
    geom_vline(xintercept = pr, linetype = 1, colour = "#66c066", size=1.5) +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0,4.5), breaks = c(2,4)) +
    theme_minimal() +
    tgg +
    xlab("")
  return(p)
}


make.2dist <- function(m1,m2,sd1,sd2,n1,n2,pr,rf,index_no,bp){
  m1 <- m1
  m2 <- m2
  sd1 <- sd1
  sd2 <- sd2 
  p1 <- (n1/(n1+n2))*2000
  p2 <- (n2/(n1+n2))*2000
  
  set.seed(1)
  a <- rnorm(n=p1, mean=m1, sd=sd1)
  b <- rnorm(n=p2, mean=m2, sd=sd2)
  
  x <- c(a,b)
  class <- c(rep('a', p1), rep('b', p2))
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  ff<-flexmix(data$x~1, k=2, model=list(mo1, mo2))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  lam <- table(clusters(ff))
  
  rawf <- data.frame(gfp = rf[[index_no]])
  ctrl <- data.frame(neg = rf[[bp]])
  
  #plot
  p <- ggplot(data, aes(x=x)) + 
    stat_function(fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#84ADAD",fill="#84ADAD",alpha=0.6,geom="polygon", size=2) +
    stat_function(fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#C96643", fill="#C96643",alpha=0.6, geom="polygon",size=2) +
    ## geom_density(data=rawf, aes(x=gfp), color="#919191",size=2) + #plot raw data
    geom_density(data=ctrl, aes(x=neg), color="#dddddd",size=2,fill="#dddddd", alpha=0.3) + #plot negative control
    geom_density(data=rawf, aes(x=gfp), color=alpha("black",0.8),size=1) + #plot raw data
    geom_hline(yintercept = 0, colour = "#dddddd",size=2) +
    geom_vline(xintercept = pr, linetype = 1, colour = "#66c066", size=1.5) +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0,4.5), breaks = c(2,4)) +
    ##scale_y_continuous(limits = c(0,4), breaks = c(2,4) +
    theme_minimal() +
    tgg +
    xlab("")
  return(p)
}


make.1dist <- function(m1, sd1, n1, pr,rf,index_no,bp){
  
  a <- rnorm(n=n1, mean=m1, sd=sd1)
  x <- a
  class <- c(rep('a', n1))
  data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
  data$class<-as.factor(data$class)
  
  rawf <- data.frame(gfp = rf[[index_no]])
  ctrl <- data.frame(neg = rf[[bp]])
  
  p <- ggplot(data, aes(x=x)) + 
    stat_function(fun=plot_mix_comps,
                  args=list(mu=m1, 
                            sigma=sd1,
                            lam=n1/n1),
                  geom = "area",
                  size=2,
                  color="#D7BCAB",
                  fill="#D7BCAB", alpha=0.6) +
    ## geom_density(data=rawf, aes(x=gfp), color="#919191",size=2) + #plot raw data
    geom_density(data=ctrl, aes(x=neg), color="#dddddd",size=2,fill="#dddddd", alpha=0.3) + #plot negative control
    geom_density(data=rawf, aes(x=gfp), color=alpha("black",0.8),size=1) + #plot raw data
    geom_hline(yintercept = 0, colour = "#dddddd",size=2) +
    geom_vline(xintercept = pr, linetype = 1, colour = "#66c066", size=1.5) +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6),expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0,4.5), breaks = c(2,4)) +
    xlab("")+
    ylab("")+
    tgg
  return(p)
}



#testing raw plotting
#reading in raw data - each csv contains 50000 rows, for each single cell sampled, and each odd no. column is the data for 
#each strain expressing a single pYTK-sfGFP combination - 
#col1=cell count
#col3=pTDH3-sfGFP
#col5=pCCW12-sfGFP
#to remove the interfering cell count columns (even no. ones), the csvs will be indexed by odd number columns
#091120 dataset was used as the raw dataset for each plot - these files are named '091120_...' in https://github.com/davinap/YTK-heterogeneity-GFP

y1_raw <- read.csv("y1.csv")[, seq(1, ncol(read.csv("y1.csv")), by = 2)]
y2_raw <- read.csv("y2.csv")[, seq(1, ncol(read.csv("y2.csv")), by = 2)]
x1_raw <- read.csv("x1.csv")[, seq(1, ncol(read.csv("x1.csv")), by = 2)]
x2_raw <- read.csv("x2.csv")[, seq(1, ncol(read.csv("x2.csv")), by = 2)]
m1_raw <- read.csv("m1.csv")[, seq(1, ncol(read.csv("m1.csv")), by = 2)]
m2_raw <- read.csv("m2.csv")[, seq(1, ncol(read.csv("m2.csv")), by = 2)]
g1_raw <- read.csv("g1.csv")[, seq(1, ncol(read.csv("g1.csv")), by = 2)]
g2_raw <- read.csv("g2.csv")[, seq(1, ncol(read.csv("g2.csv")), by = 2)]
##y1_raw[[3]] #input for plotting raw data on top of mean data in make.1/2/3.dist functions

View(g1_raw)#looks good

#testing within plot


#test above######
#note, using y1_raw[[3]] as a test for plotting, eventhough the mean and raw data are not from the same sample
#this is why the two may not look aligned at the moment
sp1<- make.1dist(x1df[1,2], x1df[1,5], x1df[1,8], rs4[1,4], y1_raw,2,21) #tdh3=1pop, in y1_raw, the 3rd column contains tdh3 data
sp1 #YESSSSSS! ####
View(head(y1_raw))

#test 2 populations
sp4<- make.2dist(x1df[4,2], x1df[4,3],
                 x1df[4,5], x1df[4,6],
                 x1df[4,8], x1df[4,9],
                 rs4[4,4],
                 y1_raw,3,21) #HHF2 2pop
sp4 #good

#test 3 populations
sp2<- make.3dist(x1df[2,2], x1df[2,3], x1df[2,4],
                 x1df[2,5], x1df[2,6],x1df[2,7],
                 x1df[2,8], x1df[2,9], x1df[2,10],
                 rs4[2,4],
                 y1_raw,3,21) #ccw12 = 3pop
sp2 #good

View(m2_raw)

#240924 - start making all plots now that functions work

####Scaling method based on GFP stationary and GFP Exponential range####
##View(met)

#min and max: calculate the range of values in both GFP Exponential and Stationary columns in met
#scale_factor: calculate the value to multiply all GFP Exponential/Stationary values by and keep them in the 5.5 to 3.5 range
#5.5-3.5 were selected as they align best with the subpopulation plots
min_val<-min(met$`GFP Exponential`,met$`GFP Stationary`, na.rm = TRUE)
max_val<-max(met$`GFP Exponential`,met$`GFP Stationary`, na.rm = TRUE)
scale_factor <- (5.5-3.5)/(max_val-min_val)
scale_factor #4.211106e-05

#make plate column with rescaled values:
#For GFP Exponential values. They were multiplied by the scale_factor and added to 3.5 to keep them in the 3.5-5.5 range
rs4<- dplyr::select(met, Promoter, Media, `GFP Exponential`) %>% 
  mutate(plate=3.5+scale_factor*(`GFP Exponential`-min_val))
##View(rs4) #good

#same as in rs4 but for 'GFP Stationary' values in met dataframe
rs5<- dplyr::select(met, Promoter, Media, `GFP Stationary`) %>% 
  mutate(plate=3.5+scale_factor*(`GFP Stationary`-min_val))
##View(rs5) #good



####EXPONENTIAL-rs4####
#Each plot was made manually - ie. the table in Fig.1B was read and then a 'make.3/2/1dist' function was
#chosen in accordance with the table to plot the subpopulations.
#Each new object here is a ggplot object for a single promoter showing the consensus values of sfGFP expression
#across 3 biological replicates - or 2 biological replicates if one didn't agree with the others
#The comment at the end of each object refers to the consensus number of subpopulations as determined by 
#Hartigan's dipstest and the BIC value in the flexmix package - see supplementary

#ypd10 exponential
##View(x1df)
##View(x1out) #x1out and x1df are the same - looks like x1out used to be called x1df. Have kept 'x1df' below
x1df<-x1out
x1df<-as.data.frame(x1df) #otherwise can't call subsets

sp1<- make.1dist(x1df[1,2], x1df[1,5], x1df[1,8], rs4[1,4], x1_raw,2,21) #tdh3=1pop

sp2<- make.3dist(x1df[2,2], x1df[2,3], x1df[2,4],
                 x1df[2,5], x1df[2,6],x1df[2,7],
                 x1df[2,8], x1df[2,9], x1df[2,10],
                 rs4[2,4], x1_raw,3,21) #ccw12 = 3pop

sp3<- make.1dist(x1df[3,2], x1df[3,5], x1df[3,8], rs4[3,4],x1_raw,4,21) #pgk1 1pop

sp4<- make.2dist(x1df[4,2], x1df[4,3],
                 x1df[4,5], x1df[4,6],
                 x1df[4,8], x1df[4,9],
                 rs4[4,4],x1_raw,5,21) #HHF2 2pop

sp5<- make.1dist(x1df[5,2], x1df[5,5], x1df[5,8], rs4[5,4],x1_raw,6,21) #tef1 1pop

sp6<- make.1dist(x1df[6,2], x1df[6,5], x1df[6,8], rs4[6,4],x1_raw,7,21) #tef2 1pop

sp7<- make.2dist(x1df[7,2], x1df[7,3],
                 x1df[7,5], x1df[7,6],
                 x1df[7,8], x1df[7,9],
                 rs4[7,4],x1_raw,8,21)   #hhf1 2 pop

sp8<- make.2dist(x1df[8,2], x1df[8,3],
                 x1df[8,5], x1df[8,6],
                 x1df[8,8], x1df[8,9],
                 rs4[8,4],x1_raw,9,21) #htb2 2pop

sp9<- make.1dist(x1df[9,2], x1df[9,5], x1df[9,8], rs4[9,4],x1_raw,10,21) #rpl18b 1pop

sp10<- make.1dist(x1df[10,2], x1df[10,5], x1df[10,8], rs4[10,4],x1_raw,11,21) #ald6 1pop

sp11<- make.1dist(x1df[11,2], x1df[11,5], x1df[11,8], rs4[11,4],x1_raw,12,21) #apab1 1pop

sp12<- make.1dist(x1df[12,2], x1df[12,5], x1df[12,8], rs4[12,4],x1_raw,13,21) #ret2 1pop

sp13<- make.2dist(x1df[13,2], x1df[13,3],
                  x1df[13,5], x1df[13,6],
                  x1df[13,8], x1df[13,9],
                  rs4[13,4], x1_raw,14,21) #rnr1 2pop

sp14<- make.1dist(x1df[14,2], x1df[14,5], x1df[14,8], rs4[14,4],x1_raw,15,21) #sac6 1pop

sp15<- make.1dist(x1df[15,2], x1df[15,5], x1df[15,8], rs4[15,4],x1_raw,16,21) #rnr2 1pop

sp16<- make.1dist(x1df[16,2], x1df[16,5], x1df[16,8], rs4[16,4],x1_raw,17,21) #pop6 1pop

sp17<- make.3dist(x1df[17,2], x1df[17,3], x1df[17,4],
                  x1df[17,5], x1df[17,6],x1df[17,7],
                  x1df[17,8], x1df[17,9], x1df[17,10],
                  rs4[17,4],x1_raw,18,21) #rad27 = 3pop

sp18<- make.1dist(x1df[18,2], x1df[18,5], x1df[18,8], rs4[18,4],x1_raw,19,21) #psp2 1pop

sp19<- make.1dist(x1df[19,2], x1df[19,5], x1df[19,8], rs4[19,4],x1_raw,20,21) #rev1 1pop



#plot all ypd10 exponential phase distributions together
x1plots<- ggarrange(sp1,sp2,sp3,sp4,sp5,sp6,sp7,sp8,sp9,sp10,sp11,sp12,sp13,sp14,sp15,sp16,sp17,sp18,sp19,
                    ncol = 1, nrow = 19)

##x1plots




####ypd exponential
##View(y1out)
##class(y1out)
##class(x1df)
y1out <- as.data.frame(y1out)

##y1out[1,3]

yp1<- make.1dist(y1out[1,2], y1out[1,5], y1out[1,8], rs4[21,4], y1_raw,2,21) #tdh3=1pop

yp2<- make.1dist(y1out[2,2], y1out[2,5], y1out[2,8], rs4[22,4],y1_raw,3,21) #ccw12 = 1pop

yp3<- make.1dist(y1out[3,2], y1out[3,5], y1out[3,8], rs4[23,4],y1_raw,4,21) #pgk1 1pop

yp4<- make.3dist(y1out[4,2], y1out[4,3], y1out[4,4],
                 y1out[4,5], y1out[4,6],y1out[4,7],
                 y1out[4,8], y1out[4,9], y1out[4,10],
                 rs4[24,4],y1_raw,5,21) #HHF2 3pop

yp5<- make.1dist(y1out[5,2], y1out[5,5], y1out[5,8], rs4[25,4],y1_raw,6,21) #tef1 1pop

yp6<- make.1dist(y1out[6,2], y1out[6,5], y1out[6,8], rs4[26,4],y1_raw,7,21) #tef2 1pop

yp7<- make.3dist(y1out[7,2], y1out[7,3], y1out[7,4],
                 y1out[7,5], y1out[7,6],y1out[7,7],
                 y1out[7,8], y1out[7,9], y1out[7,10],
                 rs4[27,4],y1_raw,8,21)  #hhf1 3 pop

yp8<- make.3dist(y1out[8,2], y1out[8,3], y1out[8,4],
                 y1out[8,5], y1out[8,6],y1out[8,7],
                 y1out[8,8], y1out[8,9], y1out[8,10],
                 rs4[28,4],y1_raw,9,21) #htb2 3pop

yp9<-  make.1dist(y1out[9,2], y1out[9,5], y1out[9,8], rs4[29,4],y1_raw,10,21) #rpl18b 1pop

yp10<- make.1dist(y1out[10,2], y1out[10,5], y1out[10,8], rs4[30,4],y1_raw,11,21) #ald6 1pop

yp11<- make.1dist(y1out[11,2], y1out[11,5], y1out[11,8], rs4[31,4],y1_raw,12,21) #apab1 1pop

yp12<- make.1dist(y1out[12,2], y1out[12,5], y1out[12,8], rs4[32,4],y1_raw,13,21) #ret2 1pop

yp13<- make.3dist(y1out[13,2], y1out[13,3], y1out[13,4],
                  y1out[13,5], y1out[13,6],y1out[13,7],
                  y1out[13,8], y1out[13,9], y1out[13,10],
                  rs4[33,4],y1_raw,14,21) #rnr1 3pop

yp14<- make.1dist(y1out[14,2], y1out[14,5], y1out[14,8], rs4[34,4],y1_raw,15,21) #sac6 1pop

yp15<- make.1dist(y1out[15,2], y1out[15,5], y1out[15,8], rs4[35,4],y1_raw,16,21) #rnr2 1pop

yp16<- make.1dist(y1out[16,2], y1out[16,5], y1out[16,8], rs4[36,4],y1_raw,17,21) #pop6 1pop

yp17<- make.1dist(y1out[17,2], y1out[17,5], y1out[17,8], rs4[37,4],y1_raw,18,21) #rad27 = 1pop

yp18<- make.1dist(y1out[18,2], y1out[18,5], y1out[18,8], rs4[38,4],y1_raw,19,21) #pyp2 1pop

yp19<- make.1dist(y1out[19,2], y1out[19,5], y1out[19,8], rs4[39,4],y1_raw,20,21) #rev1 1pop




y1plots<- ggarrange(yp1,yp2,yp3,yp4,yp5,yp6,yp7,yp8,yp9,yp10,yp11,yp12,yp13,yp14,yp15,yp16,yp17,yp18,yp19,
                    ncol = 1, nrow = 19)

##y1plots


####ynb exponential
m1out <- as.data.frame(m1out)


mp1<- make.1dist(m1out[1,2], m1out[1,5], m1out[1,8], rs4[41,4],m1_raw,2,21) #tdh3=1pop

mp2<- make.3dist(m1out[2,2], m1out[2,3], m1out[2,4],
                 m1out[2,5], m1out[2,6],m1out[2,7],
                 m1out[2,8], m1out[2,9], m1out[2,10],
                 rs4[42,4],m1_raw,3,21) #ccw12 = 3pop

mp3<- make.1dist(m1out[3,2], m1out[3,5], m1out[3,8], rs4[43,4],m1_raw,4,21) #pgk1 1pop

mp4<- make.3dist(m1out[4,2], m1out[4,3], m1out[4,4],
                 m1out[4,5], m1out[4,6],m1out[4,7],
                 m1out[4,8], m1out[4,9], m1out[4,10],
                 rs4[44,4],m1_raw,5,21) #HHF2 3pop

mp5<- make.1dist(m1out[5,2], m1out[5,5], m1out[5,8], rs4[45,4],m1_raw,6,21) #tef1 1pop

mp6<- make.1dist(m1out[6,2], m1out[6,5], m1out[6,8], rs4[46,4],m1_raw,7,21) #tef2 1pop

mp7<- make.3dist(m1out[7,2], m1out[7,3], m1out[7,4],
                 m1out[7,5], m1out[7,6],m1out[7,7],
                 m1out[7,8], m1out[7,9], m1out[7,10],
                 rs4[47,4],m1_raw,8,21)  #hhf1 3 pop

mp8<- make.3dist(m1out[8,2], m1out[8,3], m1out[8,4],
                 m1out[8,5], m1out[8,6],m1out[8,7],
                 m1out[8,8], m1out[8,9], m1out[8,10],
                 rs4[48,4],m1_raw,9,21) #htb2 3pop

mp9<-  make.1dist(m1out[9,2], m1out[9,5], m1out[9,8], rs4[49,4],m1_raw,10,21) #rpl18b 1pop

mp10<- make.1dist(m1out[10,2], m1out[10,5], m1out[10,8], rs4[50,4],m1_raw,11,21) #ald6 1pop

mp11<- make.1dist(m1out[11,2], m1out[11,5], m1out[11,8], rs4[51,4],m1_raw,12,21) #apab1 1pop

mp12<- make.1dist(m1out[12,2], m1out[12,5], m1out[12,8], rs4[52,4],m1_raw,13,21) #ret2 1pop

mp13<- make.3dist(m1out[13,2], m1out[13,3], m1out[13,4],
                  m1out[13,5], m1out[13,6],m1out[13,7],
                  m1out[13,8], m1out[13,9], m1out[13,10],
                  rs4[53,4],m1_raw,14,21) #rnr1 3pop

mp14<- make.1dist(m1out[14,2], m1out[14,5], m1out[14,8], rs4[54,4],m1_raw,15,21) #sac6 1pop

mp15<- make.1dist(m1out[15,2], m1out[15,5], m1out[15,8], rs4[55,4],m1_raw,16,21) #rnr2 1pop

mp16<- make.1dist(m1out[16,2], m1out[16,5], m1out[16,8], rs4[56,4],m1_raw,17,21) #pop6 1pop

mp17<- make.1dist(m1out[17,2], m1out[17,5], m1out[17,8], rs4[57,4],m1_raw,18,21) #rad27 = 1pop

mp18<- make.1dist(m1out[18,2], m1out[18,5], m1out[18,8], rs4[58,4],m1_raw,19,21) #pmp2 1pop

mp19<- make.1dist(m1out[19,2], m1out[19,5], m1out[19,8], rs4[59,4],m1_raw,20,21) #rev1 1pop



m1plots<- ggarrange(mp1,mp2,mp3,mp4,mp5,mp6,mp7,mp8,mp9,mp10,mp11,mp12,mp13,mp14,mp15,mp16,mp17,mp18,mp19,
                    ncol = 1, nrow = 19)



####yepg exponential
g1out <- as.data.frame(g1out)

#gp1 =tdh3 =no growth - the number 1 was used instead as mu1, sd1, lam1 and pr arguments
#for the sole purpose of maintaining plot geometry in ggarrange
#this plot was removed in the final supplementary figure
gp1<- make.1dist(1, 1, 1, 1,1,1,1) #tdh3=no growth, 1 pop placeholder

gp2<- make.3dist(g1out[2,2], g1out[2,3], g1out[2,4],
                 g1out[2,5], g1out[2,6],g1out[2,7],
                 g1out[2,8], g1out[2,9], g1out[2,10],
                 rs4[62,4],g1_raw,3,21) #ccw12 = 3pop

gp3<- make.1dist(g1out[3,2], g1out[3,5], g1out[3,8], rs4[63,4],g1_raw,4,21) #pgk1 1pop

gp4<- make.2dist(g1out[4,2], g1out[4,3],
                 g1out[4,5], g1out[4,6],
                 g1out[4,8], g1out[4,9],
                 rs4[64,4],g1_raw,5,21) #HHF2 2pop

gp5<- make.1dist(g1out[5,2], g1out[5,5], g1out[5,8], rs4[65,4],g1_raw,6,21) #tef1 1pop

gp6<- make.1dist(g1out[6,2], g1out[6,5], g1out[6,8], rs4[66,4],g1_raw,7,21) #tef2 1pop

gp7<- make.2dist(g1out[7,2], g1out[7,3],
                 g1out[7,5], g1out[7,6],
                 g1out[7,8], g1out[7,9],
                 rs4[67,4],g1_raw,8,21)  #hhf1 2 pop

gp8<- make.2dist(g1out[8,2], g1out[8,3],
                 g1out[8,5], g1out[8,6],
                 g1out[8,8], g1out[8,9],
                 rs4[68,4],g1_raw,9,21) #htb2 2pop

gp9<-  make.1dist(g1out[9,2], g1out[9,5], g1out[9,8], rs4[69,4],g1_raw,10,21) #rpl18b 1pop

gp10<- make.1dist(g1out[10,2], g1out[10,5], g1out[10,8], rs4[70,4],g1_raw,11,21) #ald6 1pop

gp11<- make.1dist(g1out[11,2], g1out[11,5], g1out[11,8], rs4[71,4],g1_raw,12,21) #pab1 1pop

gp12<- make.1dist(g1out[12,2], g1out[12,5], g1out[12,8], rs4[72,4],g1_raw,13,21) #ret2 1pop

gp13<- make.3dist(g1out[13,2], g1out[13,3], g1out[13,4],
                  g1out[13,5], g1out[13,6],g1out[13,7],
                  g1out[13,8], g1out[13,9], g1out[13,10],
                  rs4[73,4],g1_raw,14,21) #rnr1 3pop

gp14<- make.1dist(g1out[14,2], g1out[14,5], g1out[14,8], rs4[74,4],g1_raw,15,21) #sac6 1pop

gp15<- make.1dist(g1out[15,2], g1out[15,5], g1out[15,8], rs4[75,4],g1_raw,16,21) #rnr2 1pop

gp16<- make.1dist(g1out[16,2], g1out[16,5], g1out[16,8], rs4[76,4],g1_raw,17,21) #pop6 1pop

gp17<- make.1dist(g1out[17,2], g1out[17,5], g1out[17,8], rs4[77,4],g1_raw,18,21) #rad27 = 1pop

gp18<- make.1dist(g1out[18,2], g1out[18,5], g1out[18,8], rs4[78,4],g1_raw,19,21) #psp2 1pop

gp19<- make.1dist(g1out[19,2], g1out[19,5], g1out[19,8], rs4[79,4],g1_raw,20,21) #rev1 1pop



g1plots<- ggarrange(gp1,gp2,gp3,gp4,gp5,gp6,gp7,gp8,gp9,gp10,gp11,gp12,gp13,gp14,gp15,gp16,gp17,gp18,gp19,
                    ncol = 1, nrow = 19)





####STATIONARY - rs5####
#ypd10 stationary
x2out<-as.data.frame(x2out)

a1<- make.2dist(x2out[1,2], x2out[1,3],
                x2out[1,5], x2out[1,6],
                x2out[1,8], x2out[1,9],
                rs5[1,4],x2_raw,2,21) #tdh3=2pop

a2<- make.2dist(x2out[2,2], x2out[2,3],
                x2out[2,5], x2out[2,6],
                x2out[2,8], x2out[2,9],
                rs5[2,4],x2_raw,3,21) #ccw12 = 2pop

a3<- make.1dist(x2out[3,2], x2out[3,5], x2out[3,8], rs5[3,4],x2_raw,4,21) #pgk1 1pop

a4<- make.2dist(x2out[4,2], x2out[4,3],
                x2out[4,5], x2out[4,6],
                x2out[4,8], x2out[4,9],
                rs5[4,4],x2_raw,5,21) #HHF2 2pop

a5<- make.1dist(x2out[5,2], x2out[5,5], x2out[5,8], rs5[5,4],x2_raw,6,21) #tef1 1pop

a6<- make.1dist(x2out[6,2], x2out[6,5], x2out[6,8], rs5[6,4],x2_raw,7,21) #tef2 1pop

a7<- make.2dist(x2out[7,2], x2out[7,3],
                x2out[7,5], x2out[7,6],
                x2out[7,8], x2out[7,9],
                rs5[7,4],x2_raw,8,21)   #hhf1 2 pop

a8<- make.3dist(x2out[8,2], x2out[8,3], x2out[8,4],
                x2out[8,5], x2out[8,6],x2out[8,7],
                x2out[8,8], x2out[8,9], x2out[8,10],
                rs5[8,4],x2_raw,9,21) #htb2 3pop

a9<- make.3dist(x2out[9,2], x2out[9,3], x2out[9,4],
                x2out[9,5], x2out[9,6],x2out[9,7],
                x2out[9,8], x2out[9,9], x2out[9,10],
                rs5[9,4],x2_raw,10,21) #rpl18b 3pop

a10<- make.1dist(x2out[10,2], x2out[10,5], x2out[10,8], rs5[10,4],x2_raw,11,21) #ald6 1pop

a11<- make.1dist(x2out[11,2], x2out[11,5], x2out[11,8], rs5[11,4],x2_raw,12,21) #pab1 1pop

a12<- make.2dist(x2out[12,2], x2out[12,3],
                 x2out[12,5], x2out[12,6],
                 x2out[12,8], x2out[12,9],
                 rs5[12,4],x2_raw,13,21) #ret2 2pop

a13<- make.3dist(x2out[13,2], x2out[13,3], x2out[13,4],
                 x2out[13,5], x2out[13,6],x2out[13,7],
                 x2out[13,8], x2out[13,9], x2out[13,10],
                 rs5[13,4],x2_raw,14,21) #rnr1 3pop

a14<- make.1dist(x2out[14,2], x2out[14,5], x2out[14,8], rs5[14,4],x2_raw,15,21) #sac6 1pop

a15<- make.1dist(x2out[15,2], x2out[15,5], x2out[15,8], rs5[15,4],x2_raw,16,21) #rnr2 1pop

a16<- make.1dist(x2out[16,2], x2out[16,5], x2out[16,8], rs5[16,4],x2_raw,17,21) #pop6 1pop

a17<- make.2dist(x2out[17,2], x2out[17,3],
                 x2out[17,5], x2out[17,6],
                 x2out[17,8], x2out[17,9],
                 rs5[17,4],x2_raw,18,21) #rad27 = 2pop

a18<- make.1dist(x2out[18,2], x2out[18,5], x2out[18,8], rs5[18,4],x2_raw,19,21) #psp2 1pop

a19<- make.1dist(x2out[19,2], x2out[19,5], x2out[19,8], rs5[19,4],x2_raw,20,21) #rev1 1pop


x2plots_rs5<- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,
                        ncol = 1, nrow = 19)


####ypd stationary
y2out<-as.data.frame(y2out)

b1<- make.1dist(y2out[1,2], y2out[1,5], y2out[1,8], rs5[21,4],y2_raw,2,21) #tdh3=1pop

b2<- make.3dist(y2out[2,2], y2out[2,3], y2out[2,4],
                y2out[2,5], y2out[2,6],y2out[2,7],
                y2out[2,8], y2out[2,9], y2out[2,10],
                rs5[22,4],y2_raw,3,21) #ccw12 = 3pop

b3<- make.1dist(y2out[3,2], y2out[3,5], y2out[3,8], rs5[23,4],y2_raw,4,21) #pgk1 1pop

b4<- make.3dist(y2out[4,2], y2out[4,3], y2out[4,4],
                y2out[4,5], y2out[4,6],y2out[4,7],
                y2out[4,8], y2out[4,9], y2out[4,10],
                rs5[24,4],y2_raw,5,21) #HHF2 3pop

b5<- make.1dist(y2out[5,2], y2out[5,5], y2out[5,8], rs5[25,4],y2_raw,6,21) #tef1 1pop
  
b6<- make.1dist(y2out[6,2], y2out[6,5], y2out[6,8], rs5[26,4],y2_raw,7,21) #tef2 1pop

b7<- make.3dist(y2out[7,2], y2out[7,3], y2out[7,4],
                y2out[7,5], y2out[7,6],y2out[7,7],
                y2out[7,8], y2out[7,9], y2out[7,10],
                rs5[27,4],y2_raw,8,21)   #hhf1 3 pop

b8<- make.3dist(y2out[8,2], y2out[8,3], y2out[8,4],
                y2out[8,5], y2out[8,6],y2out[8,7],
                y2out[8,8], y2out[8,9], y2out[8,10],
                rs5[28,4],y2_raw,9,21) #htb2 3pop

b9<-  make.1dist(y2out[9,2], y2out[9,5], y2out[9,8], rs5[29,4],y2_raw,10,21)#rpl18b 1pop

b10<- make.1dist(y2out[10,2], y2out[10,5], y2out[10,8], rs5[30,4],y2_raw,11,21) #ald6 1pop

b11<- make.1dist(y2out[11,2], y2out[11,5], y2out[11,8], rs5[31,4],y2_raw,12,21) #pab1 1pop

b12<- make.1dist(y2out[12,2], y2out[12,5], y2out[12,8], rs5[32,4],y2_raw,13,21) #ret2 1pop

b13<- make.3dist(y2out[13,2], y2out[13,3], y2out[13,4],
                 y2out[13,5], y2out[13,6],y2out[13,7],
                 y2out[13,8], y2out[13,9], y2out[13,10],
                 rs5[33,4],y2_raw,14,21) #rnr1 3pop

b14<- make.1dist(y2out[14,2], y2out[14,5], y2out[14,8], rs5[34,4],y2_raw,15,21) #sac6 1pop

b15<- make.1dist(y2out[15,2], y2out[15,5], y2out[15,8], rs5[35,4],y2_raw,16,21) #rnr2 1pop

b16<- make.1dist(y2out[16,2], y2out[16,5], y2out[16,8], rs5[36,4],y2_raw,17,21) #pop6 1pop

b17<- make.2dist(y2out[17,2], y2out[17,3],
                 y2out[17,5], y2out[17,6],
                 y2out[17,8], y2out[17,9],
                 rs5[37,4],y2_raw,18,21) #rad27 = 2pop

b18<- make.1dist(y2out[18,2], y2out[18,5], y2out[18,8], rs5[38,4],y2_raw,19,21) #psp2 1pop

b19<- make.1dist(y2out[19,2], y2out[19,5], y2out[19,8], rs5[39,4],y2_raw,20,21) #rev1 1pop


y2plots_rs5<- ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,
                        ncol = 1, nrow = 19)


####ynb stationary
m2out<-as.data.frame(m2out)

c1<- make.1dist(m2out[1,2], m2out[1,5], m2out[1,8], rs5[41,4],m2_raw,2,21) #tdh3=1pop

c2<- make.3dist(m2out[2,2], m2out[2,3], m2out[2,4],
                m2out[2,5], m2out[2,6],m2out[2,7],
                m2out[2,8], m2out[2,9], m2out[2,10],
                rs5[42,4], m2_raw,3,21) #ccw12 = 3pop

c3<- make.1dist(m2out[3,2], m2out[3,5], m2out[3,8], rs5[43,4],m2_raw,4,21) #pgk1 1pop

c4<- make.3dist(m2out[4,2], m2out[4,3], m2out[4,4],
                m2out[4,5], m2out[4,6],m2out[4,7],
                m2out[4,8], m2out[4,9], m2out[4,10],
                rs5[44,4],m2_raw,5,21) #HHF2 3pop

c5<- make.1dist(m2out[5,2], m2out[5,5], m2out[5,8], rs5[45,4],m2_raw,6,21) #tef1 1pop

c6<- make.1dist(m2out[6,2], m2out[6,5], m2out[6,8], rs5[46,4],m2_raw,7,21) #tef2 1pop

c7<- make.3dist(m2out[7,2], m2out[7,3], m2out[7,4],
                m2out[7,5], m2out[7,6],m2out[7,7],
                m2out[7,8], m2out[7,9], m2out[7,10],
                rs5[47,4],m2_raw,8,21)   #hhf1 3 pop

c8<- make.3dist(m2out[8,2], m2out[8,3], m2out[8,4],
                m2out[8,5], m2out[8,6], m2out[8,7],
                m2out[8,8], m2out[8,9], m2out[8,10],
                rs5[48,4],m2_raw,9,21) #htb2 3pop

c9<-  make.1dist(m2out[9,2], m2out[9,5], m2out[9,8], rs5[49,4],m2_raw,10,21)#rpl18b 1pop

c10<- make.1dist(m2out[10,2], m2out[10,5], m2out[10,8], rs5[50,4],m2_raw,11,21) #ald6 1pop

c11<- make.1dist(m2out[11,2], m2out[11,5], m2out[11,8], rs5[51,4],m2_raw,12,21) #pab1 1pop

c12<- make.1dist(m2out[12,2], m2out[12,5], m2out[12,8], rs5[52,4],m2_raw,13,21) #ret2 1pop

c13<- make.3dist(m2out[13,2], m2out[13,3], m2out[13,4],
                 m2out[13,5], m2out[13,6],m2out[13,7],
                 m2out[13,8], m2out[13,9], m2out[13,10],
                 rs5[53,4],m2_raw,14,21) #rnr1 3pop

c14<- make.1dist(m2out[14,2], m2out[14,5], m2out[14,8], rs5[54,4],m2_raw,15,21) #sac6 1pop

c15<- make.1dist(m2out[15,2], m2out[15,5], m2out[15,8], rs5[55,4],m2_raw,16,21) #rnr2 1pop

c16<- make.1dist(m2out[16,2], m2out[16,5], m2out[16,8], rs5[56,4],m2_raw,17,21) #pop6 1pop

c17<- make.3dist(m2out[17,2], m2out[17,3], m2out[17,4],
                 m2out[17,5], m2out[17,6],m2out[17,7],
                 m2out[17,8], m2out[17,9], m2out[17,10],
                 rs5[57,4],m2_raw,18,21) #rad27 = 3pop

c18<- make.1dist(m2out[18,2], m2out[18,5], m2out[18,8], rs5[58,4],m2_raw,19,21) #psp2 1pop

c19<- make.1dist(m2out[19,2], m2out[19,5], m2out[19,8], rs5[59,4],m2_raw,20,21) #rev1 1pop


m2plots_rs5<- ggarrange(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,
                        ncol = 1, nrow = 19)


####yepg stationary
g2out<-as.data.frame(g2out)

d1<- make.1dist(1, 1, 1, 1,1,1,1) #tdh3=1pop #no growth in this condition, so placeholder value '1' put instead

d2<- make.2dist(g2out[2,2], g2out[2,3],
                g2out[2,5], g2out[2,6],
                g2out[2,8], g2out[2,9],
                rs5[62,4],g2_raw,3,21) #ccw12 = 2pop

d3<- make.1dist(g2out[3,2], g2out[3,5], g2out[3,8], rs5[63,4],g2_raw,4,21) #pgk1 1pop

d4<- make.2dist(g2out[4,2], g2out[4,3],
                g2out[4,5], g2out[4,6],
                g2out[4,8], g2out[4,9],
                rs5[64,4],g2_raw,5,21) #HHF2 2pop

d5<- make.1dist(g2out[5,2], g2out[5,5], g2out[5,8], rs5[65,4],g2_raw,6,21) #tef1 1pop

d6<- make.1dist(g2out[6,2], g2out[6,5], g2out[6,8], rs5[66,4],g2_raw,7,21) #tef2 1pop

d7<- make.2dist(g2out[7,2], g2out[7,3],
                g2out[7,5], g2out[7,6],
                g2out[7,8], g2out[7,9],
                rs5[67,4],g2_raw,8,21)   #hhf1 2 pop

d8<- make.2dist(g2out[8,2], g2out[8,3],
                g2out[8,5], g2out[8,6],
                g2out[8,8], g2out[8,9],
                rs5[68,4],g2_raw,9,21) #htb2 2pop

d9<-  make.1dist(g2out[9,2], g2out[9,5], g2out[9,8], rs5[69,4],g2_raw,10,21)#rpl18b 1pop

d10<- make.1dist(g2out[10,2], g2out[10,5], g2out[10,8], rs5[70,4],g2_raw,11,21) #ald6 1pop

d11<- make.1dist(g2out[11,2], g2out[11,5], g2out[11,8], rs5[71,4],g2_raw,12,21) #pab1 1pop

d12<- make.1dist(g2out[12,2], g2out[12,5], g2out[12,8], rs5[72,4],g2_raw,13,21) #ret2 1pop

d13<- make.3dist(g2out[13,2], g2out[13,3], g2out[13,4],
                 g2out[13,5], g2out[13,6],g2out[13,7],
                 g2out[13,8], g2out[13,9], g2out[13,10],
                 rs5[73,4],g2_raw,14,21) #rnr1 3pop

d14<- make.1dist(g2out[14,2], g2out[14,5], g2out[14,8], rs5[74,4],g2_raw,15,21) #sac6 1pop

d15<- make.1dist(g2out[15,2], g2out[15,5], g2out[15,8], rs5[75,4],g2_raw,16,21) #rnr2 1pop

d16<- make.1dist(g2out[16,2], g2out[16,5], g2out[16,8], rs5[76,4],g2_raw,17,21) #pop6 1pop

d17<- make.1dist(g2out[17,2], g2out[17,5], g2out[17,8], rs5[77,4],g2_raw,18,21) #rad27 1pop

d18<- make.1dist(g2out[18,2], g2out[18,5], g2out[18,8], rs5[78,4],g2_raw,19,21) #psp2 1pop

d19<- make.1dist(g2out[19,2], g2out[19,5], g2out[19,8], rs5[79,4],g2_raw,20,21) #rev1 1pop


g2plots_rs5<- ggarrange(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,
                        ncol = 1, nrow = 19)





####plotting with rs4 and rs5####
#put together all the exponential phase plots:
expplots_rs4<- ggarrange(x1plots,y1plots,m1plots,g1plots,
                         ncol=4, nrow=1)

#put together all the stationary phase plots:
statplots_rs5<- ggarrange(x2plots_rs5, y2plots_rs5,m2plots_rs5,g2plots_rs5,
                          ncol = 4, nrow=1)

expplots_rs4 #for supplementary
statplots_rs5 #for supplementary

x1plots #for Figure 1






