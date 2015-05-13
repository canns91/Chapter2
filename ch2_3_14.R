#wd
setwd("C:/Users/Cara/Documents/CBL/GitHub/Chapter2")

###packages
library(OpenStreetMap)
library(ggplot2)
library(dplyr)
library(tidyr)

#functions

binar<- function(x)
{
  ifelse(is.na(x),0,1)
}

prop<- function(x)
{
  sum(x)/length(x)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#files
larv<-read.csv("larv_6_4_14.csv", header=T)
loz<-read.csv("lozano.csv", header=T)

#fixing up the larval file

#remove regions of no occurace
#identify regions of low occurance 
#regs<-aggregate(ABUNDANCE~REGION,larv,sum)

#larva with only regions 1 and 2, MAB and SNE
larv<-larv[larv$REGION<3,]

#group months in pairs to match survey design
larv$MEAN_MOs<-factor(ifelse(larv$MEAN_MO==1 | larv$MEAN_MO==2, "1-2",
                             ifelse(larv$MEAN_MO==3|larv$MEAN_MO==4, "3-4",
                                    ifelse(larv$MEAN_MO==5|larv$MEAN_MO==6, "5-6",
                                           ifelse(larv$MEAN_MO==7|larv$MEAN_MO==8, "7-8",
                                                  ifelse(larv$MEAN_MO==9|larv$MEAN_MO==10, "9-10","11-12"))))))

###found an error!!! The 37 reported in the above cruise is most likely meant to be a 3.7 
#replace in larv_rms
larv$LENGTH[which(larv$LENGTH==37)]<-3.7


#the largest should be removed because broadening begins to 
#occur around 30 mm and I can't reliably predict age for the 36 mm
larv$LENGTH[which(larv$UNIQUE=="EK8001_31")]<-NA
larv$LENGTH_ABUNDANCE[which(larv$UNIQUE=="EK8001_31")]<-NA
larv$TOTAL_COUNT[which(larv$UNIQUE=="EK8001_31")]<-0
larv$ABUNDANCE[which(larv$UNIQUE=="EK8001_31")]<-0
larv$COUNT_AT_LENGTH[which(larv$UNIQUE=="EK8001_31")]<-0
SEASON<-NULL
for (x in 1:nrow(larv))
{
  h<-ifelse(as.numeric(as.character(larv$MONTH[x]))>7, as.numeric(as.character(larv$YEAR[x]))+1
            ,as.numeric(as.character(larv$YEAR[x])))
  SEASON<-c(SEASON,h)
}
larv$SEASON<-as.factor(SEASON)

#extract age/length data from carlos's catch
lenage<-loz[names(loz) %in% c("age","bclength")]
lenage<-aggregate(bclength~age,lenage,mean)
lenage<-rbind(c(0,0),lenage)
lenage<-lenage[c(1:69,76),]
lenage<-lenage[c(1,3:70),]


#add length info 
larv_no0<-larv[larv$ABUNDANCE!=0,]
#larv_no0$UNIQUE[which(is.na(larv_no0$LENGTH))]
#there are 5 tows where no menhaden are measured, although some are observed
#****I could estimate the length from similar timed cruises
ages_c<-NULL
test<-larv_no0$LENGTH
for (i in test)
{
  h<-min(which(lenage$bclength>=i))
  x_1<-lenage[h-1,]
  x_2<-lenage[h,]
  m<-(x_2[1]-x_1[1])/(x_2[2]-x_1[2])
  y<-m*(i-x_1[2])+x_1[1]
  ages_c<-c(ages_c,y)
}
ages_c<-unlist(ages_c)
larv_no0$AGE<-ages_c

larv_w0<-merge(larv,larv_no0, all=T)

larv_w0$MEAN_YR<-as.factor(larv_w0$MEAN_YR)
larv_w0$MEAN_MO<-as.factor(larv_w0$MEAN_MO)
larv_w0$PLANKTON_STRATUM<-as.factor(larv_w0$PLANKTON_STRATUM)

larv<-larv_w0

#to start, lets remove all individuals with estimated ages above a week
  #I need to be careful not to remove the tows 

larve<-larv[names(larv) %in% c("CRUISE_NAME","LATITUDE","LONGITUDE","PLANKTON_STRATUM","EVENT_JUL","TOTAL_COUNT","ABUNDANCE","COUNT_AT_LENGTH","LENGTH","LENGTH_ABUNDANCE","UNIQUE","MEAN_MOs","SEASON","AGE")]
larve_u<-unique(larv[names(larv) %in% c("CRUISE_NAME","PLANKTON_STRATUM","LATITUDE","LONGITUDE","EVENT_JUL","TOTAL_COUNT","ABUNDANCE","MEAN_MO","MEAN_MOs","SEASON","UNIQUE")])

larve_u$MEAN_MOs<-factor(larve_u$MEAN_MOs,levels(larve_u$MEAN_MOs)[c(5,6,2,1,3,4)])

larve_u<-larve_u[order(larve_u[,11],larve_u[,10]),]


larva<-larve[larve$LENGTH<=5.3,]
larva_n0<-larva[!is.na(larva$AGE),]

larva<-merge(larve_u,larva_n0,all=T)

abuna<-NULL
for (x in 1:nrow(larve_u))
{
  h<-larva[larva$UNIQUE==larve_u$UNIQUE[x],]
  i<-c(sum(h$COUNT_AT_LENGTH),sum(h$LENGTH_ABUNDANCE))
  abuna<-rbind(abuna,i)
}
##adding the adjusted abundance and counts back into the matrix
larve_u$ADJ_ABUNa<-abuna[,2]
larve_u$ADJ_COUNTa<-abuna[,1]
larve_u$ABUNbina <- factor(ifelse(is.na(larve_u$ADJ_ABUNa) , "0",
                                  ifelse(larve_u$ADJ_ABUNa > 0 & larve_u$ADJ_ABUNa <= 10,"0-10",
                                         ifelse(larve_u$ADJ_ABUNa > 10 & larve_u$ADJ_ABUNa <= 100, '10-100',
                                                '100+'))), 
                           levels=c('0','0-10','10-100','100+'),
                           ordered=T)


#second group
larvb<-larve[larve$LENGTH>5.3 & larve$LENGTH<=8.8,]
larvb_n0<-larvb[!is.na(larvb$AGE),]

larvb<-merge(larve_u,larvb_n0,all=T)

abunb<-NULL
for (x in 1:nrow(larve_u))
{
  h<-larvb[larvb$UNIQUE==larve_u$UNIQUE[x],]
  i<-c(sum(h$COUNT_AT_LENGTH),sum(h$LENGTH_ABUNDANCE))
  abunb<-rbind(abunb,i)
}

larve_u$ADJ_ABUNb<-abunb[,2]
larve_u$ADJ_COUNTb<-abunb[,1]
larve_u$ABUNbinb <- factor(ifelse(is.na(larve_u$ADJ_ABUNb) , "0",
                                  ifelse(larve_u$ADJ_ABUNb > 0 & larve_u$ADJ_ABUNb <= 10,"0-10",
                                         ifelse(larve_u$ADJ_ABUNb > 10 & larve_u$ADJ_ABUNb <= 100, '10-100',
                                                '100+'))), 
                           levels=c('0','0-10','10-100','100+'),
                           ordered=T)
#third group
larvc<-larve[larve$LENGTH>8.8 & larve$LENGTH<=13,]
larvc_n0<-larvc[!is.na(larvc$AGE),]

larvc<-merge(larve_u,larvc_n0,all=T)

abunc<-NULL
for (x in 1:nrow(larve_u))
{
  h<-larvc[larvc$UNIQUE==larve_u$UNIQUE[x],]
  i<-c(sum(h$COUNT_AT_LENGTH),sum(h$LENGTH_ABUNDANCE))
  abunc<-rbind(abunc,i)
}

larve_u$ADJ_ABUNc<-abunc[,2]
larve_u$ADJ_COUNTc<-abunc[,1]
larve_u$ABUNbinc <- factor(ifelse(is.na(larve_u$ADJ_ABUNc) , "0",
                                  ifelse(larve_u$ADJ_ABUNc > 0 & larve_u$ADJ_ABUNc <= 10,"0-10",
                                         ifelse(larve_u$ADJ_ABUNc > 10 & larve_u$ADJ_ABUNc <= 100, '10-100',
                                                '100+'))), 
                           levels=c('0','0-10','10-100','100+'),
                           ordered=T)

#fourth
larvd<-larve[larve$LENGTH>13,]
larvd_n0<-larve[!is.na(larvd$AGE),]

larvd<-merge(larve_u,larvd_n0,all=T)

abund<-NULL
for (x in 1:nrow(larve_u))
{
  h<-larvd[larvd$UNIQUE==larve_u$UNIQUE[x],]
  i<-c(sum(h$COUNT_AT_LENGTH),sum(h$LENGTH_ABUNDANCE))
  abund<-rbind(abund,i)
}

larve_u$ADJ_ABUNd<-abund[,2]
larve_u$ADJ_COUNTd<-abund[,1]
larve_u$ABUNbind <- factor(ifelse(is.na(larve_u$ADJ_ABUNd) , "0",
                                  ifelse(larve_u$ADJ_ABUNd > 0 & larve_u$ADJ_ABUNd <= 10,"0-10",
                                         ifelse(larve_u$ADJ_ABUNd > 10 & larve_u$ADJ_ABUNd <= 100, '10-100',
                                                '100+'))), 
                           levels=c('0','0-10','10-100','100+'),
                           ordered=T)



larve<-larve_u

larve$REGION <- factor(ifelse(larve$PLANKTON_STRATUM=="2" |larve$PLANKTON_STRATUM=="3" |larve$PLANKTON_STRATUM=="5" |larve$PLANKTON_STRATUM=="6", "South",
                                 ifelse(larve$PLANKTON_STRATUM=="9" |larve$PLANKTON_STRATUM=="8" |larve$PLANKTON_STRATUM=="12" |larve$PLANKTON_STRATUM=="11"|larve$PLANKTON_STRATUM=="13","Mid",
                                        "North")), 
                          levels=c("South","Mid","North"),
                          ordered=T)


larve$POSa<-binar(larve$ADJ_ABUNa)
larve$POSb<-binar(larve$ADJ_ABUNb)
larve$POSc<-binar(larve$ADJ_ABUNc)
larve$POSd<-binar(larve$ADJ_ABUNd)

####summary stats

#overall abundance MARMAP VS ECOMON 
larve<-larve[larve$SEASON!='1977'&larve$SEASON!='1988'&larve$SEASON!='1999',]

larve$MEAN_MOs<-factor(ifelse(larve$MEAN_MO==1 | larve$MEAN_MO==2, "Jan-Feb",
                             ifelse(larve$MEAN_MO==3|larve$MEAN_MO==4, "Mar-Apr",
                                    ifelse(larve$MEAN_MO==5|larve$MEAN_MO==6, "May-Jun",
                                           ifelse(larve$MEAN_MO==7|larve$MEAN_MO==8, "Jul-Aug",
                                                  ifelse(larve$MEAN_MO==9|larve$MEAN_MO==10, "Sep-Oct","Nov-Dec"))))),
                       levels=c("Jul-Aug","Sep-Oct","Nov-Dec","Jan-Feb","Mar-Apr","May-Jun"),order=T)

agg1<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~SEASON,larve,sum,na.action=na.pass, na.rm=TRUE)
agg1m<-melt(agg1,id.vars="SEASON")
agg1p<-ggplot(agg1m, aes(x=SEASON,y=value,fill=variable))
agg1p+geom_bar(stat="identity") +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Year", y="Larval Abundance") 

mve<-cbind(agg1[,c(2:5)],factor(c(rep('MARMAP',10),rep('EcoMon',14)),levels=c('MARMAP','EcoMon'),order=T))
colnames(mve)<-c(colnames(agg1[,c(2:5)]),'Survey')
agg2<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~Survey,mve,sum)
agg2m<-melt(agg2,id.vars='Survey')

agg2p<-ggplot(agg2m, aes(x=Survey,y=value,fill=variable))
agg2p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Survey", y="Larval Abundance") 

#now look at overall size by month


agg3<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~MEAN_MOs,larve,sum,na.action=na.pass, na.rm=TRUE)
agg3m<-melt(agg3,id.vars="MEAN_MOs")
agg3p<-ggplot(agg3m, aes(x=MEAN_MOs,y=value,fill=variable))
agg3p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Month", y="Larval Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

larve1<-larve[larve$EVENT_JUL<8000,]
agg4<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~MEAN_MOs,larve1,sum,na.action=na.pass, na.rm=TRUE)
agg4m<-melt(agg4,id.vars="MEAN_MOs")
agg4p<-ggplot(agg4m, aes(x=MEAN_MOs,y=value,fill=variable))
agg4p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Month", y="Larval Abundance",title="MARMAP") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits = c(0, 13600))

#EcoMon

larve2<-larve[larve$EVENT_JUL>8000,]
agg5<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~MEAN_MOs,larve2,sum,na.action=na.pass, na.rm=TRUE)
agg5m<-melt(agg5,id.vars="MEAN_MOs")
agg5p<-ggplot(agg5m, aes(x=MEAN_MOs,y=value,fill=variable))
agg5p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Month", y="Larval Abundance",title="EcoMon") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 13600))

##Now by region

agg6<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~REGION,larve,sum,na.action=na.pass, na.rm=TRUE)
agg6m<-melt(agg6,id.vars="REGION")
agg6p<-ggplot(agg6m, aes(x=REGION,y=value,fill=variable))
agg6p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Region", y="Larval Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

agg7<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~REGION,larve1,sum,na.action=na.pass, na.rm=TRUE)
agg7m<-melt(agg7,id.vars="REGION")
agg7p<-ggplot(agg7m, aes(x=REGION,y=value,fill=variable))
agg7p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Region", y="Larval Abundance",title="MARMAP") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 14200))

#EcoMon


agg8<-aggregate(cbind(ADJ_ABUNa,ADJ_ABUNb,ADJ_ABUNc,ADJ_ABUNd)~REGION,larve2,sum,na.action=na.pass, na.rm=TRUE)
agg8m<-melt(agg8,id.vars="REGION")
agg8p<-ggplot(agg8m, aes(x=REGION,y=value,fill=variable))
agg8p+geom_bar(stat="identity",position=position_dodge()) +
  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
  theme(axis.ticks= element_line(color="black", size=1))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank()) +
  labs(x="Region", y="Larval Abundance",title="EcoMon") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits = c(0, 14200))








                aggp1<-aggregate(cbind(POSa,POSb,POSc,POSd)~SEASON,larve,prop)
                
                colnames(aggp1)<-c('SEASON','POSa','POSb','POSc','POSd')

                aggp1m<-melt(aggp1,id.vars="SEASON")
                aggp1p<-ggplot(aggp1m, aes(x=SEASON,y=value,fill=variable))
                aggp1p+geom_bar(stat="identity") +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Year", y="Proportion Positive")  +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                mve<-cbind(aggp1[,c(2:5)],factor(c(rep('MARMAP',10),rep('EcoMon',14)),levels=c('MARMAP','EcoMon'),order=T))
                colnames(mve)<-c(colnames(aggp1[,c(2:5)]),'Survey')
                aggp2<-aggregate(cbind(POSa,POSb,POSc,POSd)~Survey,mve,prop)
                aggp2m<-melt(aggp2,id.vars='Survey')
                
                aggp2p<-ggplot(aggp2m, aes(x=Survey,y=value,fill=variable))
                aggp2p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Survey", y="Proportion Positive") 
                
                #now look at overall size by month
                
                
                aggp3<-aggregate(cbind(POSa,POSb,POSc,POSd)~MEAN_MOs,larve,prop)


#woo<-aggregate(POSd~MEAN_MOs,larve,prop)
                #aggp3<-data.frame(matrix(ncol = 5, nrow =6))
                #for (x in 1:length(unique(larve$MEAN_MOs)))
                #{
               #   h<-subset(larve,MEAN_MOs==unique(larve$MEAN_MOs)[x])
               #   aggp3[x,1]<-as.character(unique(larve$MEAN_MOs)[x])
               #   aggp3[x,2:5]<-c(sum(h$POSa),sum(h$POSb),sum(h$POSc),sum(h$POSd))/nrow(h)
               # }
               # colnames(aggp3)<-c('MEAN_MOs','POSa','POSb','POSc','POSd')




                aggp3m<-melt(aggp3,id.vars="MEAN_MOs")
                aggp3p<-ggplot(aggp3m, aes(x=MEAN_MOs,y=value,fill=variable))
                aggp3p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Month", y="Propotion Positive") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                larve1<-larve[larve$EVENT_JUL<8000,]
                aggp4<-aggregate(cbind(POSa,POSb,POSc,POSd)~MEAN_MOs,larve1,prop)
                aggp4m<-melt(aggp4,id.vars="MEAN_MOs")
                aggp4p<-ggplot(aggp4m, aes(x=MEAN_MOs,y=value,fill=variable))
                aggp4p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Month", y="Proportion Positive",title="MARMAP") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                  scale_y_continuous(limits = c(0, .181))
                
                #EcoMon
                
                larve2<-larve[larve$EVENT_JUL>8000,]
                aggp5<-aggregate(cbind(POSa,POSb,POSc,POSd)~MEAN_MOs,larve2,prop)
                aggp5m<-melt(aggp5,id.vars="MEAN_MOs")
                aggp5p<-ggplot(aggp5m, aes(x=MEAN_MOs,y=value,fill=variable))
                aggp5p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Month", y="Proportion Positive",title="EcoMon") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#+
#                  scale_y_continuous(limits = c(0, .1))
                
                ##Now by region
                
                aggp6<-aggregate(cbind(POSa,POSb,POSc,POSd)~REGION,larve,prop)
                aggp6m<-melt(aggp6,id.vars="REGION")
                aggp6p<-ggplot(aggp6m, aes(x=REGION,y=value,fill=variable))
                aggp6p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Region", y="Proportion Positive") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                aggp7<-aggregate(cbind(POSa,POSb,POSc,POSd)~REGION,larve1,prop)
                aggp7m<-melt(aggp7,id.vars="REGION")
                aggp7p<-ggplot(aggp7m, aes(x=REGION,y=value,fill=variable))
                aggp7p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Region", y="Proportion Positive",title="MARMAP") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  scale_y_continuous(limits = c(0, .11))
                
                #EcoMon
                
                
                aggp8<-aggregate(cbind(POSa,POSb,POSc,POSd)~REGION,larve2,prop)
                aggp8m<-melt(aggp8,id.vars="REGION")
                aggp8p<-ggplot(aggp8m, aes(x=REGION,y=value,fill=variable))
                aggp8p+geom_bar(stat="identity",position=position_dodge()) +
                  scale_fill_manual(values=c("cyan","brown3","yellow","black"),name="Length",labels=c("1.1-5.3mm", "5.3-8.8mm", "8.8-13mm","13+")) +
                  theme(axis.title.x = element_text(size = rel(2), angle = 00, color="black", vjust=-0.5, face="bold"))+
                  theme(axis.title.y = element_text(size = rel(2), angle = 90, color="black", vjust= 1, face="bold"))+
                  theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29))+
                  theme(axis.text.x = element_text(size=15, face="bold", color="black"))+
                  theme(axis.text.y = element_text(size=15, face="bold", color="black"))+
                  theme(panel.background = element_rect(fill = "white", color = "black", size = 2))+
                  theme(axis.ticks= element_line(color="black", size=1))+
                  theme(panel.grid.major = element_blank())+
                  theme(panel.grid.minor = element_blank()) +
                  labs(x="Region", y="Proportion Positive",title="EcoMon") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                  scale_y_continuous(limits = c(0, .11))

#Make a table showing the number of cruises and tows for each month and year
larve$YRM<-paste(larve$MEAN_MOs,larve$SEASON,sep='_')

#%<% is like saying this 
  #then give the command groupby (it's like aggregates list)
h<-larve%>%
  group_by(SEASON,MEAN_MOs)%>%
    summarize(c_count=length(unique(CRUISE_NAME)))


spread(h,MEAN_MOs,c_count)

###


###
strs_loc<-larve[, names(larve) %in% c('PLANKTON_STRATUM','LATITUDE','LONGITUDE')]
#  #find mean lat and lon by strata
strs_lat<-aggregate(LATITUDE~PLANKTON_STRATUM,strs_loc,mean)
strs_lon<-aggregate(LONGITUDE~PLANKTON_STRATUM,strs_loc,mean)
strs_loc<-merge(strs_lat,strs_lon)
colnames(strs_loc)<-c('strat','LATITUDE','LONGITUDE')
#max(strs_loc$LATITUDE) [1] 41.16227
#min(strs_loc$LATITUDE) [1] 35.87076
#max(strs_loc$LONGITUDE) [1] -69.58234
#min(strs_loc$LONGITUDE) [1] -75.69831


lat <- c(41.5, 35)
lon <- c(-78, -69)

larvmap <- openmap(c(lat[1],lon[1]), c(lat[2],lon[2]), type='osm-bbike')
larvmap <- autoplot.OpenStreetMap(larvmap)


merc <- as.data.frame(projectMercator(larve$LATITUDE,larve$LONGITUDE))
names(merc) <- c('latm','lonm')
larve<-data.frame(larve,merc)



UNIQUE<-unique(larve$YRM)
setwd("C:/Users/Cara/Documents/CBL/GitHub/Chapter2/figures")
for (x in 1:length(UNIQUE))
{
  y<-larve[larve$YRM==UNIQUE[x],]
  y1<-y$MEAN_MOs[1]
  y2<-y$SEASON[1]
  a <- larvmap + geom_point(aes(x=latm, y=lonm,
                                        size=ABUNbina), 
                                    data=larve[larve$YRM==UNIQUE[x],]) +
    labs(title=paste('a',y1, y2)) +
    scale_x_continuous("", breaks = NULL) + 
    scale_y_continuous("", breaks = NULL) +
    theme(legend.position="none") 
  
  b <- larvmap + geom_point(aes(x=latm, y=lonm,
                                size=ABUNbinb), 
                            data=larve[larve$YRM==UNIQUE[x],]) +
    labs(title=paste('b',y1, y2)) +
    scale_x_continuous("", breaks = NULL) + 
    scale_y_continuous("", breaks = NULL) +
    theme(legend.position="none") 
  
  c <- larvmap + geom_point(aes(x=latm, y=lonm,
                                size=ABUNbinc), 
                            data=larve[larve$YRM==UNIQUE[x],]) +
    labs(title=paste('c',y1, y2)) +
    scale_x_continuous("", breaks = NULL) + 
    scale_y_continuous("", breaks = NULL) +
    theme(legend.position="none") 
  
  d <- larvmap + geom_point(aes(x=latm, y=lonm,
                                size=ABUNbind), 
                            data=larve[larve$YRM==UNIQUE[x],]) +
    labs(title=paste('d',y1, y2)) +
    scale_x_continuous("", breaks = NULL) + 
    scale_y_continuous("", breaks = NULL) +
    theme(legend.position="none")
  
  jpeg(filename = paste(x,'map',y1,y2,'.jpg',sep=''), pointsize =12, quality = 200, bg = "white", res = NA, restoreConsole = TRUE)
  multiplot(a,b,c,d, cols=2)
  dev.off()
}

