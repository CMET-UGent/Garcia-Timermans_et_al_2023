library(Phenoflow)
library(ggcyto)
library(tidyverse)
library(dplyr)
library(tidyr)
library(kableExtra)
library(ggplot2)
library(lubridate) # better time handling
library(vegan)
library(cowplot)
library(flowClean)
library(flowCore)
library(ggrepel)
library(gganimate)
library(plotly)
library(gapminder)
knitr::opts_chunk$set(echo = TRUE)
options("yaml.eval.expr" = TRUE)

 
Datapath <- c("/Projects2/ThomasPluym/Pilot/Pilot_RepExp/FCS_files_RepExp") # SHINY
Datapath <- c("/Projects2/Cristina/Pilot/FCS_files_RepExp/") 
fcsfiles <- list.files(path=Datapath,recursive=TRUE,pattern=".fcs",
                       full.names=TRUE)
Flowdata <- read.flowSet(files = fcsfiles, transformation = FALSE)

flowData_transformed <- transform(Flowdata,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H","FL3-H","SSC-H","FSC-H")


sqrcutWater<- matrix(c(9.75, 9.75, 15, 15,
                       2.5, 10, 13, 2.5),
                     ncol = 2,
                     nrow = 4)
colnames(sqrcutWater) <- c("FL1-H", "FL3-H")
polyGate1 <- polygonGate(.gate = sqrcutWater, filterId = "Total Cells")
p_scatter1 <- ggcyto::ggcyto(flowData_transformed[sample(20:length(flowData_transformed), 6)], aes(x = `FL1-H`, y = `FL3-H`)) +
  geom_hex(bins = 300) +
  theme_bw() + labs(x = "FL1-H", y = "FL3-H") + geom_gate (polyGate1)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "white")) + coord_cartesian(xlim = c(5, 15), ylim = c(2, 15))
p_scatter1

Flowdata_cells <- (Subset(flowData_transformed, polyGateSG))
Flowdata_death <- (Subset(flowData_transformed, polyGatePI))

Metadata_pilot<-
  data.frame(Sample_names = flowCore::sampleNames (Flowdata),
             do.call(rbind, lapply(strsplit(flowCore::sampleNames (Flowdata), "_"), rbind
             )))
colnames (Metadata_pilot) <- c("Sample_names", "Analysis", "Day", "Time", "Loop", "Sample_point")

A <- flowCore::filter(flowData_transformed, polyGate1)
Live_count <- summary(A);Live_count <-toTable(Live_count)


vol <-
  as.numeric(fsApply(
    flowData_transformed,
    FUN = function(flowData_transformed)
      flowData_transformed@description$`$VOL`
  ))/1000


volume.frame <- data.frame(sample_names = flowCore::sampleNames(flowData_transformed), vol = 
                             as.numeric(fsApply(flowData_transformed,
                                                FUN = function(flowData_transformed)
                                                  flowData_transformed@description$`$VOL`
                             ))/1000)

Live_count_combined <- Live_count %>% 
  # dplyr :: filter(Pressure =="20", Condition == "HPG") %>% 
  mutate(Live_count = true) %>% 
  select(sample, Live_count) 

Metadata_total<- left_join(Live_count_combined, Metadata_pilot, by = c("sample" = "Sample_names"))

# Merge
Metadata <-  left_join (volume.frame, Metadata_total, by = c("sample_names"="sample"))
Metadata$Day[Metadata$Day == "D0"] <- "2022-03-15"
Metadata$Day[Metadata$Day == "D1"] <- "2022-03-16"
Metadata$Day[Metadata$Day == "D2"] <- "2022-03-17"
Metadata$Day[Metadata$Day == "D3"] <- "2022-03-18"
Metadata$Day[Metadata$Day == "D4"] <- "2022-03-19"
Metadata$Day[Metadata$Day == "D5"] <- "2022-03-20"
Metadata$Day[Metadata$Day == "D6"] <- "2022-03-21"
Metadata$Day[Metadata$Day == "D7"] <- "2022-03-22"
Metadata$Time[Metadata$Time == "V"] <- "11:00:00"
Metadata$Time[Metadata$Time == "N"] <- "16:00:00"
Metadata$Sample_point[Metadata$Sample_point == "B.fcs"] <- "Beginning of the loop"
Metadata$Sample_point[Metadata$Sample_point == "M.fcs"] <- "Middle of the loop"
Metadata$Sample_point[Metadata$Sample_point == "E.fcs"] <- "End of the loop"
Metadata$Datetime <- paste (Metadata$Day, Metadata$Time)
Metadata %>% mutate(Datetime = as.POSIXct(Datetime, format= "%Y-%m-%d %H:%M:%OS",tz = "UTC"))
Metadata %>% group_by(Datetime,Loop, vol) %>% 
  summarize(Loopmean = mean(Live_count))


Metadata_loop1 <-Metadata %>% dplyr::filter(Metadata$Loop=="1")
Metadata_loop2 <- Metadata %>% dplyr::filter(Metadata$Loop=="2")
Metadata_loop3 <-Metadata %>% dplyr::filter(Metadata$Loop=="3")

##fancy y axis
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific= TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


L <-  Metadata %>% group_by(Datetime,Loop) %>% 
  summarize(Loopmean = mean(Live_count)) %>% ggplot(aes(x = Datetime, y = Loopmean/25.1*1000)) + geom_line(aes(group=Loop, color=Loop))+
  theme_bw()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9)) +
  geom_point(size = 2, shape = 19, aes(color=`Loop`))+
  scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))+labs(title = "All loops")+
  ylab("Living bacterial cells/mL")+xlab("Day")+ scale_y_continuous(trans="log10")+theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_y_continuous(labels=fancy_scientific)+
  scale_x_discrete(labels= c("Day 1-16h", "Day 2-11h","Day 2-16h", "Day 3-11h","Day 3-16h","Day 4-11h","Day 4-16h","Day 5-11h","Day 5-16h",
                             "Day 6-11h","Day 6-16h","Day 7-11h","Day 7-16h","Day 8-11h","Day 8-16h"))

L +geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                 position=position_dodge(.9))

L1 <- ggplot(Metadata_loop1, aes(x = Datetime, y = Live_count/vol*1000)) + geom_line(aes(group=Sample_point, color=Sample_point))+
  theme_bw()+
  geom_point(size = 2, shape = 19,aes(color=`Sample_point`))+
  scale_color_manual(name= "Sample point", values=c("#105757","#009999", "#84d1d1"))+labs(title = "Loop 1")+
  scale_y_continuous(labels=fancy_scientific)+
  ylab("Living bacterial cells/mL")+xlab("Day")+ scale_y_continuous(trans="log10")+theme(axis.text.x = element_text(angle = 45, hjust=1))
  #facet_grid(`Sample_point`~.)
L1

L2 <- ggplot(Metadata_loop2, aes(x = Datetime, y = Live_count/vol*1000)) + geom_line(aes(group=Sample_point, color=Sample_point))+
  theme_bw()+
  geom_point(size = 2, shape = 19,aes(color=`Sample_point`))+
  scale_color_manual(name= "Sample point", values=c("#6e0b33","#e2005b", "#e882ab"))+labs(title = "Loop 2")+
  scale_y_continuous(labels=fancy_scientific)+
  ylab("Living bacterial cells/mL")+xlab("Day")+ scale_y_continuous(trans="log10")+theme(axis.text.x = element_text(angle = 45, hjust=1))
#facet_grid(`Sample_point`~.)
L2

L3 <- ggplot(Metadata_loop3, aes(x = Datetime, y = Live_count/vol*1000)) + geom_line(aes(group=Sample_point, color=Sample_point))+
  theme_bw()+
  geom_point(size = 2, shape = 19,aes(color=`Sample_point`))+
  scale_color_manual(name= "Sample point", values=c("#054966","#3a7f9c", "#78ccee"))+labs(title = "Loop 3")+
  ylab("Living bacterial cells/mL")+xlab("Day")+ scale_y_continuous(trans="log10")+#theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_y_continuous(labels=fancy_scientific)+
  scale_x_discrete(labels= c("Day 1-16h", "Day 2-11h","Day 2-16h", "Day 3-11h","Day 3-16h","Day 4-11h","Day 4-16h","Day 5-11h","Day 5-16h",
                               "Day 6-11h","Day 6-16h","Day 7-11h","Day 7-16h","Day 8-11h","Day 8-16h"))

#facet_grid(`Sample_point`~.)
L3

Inter1 <- ggplotly(L1)
Inter2 <- ggplotly(L2)
Inter3 <- ggplotly(L3)
htmlwidgets::saveWidget(Inter1, "Loop1_celldensities.html")
htmlwidgets::saveWidget(Inter2, "Loop2_celldensities.html")
htmlwidgets::saveWidget(Inter3, "Loop3_celldensities.html")



summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,"FL1-H"]) #Replace with the column representing the green fluorescence channel (e.g. "FITC-H")
mytrans <- function(x) x/maxval
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)


Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)
beta.div <- (beta_div_fcm(fbasis, ord.type="PCoA"))
#PCA <- (beta_div_fcm(fbasis, ord.type="PCA"))
#NMDS <- (beta_div_fcm(fbasis, ord.type = "NMDS"))
mod_all <- flowFPModel(FCS_resample(flowData_transformed), parameters=c("FL1-H", "FL3-H"), nRecursions=5)
flowfp_all <- flowFP(FCS_resample(flowData_transformed), mod_all)

beta_div_fcm_flowFP <- function(x, d = 4, dist = "bray", k = 2, iter = 100, 
                                ord.type = c("NMDS", "PCoA"), INDICES=NULL,
                                binary = FALSE) {
  x <- x@counts/apply(x@counts, 1, max)
  x <- base::round(x, d)
  if(!is.null(INDICES)){
    x <- by(x, INDICES = INDICES, FUN = colMeans)
    x <- do.call(rbind, x)
  }
  input.dist <- vegan::vegdist(x, method = dist, binary = binary)
  if (ord.type == "PCoA"){ 
    mds.fbasis <- vegan::metaMDS(input.dist, autotransform = FALSE, k,
                                 trymax = iter) 
  } else {
    mds.fbasis <- stats::cmdscale(input.dist, k = 2, eig = TRUE, add = TRUE)
  }
  return(mds.fbasis)
}
# FlowFP fingerprint


mod_all <- flowFPModel(FCS_resample(flowData_transformed), parameters=c("FL1-H", "FL3-H"), nRecursions=5)
flowfp_all <- flowFP(FCS_resample(flowData_transformed), mod_all)
beta.div.flowFP <- beta_div_fcm_flowFP(flowfp_all, ord.type="PCOA")
# Plot beta diversity based on FlowFP
beta.div.co <- data.frame(beta.div.flowFP[["points"]])
colnames(beta.div.co) <- c("Axis1","Axis2")

var <-  round(100*beta.div.flowFP$eig/(sum(beta.div.flowFP$eig)),1)

Metadata_total <- left_join(Metadata, Diversity.fbasis, by = c("sample_names" = "Sample_name"))
Metadata_pilot <- cbind(Metadata_total, beta.div.co)
Metadata_pilot$datetime <- as.POSIXct(paste(Metadata_pilot$Day, Metadata_pilot$Time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC") 
Metadata_beta_loop1 <-Metadata_pilot %>% dplyr::filter(Metadata_pilot$Loop=="1")
Metadata_beta_loop2 <- Metadata_pilot %>% dplyr::filter(Metadata_pilot$Loop=="2")
Metadata_beta_loop3 <-Metadata_pilot %>% dplyr::filter(Metadata_pilot$Loop=="3")


Beta.loop1 <- Metadata_beta_loop1 %>% ggplot() +
  geom_point(aes(x=Axis1, y=Axis2, colour = `Sample_point`),
             size=4,alpha=1, shape=20, show.legend = TRUE) +
  labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),title= "PCoA analysis Loop 1 15/03 - 22/03")+
  coord_fixed(ratio = 1) +
  scale_color_manual(name="Sample point", values= c("#9103AB", "#1A03AB", "#0CBDC6"))+
  #scale_fill_fermenter(name="RO01-TT01 RO Temperatuur in (°C)",
  # palette="YlOrRd",direction = 1)+
  scale_shape_manual(name="Sample point", values= c(21,22,23))+
  theme_bw(base_size = 16, base_family = "")+
  theme(legend.position = "right", panel.grid.major = element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        plot.margin=unit(c(0,0.1,0,0.1),"cm"))+
  coord_fixed(ratio = 1) +
  theme_bw()+
  theme(strip.text.y = element_text(size = 14), axis.text = element_text(size = 14),legend.title = element_text(size = 14),
        legend.text = element_text(size =14), axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"))+theme(legend.position = "right", panel.grid.major = element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
                                                           plot.margin=unit(c(0,0.1,0,0.1),"cm"))+
  guides(fill = guide_legend(override.aes = list(size=4)))

Beta.loop1

Beta.loop1.time <- Metadata_beta_loop1 %>% ggplot(aes(x=Axis1,y=Axis2,fill=`Sample_point`)) +
  geom_point(shape=21,size= 2,alpha=0.6) + labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),title= "PCoA analysis Loop 1")+
  geom_label_repel(aes(label = Datetime,
                       fill = `Day`), color = "white",
                   size = 3.5) +
  theme(legend.position = "bottom")

Beta.loop1.time


Beta.loop2.time <- Metadata_beta_loop2 %>% ggplot(aes(x=Axis1,y=Axis2,fill=`Sample_point`)) +
  geom_point(shape=21,size= 2,alpha=0.6) + labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),title= "PCoA analysis Loop 2")+
  geom_label_repel(aes(label = Datetime,
                       fill = `Day`), color = "white",
                   size = 3.5) +
  theme(legend.position = "bottom")

Beta.loop2.time


Beta.loop3.time <- Metadata_beta_loop3 %>% ggplot(aes(x=Axis1,y=Axis2,fill=`Sample_point`)) +
  geom_point(shape=21,size= 2,alpha=0.6) + labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),title= "PCoA analysis Loop 3")+
  geom_label_repel(aes(label = Datetime,
                       fill = `Day`), color = "white",
                   size = 2.5,segment.size  = 0.2 ) +
  theme(legend.position = "bottom")
Beta.loop3.time

Pilot2 <- Metadata_pilot %>% ggplot() +
  geom_point(aes(x=Axis1, y=Axis2, fill =`Loop`),
             size=3,, shape= 21, alpha=1, show.legend = TRUE) + stat_ellipse(aes(x= `Axis1`, y= `Axis2`, colour=`Day`))+
  labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),title= "PCoA analysis Pilot 15/03 - 22/03")+
  coord_fixed(ratio = 1) +
  scale_fill_manual(name="Loop", values= c("#9103AB", "#1A03AB", "#0CBDC6"))+
  #scale_shape_manual(name="Sample point", values= c(21,22,23))+
  theme_bw(base_size = 16, base_family = "")+
  theme(legend.position = "right", panel.grid.major = element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        plot.margin=unit(c(0,0.1,0,0.1),"cm"))+
  coord_fixed(ratio = 1) +
  theme_bw()+
  theme(strip.text.y = element_text(size = 14), axis.text = element_text(size = 14),legend.title = element_text(size = 14),
        legend.text = element_text(size =14), axis.title = element_text(size=14),
        axis.ticks = element_line(colour = "black"))+theme(legend.position = "right", panel.grid.major = element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
                                                           plot.margin=unit(c(0,0.1,0,0.1),"cm"))+
  guides(fill = guide_legend(override.aes = list(size=4)))

Pilot2
ggplotly(Pilot2)

Bacdiv <-  Metadata_pilot %>% 
  ggplot(aes(x = Datetime, y = D2, shape=Sample_point, fill=Loop))+
  geom_point(size = 3, shape=21)+ # scale_shape_manual(name="Sample point",values = c(21,22,23))+
  scale_fill_manual(name="Loop", values = c("#59B4ED","#9459ED", "red"))+ 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), trans="log10", limits = c(1300,3500))+
  labs(y = "Diversity index (A.U.)", x = "")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("")
  
Bacdiv

BacdivLoop1 <-  Metadata_beta_loop1 %>% 
  ggplot(aes(x = Datetime, y = D2, fill=Sample_point))+
  geom_point(size = 3, shape= 21)+ scale_fill_manual(name="Loop", values = c("#59B4ED","#9459ED", "red"))+ 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), trans="log10", limits = c(1300,3500))+
  labs(y = "Diversity index (A.U.)", x = "")+labs(title = "Bacterial diversity Loop 1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

BacdivLoop1

BacdivLoop2 <-  Metadata_beta_loop2 %>% 
  ggplot(aes(x = Datetime, y = D2, fill=Sample_point))+
  geom_point(size = 3, shape= 21)+ scale_fill_manual(name="Loop", values = c("#59B4ED","#9459ED", "red"))+ 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), trans="log10", limits = c(1300,3500))+
  labs(y = "Diversity index (A.U.)", x = "")+ labs(title = "Bacterial diversity Loop 2")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

BacdivLoop2

BacdivLoop3 <-  Metadata_beta_loop3 %>% 
  ggplot(aes(x = Datetime, y = D2, fill=Sample_point))+
  geom_point(size = 3, shape= 21)+ scale_fill_manual(name="Loop", values = c("#59B4ED","#9459ED", "red"))+ 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), trans="log10", limits = c(1300,3500))+
  labs(y = "Diversity index (A.U.)", x = "", title = "Bacterial diversity Loop 3")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

BacdivLoop3


pH <- read.csv(file="Ph_Loops0.csv",sep=";" )
head(pH)

pH$TimeString <- as.POSIXct(pH$TimeString,format="%d.%m.%Y %H:%M:%S",tz=Sys.timezone())
pH$VarValue <- gsub(",",".", pH$VarValue)
pH$VarValue <- as.numeric(pH$VarValue)

pH_exp <-pH %>% dplyr::filter(pH$TimeString=="")

pHgraph <- pH %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T16:00:04+00:00")) %>% 
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  theme(text = element_text(size = 10))+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))+ #labs(title = "pH over time in the three loops")+
  ylab("pH")+xlab("Time")+ scale_y_continuous(limits=(c(7,9)))+theme(axis.text.x = element_text(angle = 45, hjust=1))
pHgraph

Flow <- read.csv(file="Flow_Loops0.csv",sep=";" )
Flow$VarValue <- gsub(",",".", Flow$VarValue)
Flow$VarValue <- as.numeric(Flow$VarValue)
Flow$TimeString <- as.POSIXct(Flow$TimeString,format="%d.%m.%Y %H:%M:%S",tz=Sys.timezone())
meanFlow <- Flow %>% group_by(VarName)%>%summarise(mean_val=mean(VarValue))

Flow$label <- as.factor(Flow$VarName)
Flow$label<-forcats::fct_collapse(Flow$label, L1 = c("Loop 1_Flow", "Loop 2_Flow", "Loop 3_Flow"))
Fsum <-Rmisc::summarySE(Flow, measurevar="VarValue", groupvars=c("label"))

Flowgraph <- Flow %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T16:00:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))+ #labs(title = "Flow over time in the three loops")+
  ylab("Flowrate (m/s)")+xlab("Time")+ scale_y_continuous(limits=(c(21,28)))+theme(axis.text.x = element_text(angle = 45, hjust=1))
Flowgraph

Cond <- read.csv(file="EC_Loops0.csv",sep=";" )
Cond$VarValue <- gsub(",",".", Cond$VarValue)
Cond$VarValue <- as.numeric(Cond$VarValue)
Condgraph <- Cond %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T16:00:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  #theme(text = element_text(size = 20))+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))+ #labs(title = "Conductivity")+
  ylab("Conductivity µS/cm")+xlab("Time")+ scale_y_continuous(limits=(c(400, 450)))+theme(axis.text.x = element_text(angle = 45, hjust=1))
Condgraph 

Cond$label <- as.factor(Cond$VarName)
Cond$label<-forcats::fct_collapse(Cond$label, all = c("Loop 1_Conductivity", "Loop 2_Conductivity", "Loop 3_Conductivity"))
Csum <-Rmisc::summarySE(Cond, measurevar="VarValue", groupvars=c("label"))


Temp <- read.csv(file="Temperature_Loops0.csv",sep=";" )
Temp$VarValue <- gsub(",",".", Temp$VarValue)
Temp$VarValue <- as.numeric(Temp$VarValue)
Temp$TimeString <- as.POSIXct(Temp$TimeString,format="%d.%m.%Y %H:%M:%S",tz=Sys.timezone())
Temp1 <- Temp %>% dplyr::filter(Temp$VarName=="Temp 1 Loop 1"| Temp$VarName =="Loop 1_Temperature_002"| Temp$VarName =="Temp 3 Loop 1")
Temp2 <- Temp %>% dplyr::filter(Temp$VarName=="Temp 1 Loop 2"| Temp$VarName =="Loop 2_Temperature_002"| Temp$VarName =="Temp 3 Loop 2")
Temp3 <- Temp %>% dplyr::filter(Temp$VarName=="Temp 1 Loop 3"| Temp$VarName =="Loop 3_Temperature_002"| Temp$VarName =="Temp 3 Loop 3")

tgc <- Rmisc::summarySE(tg, measurevar="VarValue", groupvars=c("supp","dose"))
Tsum <-Rmisc::summarySE(Temp, measurevar="VarValue", groupvars=c("VarName"))

Temp$label <- as.factor(Temp$VarName)

Temp$label<-forcats::fct_collapse(Temp$label, L1 = c("Temp 1 Loop 1","Loop 1_Temperature_002","Temp 3 Loop 1"), 
                      L2 = c("Temp 1 Loop 2", "Loop 2_Temperature_002", "Temp 3 Loop 2"),
                      L3 = c("Temp 1 Loop 3","Loop 3_Temperature_002","Temp 3 Loop 3"))

Tsum <-Rmisc::summarySE(Temp, measurevar="VarValue", groupvars=c("label"))

Temp$label<-forcats::fct_collapse(Temp$label, all = c("L1", "L2", "L3"))
Tsum <-Rmisc::summarySE(Temp, measurevar="VarValue", groupvars=c("label"))


Tempgraph <- Temp %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T14:55:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  #scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))
  labs(title = "Temperature")+
  ylab("Temperature °C")+xlab("Time (days)")+ scale_y_continuous(limits=(c(12.5, 20)))+theme(axis.text.x = element_text(angle = 45, hjust=1))
Tempgraph

Temp %>% group_by(VarName,Loop) %>% 
  summarize(Loopmean = mean(Live_count)) %>% ggplot(aes(x = Datetime, y = Loopmean/25.1*1000)) + geom_line(aes(group=Loop, color=Loop))+
  theme_bw()+

TLoop1<- Temp1 %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T14:55:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) +# geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  theme(text = element_text(size = 20),axis.title = element_blank())+
  geom_point(size = 0.1, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days", date_labels=c("0", "1", "2", "3", "4", "5", "6", "7"))+
  scale_color_manual(name= "Measuring point", values=c("#6e0b33","#e2005b", "#e882ab"))+
  labs(title = "Loop 1", x="Day", y= "Temperature °C")+ 
  scale_y_continuous(limits=(c(14, 18))) #+theme(axis.text.x = element_text(angle = 45, hjust=1))
TLoop1 

TLoop2<- Temp2 %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T14:55:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + #geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  theme(text = element_text(size = 20))+
  geom_point(size = 0.1, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days", date_labels=c("0", "1", "2", "3", "4", "5", "6", "7"))+
  scale_color_manual(name= "Measuring point", values=c("#105757","#009999", "#84d1d1"))+
  labs(title = "Loop 2", x="Day", y= "Temperature °C")+ 
  scale_y_continuous(limits=(c(14, 18))) #+theme(axis.text.x = element_text(angle = 45, hjust=1))
TLoop2

TLoop3<- Temp3 %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T14:55:04+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) +# geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  theme(text = element_text(size = 20))+
  geom_point(size = 0.25, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days", date_labels=c("0", "1", "2", "3", "4", "5", "6", "7"))+
  scale_color_manual(name= "Measuring point", values=c("#054966","#3a7f9c", "#78ccee"))+
  labs(title = "Loop 3", x="Day", y= "Temperature °C")+ 
  scale_y_continuous(limits=(c(14, 18))) #+theme(axis.text.x = element_text(angle = 45, hjust=1))
TLoop3

Temperature <- cowplot::plot_grid(TLoop1,TLoop2,TLoop3, ncol=1, align="h")
Temperature

Press <- read.csv(file="Pressure_Loops0.csv",sep=";" )
Press$VarValue <- gsub(",",".", Press$VarValue)
Press$VarValue <- as.numeric(Press$VarValue)
Press$TimeString <- as.POSIXct(Press$TimeString,format="%d.%m.%Y %H:%M:%S",tz=Sys.timezone())

Pressgraph <- Press %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T16:00:00+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  theme_bw()+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  scale_color_manual(name= "Loop", values=c("#e2005b","#009999", "#78ccee"))+labs(title = "Pressure over time in the three loops")+
  ylab("Pressure (bar)")+xlab("Time")+ scale_y_continuous(limits=(c(0.6, 0.9)))+theme(axis.text.x = element_text(angle = 45, hjust=1))
Pressgraph
Pressgraphlines <- Press %>% filter(TimeString > as_datetime("2022-03-15T16:00:04+00:00") & TimeString < as_datetime("2022-03-22T16:00:00+00:00")) %>%
  ggplot(aes(x = TimeString, y = VarValue)) + geom_line(aes(group=VarName, color=VarName))+
  geom_point(size = 0.5, shape = 19,aes(color=`VarName`))+scale_x_datetime(date_breaks = "days")+
  scale_color_manual(name= "Loop", values=c("#C70039","#FF5733", "#FFC300"))+ #labs(title = "Pressure over time in the three loops")+
  ylab("Pressure (bar)")+xlab("Time")+ scale_y_continuous(limits=(c(0.6, 0.9)))+theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  geom_vline(xintercept =as_datetime("2022-03-17T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-17T15:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-19T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-19T15:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-18T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-18T16:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-16T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-16T15:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-20T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-20T15:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-21T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-21T15:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-22T11:00:00+00:00"))+
  geom_vline(xintercept =as_datetime("2022-03-22T15:00:00+00:00"))
  Pressgraphlines




ggsave(plot= pHgraph, filename ="pH_RepExp.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Flowgraph, filename ="Flow_RepExp.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= L, filename ="FCM_allLoops.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= L1, filename ="FCM_Loop1.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= L2, filename ="FCM_Loop2.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= L3, filename ="FCM_Loop3.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Bacdiv, filename ="Bacterialdiversity.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= BacdivLoop1, filename ="BacterialLoop1.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= BacdivLoop2, filename ="BacterialLoop2.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= BacdivLoop3, filename ="BacterialLoop3.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Beta.loop1.time, filename ="PCoALoop1.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Beta.loop2.time, filename ="PCoALoop2.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Beta.loop3.time, filename ="PCoALoop3.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
ggsave(plot= Pilot2, filename ="PCoA.jpeg", device="jpeg", width = 12, height = 7, dpi = 300)
