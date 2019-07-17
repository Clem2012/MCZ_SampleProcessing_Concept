
rm(list=ls())
source("C:/Users/cg05/Documents/R/Other Scripts/Functions/Multiplot function.R")
library(broom)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggforce)
library(ggfortify)
library(ggrepel)
library(kdensity)
library(ks)
library(shape)
library(RColorBrewer)

#******************************************************************************************
## How much information needed for number of replicates
#******************************************************************************************
ABN<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/ComABN.csv', na.strings= c('',"", NA))
BIO<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/ComBIO.csv', na.strings= c('',"", NA))

## Abundance
ABNdf<-ABN %>% 
  gather(Site, A, WSL.H_A_1_O:HP35_B_3_N) %>% 
  separate(Site, c("Station","Rep", "RepNb", "Design"), "_") %>% 
  unite("Replicate",c("Rep", "RepNb")) %>% 
  filter(Design == "N")

## Biomass
BIOdf<-BIO %>% 
  gather(Site, B, WSL.H_A_1_O:HP35_B_3_N) %>% 
  separate(Site, c("Station","Rep", "RepNb", "Design"), "_") %>% 
  unite("Replicate",c("Rep", "RepNb")) %>% 
  filter(Design == "N")

## Combining
COMdf<-cbind(ABNdf, BIOdf[, 'B'])
colnames(COMdf)[6]<-'B'

COMdf[COMdf$Replicate %in% "A_1", 'Replicate']<- 'R1'
COMdf[COMdf$Replicate %in% "A_2", 'Replicate']<- 'R2'
COMdf[COMdf$Replicate %in% "A_3", 'Replicate']<- 'R3'
COMdf[COMdf$Replicate %in% "B_1", 'Replicate']<- 'R4'
COMdf[COMdf$Replicate %in% "B_2", 'Replicate']<- 'R5'
COMdf[COMdf$Replicate %in% "B_3", 'Replicate']<- 'R6'

## Simulating data based on the real community
## Probability of sampling for each of the taxa per station
PA<- ABNdf %>% 
  group_by(Station) %>% 
  summarise(SumA = sum(A)) %>% 
  left_join(ABNdf, by='Station') %>% 
  group_by(Station, Taxon, SumA) %>% 
  summarise(SumTaxa=sum(A)) %>% 
  mutate(PA=SumTaxa/SumA)

PA<-as.data.frame(PA)

## Number of individual caught on average per replicate 
REP<- ABNdf %>% 
  group_by(Station, Replicate) %>% 
  summarise(sumRep= sum(A))

## Simulation of data for each station
## Choice of the station
for (k in unique(REP$Station)){

## Output dataframe    
  StationSim<-expand.grid(Simulation = 1: 100, Taxon =unique(ABNdf$Taxon), 
                          Replicate = paste0('R',1:6), simA = 0 )

## Number of simulation
    for (l in 1:100){
  cat("sim_count: ", l, "\n")
## Selecting the station-specific probability of sampling taxa  
    PROBA<-PA[PA$Station %in% k,]
      
## Selecting the range of number of individuals that each replicate can collect
    PREP<-data.frame(R = paste0('R', 1:6), REPvalue = round(runif(6, min=min(REP[REP$Station%in% k, 'sumRep']),
                                                                  max(REP[REP$Station%in% k, 'sumRep']))))
                                                                    
## For each replicate
    for (i in PREP$R){
## And the defined number of individuals sampled by replicate
      for (j in 1:PREP[PREP$R %in% i, 'REPvalue']){
## Draw the species according to their probability
        SELECT<-which(cumsum(PA[PA$Station %in% k, 'PA'])>=runif(1))[1]
        Species<-PA[SELECT,'Taxon']
## Individuals are added in the new simulated data
        StationSim[StationSim$Simulation %in% l & 
                     StationSim$Taxon %in% Species &
                     StationSim$Replicate %in% i, 'simA']<- StationSim[StationSim$Simulation %in% l &
                                                                         StationSim$Taxon %in% Species &
                                                                         StationSim$Replicate %in% i, 'simA']+
          length(Species)
        
      }
    }
    }
## That would be too big of a file so we're exporting each simulation
 write.csv(StationSim, paste0('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/output/SIM',k,'.csv'))
}

## Tried with only one simulation to compare with real data
## Simulated data
sumSIM<-StationSim %>% 
  group_by(Station, Taxon) %>% 
  summarise(sumSIMA=sum(simA)) %>% 
  as.data.frame()

colnames(sumSIM)[3]<-'SumTaxa'
sumSIM$Mode<-'Simulation'

## Real data
sumREAL<-PA[, c('Station', 'Taxon', 'SumTaxa')]
sumREAL$Mode<-'Real'

## Comparison
Comp<-rbind(sumSIM, sumREAL)

ggplot(Comp, aes(Taxon, log(SumTaxa), fill=factor(Mode))) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  facet_wrap(~Station)

## Calculating the Shannon entropie difference betweem maximum uncertainty replicate information
## Shannon maximum uncertainty is when each species is equally represented
log2(length(unique(ABNdf$Taxon)))

## Average decrease of uncertainty that the replicate 1 provides
## Selection of replicate 1
ABNrep1<-ABNdf %>% 
  select(Taxon, Station, Replicate, A) %>% 
  filter(Replicate %in% 'R1') %>% 
  spread(Taxon, A)

## formatting
rownames(ABNrep1)<-ABNrep1$Station
ABNrep1[, c('Station', 'Replicate')]<-NULL

## Shannon entropy from vegan
H <- diversity(ABNrep1, base=2)

## Difference max and replicate 1
log2(length(unique(ABNdf$Taxon))) - H

## Evolution of uncertainty after each replicate
COMBIT<-expand.grid(Station = unique(ABNdf$Station), Replicate = paste0('R', 0:6), 
                 NbInd=NA, H = NA, Bits = NA, BitsDatum = NA)

## Matching order with the output
COMBIT<-COMBIT[with(COMBIT, order(as.character(Replicate), as.character(Station))), , drop = FALSE ]

ABNrep<-ABNrep1
ABNrep[,]<-0

REP<-paste0('R', 1:6)

for (i in REP){
  ABNrepTemp<-ABNdf %>% 
    select(Taxon, Station, Replicate, A) %>% 
    filter(Replicate %in% i) %>% 
    spread(Taxon, A)
  rownames(ABNrepTemp)<-ABNrepTemp$Station
  ABNrepTemp[, c('Station', 'Replicate')]<-NULL
  
  ABNrep<-ABNrep+ABNrepTemp
  
  H <- diversity(ABNrep, base=2)
  
  COMBIT[COMBIT$Replicate == i, 'NbInd']<- apply(ABNrep, 1, sum)
  COMBIT[COMBIT$Replicate == i, 'H']<-H
  COMBIT[COMBIT$Replicate == i, 'Bits']<-log2(length(unique(ABNdf$Taxon))) - H
  COMBIT[COMBIT$Replicate == i, 'BitsDatum']<-log2(length(unique(ABNdf$Taxon))) / apply(ABNrep, 1, sum)

}

#*************************************************************************
## Other method but it won't work as there're 'new' categories of species appearing
## There are zeros before which mess up the calculation in the same way as the biomass categories
## Kept it for the options
for (i in REP){
  ABNrepTemp<-ABNdf %>% 
    select(Taxon, Station, Replicate, A) %>% 
    filter(Replicate %in% i) %>% 
    group_by(Station, Replicate) %>% 
    mutate(pi = A/sum(A)) %>% 
    mutate(Hi = ifelse(pi>0, -pi * log2(pi), 0)) %>% 
    group_by(Station) %>% 
    summarise(NbInd = sum(A))
    
  ABNrep<-ABNrep+ABNrepTemp
  
  H <- diversity(ABNrep, base=2)
  
  COMBIT[COMBIT$Replicate == i, 'NbInd']<- apply(ABNrep, 1, sum)
  COMBIT[COMBIT$Replicate == i, 'H']<-H
  COMBIT[COMBIT$Replicate == i, 'Bits']<-log2(length(unique(ABNdf$Taxon))) - H
  COMBIT[COMBIT$Replicate == i, 'BitsDatum']<-log2(length(unique(ABNdf$Taxon))) / apply(ABNrep, 1, sum)
  
}
#*************************************************************************************************

## Defining a starting point
COMBIT[COMBIT$Replicate == "R0", 'NbInd']<-0
COMBIT[COMBIT$Replicate == "R0", 'H']<-0
COMBIT[COMBIT$Replicate == "R0", 'Bits']<-log2(length(unique(ABNdf$Taxon)))
COMBIT[COMBIT$Replicate == "R0", 'BitsDatum']<-0


COMBITdf<-COMBIT %>% 
  gather(Uncertainty, Value, H:BitsDatum)

ggplot(COMBITdf[COMBITdf$Uncertainty %in% c('Bits', 'BitsDatum'),], aes(Replicate, Value, colour=Station, group=Station)) + 
  geom_line()+
  #geom_jitter(colour="grey",alpha=0.1)+
  #stat_density2d(aes(fill=..level.., alpha=..level..),geom='polygon', colour='black', bins=7) + 
  #scale_fill_gradientn(colours=rev(brewer.pal(7,"Spectral")))+
  #geom_smooth(method='gam',alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)+
  scale_y_continuous(name = 'uncertainty')+
  scale_x_discrete(name = 'increasing replicates')+
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Uncertainty, scales = 'free')

## With the simulations from previously
## Need to upload station one by one and calculate the uncertainty value


## Evolution of uncertainty after each replicate
COMBIT<-expand.grid(Simulation = 1:100, Station = unique(ABNdf$Station), Replicate = paste0('R', 0:6), 
                    NbInd=NA, H = NA, Bits = NA, BitsDatum = NA)

## Matching order with the output
COMBIT<-COMBIT[with(COMBIT, order(as.character(Replicate), as.character(Station))), , drop = FALSE ]


REP<-paste0('R', 1:6)

for(k in unique(ABNdf$Station)){
  
  SIMTEMP<-read.csv(paste0('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/output/SIM',k,'.csv'))
  
  for (j in 1:100){
    ABNrep<-ABNrep1[rownames(ABNrep1)%in% k,]
    ABNrep[,]<-0
    
    for (i in REP){
      ABNrepTemp<-SIMTEMP %>% 
        select(Simulation, Taxon, Replicate, simA) %>% 
        filter(Simulation %in% j, Replicate %in% i) %>% 
        spread(Taxon, simA)
      ABNrepTemp[, c('Station', 'Replicate', 'Simulation')]<-NULL
      
      ABNrep<-ABNrep+ABNrepTemp
      
      H <- diversity(ABNrep, base=2)
      
      COMBIT[COMBIT$Station %in% k &
               COMBIT$Simulation %in% j &
               COMBIT$Replicate == i, 'NbInd']<- apply(ABNrep, 1, sum)
      COMBIT[COMBIT$Station %in% k &
               COMBIT$Simulation %in% j &
               COMBIT$Replicate == i, 'H']<-H
      COMBIT[COMBIT$Station %in% k &
               COMBIT$Simulation %in% j &
               COMBIT$Replicate == i, 'Bits']<-log2(length(unique(ABNdf$Taxon))) - H
      COMBIT[COMBIT$Station %in% k &
               COMBIT$Simulation %in% j &
               COMBIT$Replicate == i, 'BitsDatum']<-log2(length(unique(ABNdf$Taxon))) / apply(ABNrep, 1, sum)
      
    }
    
  }
  
}
## Defining a starting point
COMBIT[COMBIT$Replicate == "R0", 'NbInd']<-0
COMBIT[COMBIT$Replicate == "R0", 'H']<-0
COMBIT[COMBIT$Replicate == "R0", 'Bits']<-log2(length(unique(ABNdf$Taxon)))
COMBIT[COMBIT$Replicate == "R0", 'BitsDatum']<-0

COMBITdf<-COMBIT %>% 
  gather(Uncertainty, Value, H:BitsDatum)

ggplot(COMBITdf[COMBITdf$Uncertainty %in% c('Bits'),], aes(Replicate, Value, group=1)) + 
  geom_jitter(colour="grey", alpha=0.1)+
  stat_density2d(aes(fill=..level.., alpha=..level..),geom='polygon', colour='black', bins=7) + 
  scale_fill_gradientn(colours=rev(brewer.pal(7,"Spectral")))+
  geom_smooth(method='loess',alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)+
  scale_y_continuous(name = 'uncertainty')+
  scale_x_discrete(name = 'number of replicates (%)')+
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Station, scales = 'free')



#******************************************************************************************
## How much information needed for individual biomass
#******************************************************************************************
### Importing all bivalve data
BACIA<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACIA.csv', na.strings= c('',"", NA))
BACIB<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACIB.csv', na.strings= c('',"", NA))
BACIC<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACIC.csv', na.strings= c('',"", NA))
BACID<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACID.csv', na.strings= c('',"", NA))
BACIE<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACIE.csv', na.strings= c('',"", NA))

BACI<-rbind(BACIA, BACIB, BACIC, BACID, BACIE)
summary(BACI)

TOTAL<-BACI %>% 
  filter(State == 'Intact' & !is.na(Species)) %>% 
  group_by(Species) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(as.character(Species)) %>% 
  print(n = 92, width = Inf)
  

### Intact Individuals (2228) per bivalve species (76)
ggplot(TOTAL,  aes(x=reorder(Species, -n), y=n))+ 
  geom_bar(stat="identity") +
  scale_x_discrete(name="species")+
  scale_y_continuous(name="number of individuals")+
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  annotate('text', x = 50, y = 400, label = "76 species - 2228 individuals")

ggsave('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Plot/Exploration/BivalveNb.png',
       width = 25, height = 17, units = "cm", dpi = 600)

### Tellina fabula - 587 individuals
INTACT<-BACI %>% 
  filter(State == 'Intact' & !is.na(Species))

Tfab<-INTACT %>% 
  filter(Species == 'Tellina fabula' &  !is.na(Biomass))

Tfab[is.na(Tfab$Biomass)]

ggplot(Tfab, aes(x=Biomass)) + 
  geom_histogram(aes(y=..density..),
                 binwidth=.002,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

#*********************************************************************
# Exploration for loop

TfabG<-Tfab %>% 
  mutate(gr=cut(Biomass, breaks= seq(0, max(Biomass), by = 0.0005))) %>% 
  group_by(gr) %>% 
  summarise(n=n()) %>% 
  mutate(pi = n/587) %>% 
  mutate(Hi = -pi * log2(pi))

TfabG %>% 
  print(n = 51)

log2(105)
# Hmax = 6.714246
sum(TfabG$Hi)

## With the measurement of 587 individuals the uncertainty becomes 4.8524
## This means that the decrease of average uncertainty is (6.714246-4.852382)
log2(105)-sum(TfabG$Hi)

# 1.861864 bits of information

# or

(log2(105)-sum(TfabG$Hi))/587

# 0.0032 bits per datum
#**********************************************************************

## Loop preparation
set.seed(2)

## Output tibble
BIT<-expand.grid(PERC = c(1:100), Simulation = 1:500, NbInd=NA, sumHi = NA, Bits = NA, BitsDatum = NA)

## Range of biomass value (assuming maximum information)
## Not to be changed throughout the simulation (105 bins)
RANGE<-Tfab %>% 
  sample_n(size= round((100*nrow(Tfab))/100)) %>% 
  mutate(gr=cut(Biomass, breaks= seq(0, max(Tfab$Biomass)+0.0004, by = 0.0005))) %>% 
  group_by(gr) %>% 
  summarise(n=n())

## Loop from 100 down to 1% of the total amount of data
for (j in 1:500){
  cat("sim_count: ", j, "\n")
  
  for (i in 1:100){
  Tfab %>% 
    # Sample the % of data
    sample_n(size= round((i*nrow(Tfab))/100), replace=T) %>%  
    
    # Allocate the correct bins of the subset
    mutate(gr=cut(Biomass,                         
                  breaks= seq(0, max(Tfab$Biomass)+0.0004, 
                              by = 0.0005))) -> TEMP
  
  # Relate the fixed number of bins
  TEMP<-full_join(RANGE[, 'gr'], TEMP, by='gr')    
  TEMP %>% 
    # Compute number of rows per bins (if individuals -> 1, if NA -> 0)
      mutate(n = ifelse(is.na(Biomass) == T, 0, 1)) %>%
      group_by(gr) %>% 
      summarise(N=sum(n)) %>%
    # Compute pi (probability of outcome) value for each bin
      mutate(pi = N/round((i*nrow(Tfab))/100)) %>% 
    # Compute Hi (degree of uncertainty) value for each bin
      mutate(Hi = ifelse(pi>0, -pi * log2(pi), 0))->TEMP

  RANGE[,paste0(i, 'Perc')]<- TEMP$pi
  BIT[BIT$PERC == i & BIT$Simulation == j, 'NbInd']<- round((i*nrow(Tfab))/100)
  BIT[BIT$PERC == i & BIT$Simulation == j, 'sumHi']<-sum(TEMP$Hi)
  BIT[BIT$PERC == i & BIT$Simulation == j, 'Bits']<-(log2(105)) - sum(TEMP$Hi)
  BIT[BIT$PERC == i & BIT$Simulation == j, 'BitsDatum']<-((log2(105)) - sum(TEMP$Hi))/sum(TEMP$N)
  } 
}

BITdf<-BIT %>% 
  gather(Uncertainty, Value, sumHi:BitsDatum)

ggplot(BITdf[BITdf$Uncertainty %in% c('Bits', 'BitsDatum'),], aes(PERC, Value)) + 
  geom_jitter(colour="grey",alpha=0.1)+
  stat_density2d(aes(fill=..level.., alpha=..level..),geom='polygon', colour='black', bins=7) + 
  scale_fill_gradientn(colours=rev(brewer.pal(7,"Spectral")))+
  #geom_smooth(method='gam',alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)+
  scale_y_continuous(name = 'Uncertainty')+
  scale_x_continuous(name = 'Number of Individuals (%)')+
  theme_bw()+
  theme(legend.position="none")+
  facet_wrap(~Uncertainty, scales = 'free')

#********************************************************************************************************
## What to do with information on biomass
#********************************************************************************************************
## Finding the biovolume
# Bivalve - assuming cylinder with an elliptic shape
# Ellipse area * length
# (pi * depth * width) * length
Tfab$Biovolume<- (Tfab$Width*Tfab$Depth*pi)*Tfab$Length

ggplot(Tfab, aes(Biomass, Biovolume)) + 
  geom_point(colour="grey")+
  stat_smooth()

bioVmod<-lm(Biomass~Biovolume, data=Tfab)
summary(bioVmod)

## Using the intra-specific variation
BACIA<-read.csv('C:/Users/cg05/Documents/Science/Project - Seedcorn/Proposal/MCZ Trait measurement/Concept/Example data/input/BACIA.csv', na.strings= c('',"", NA))

apply(table(BACIA$StationName, BACIA$Species), 1, sum)

## Station A012 - Tellina fabula
A012<-BACIA[BACIA$StationName %in% "A012",]
summary(A012)

##************************************************************************************
# NOT RUNT - EXAMPLE OF EXTRACTING DATA FROM a ggplot object
g1<-ggplot(A012, aes(x=Biomass, color=Species, fill=Species)) + 
  geom_histogram()+
  geom_density(alpha=0.5)

g2<-ggplot(A012, aes(x=biomass, color=sex, fill=sex)) +
  geom_density(position = "stack", alpha=0.1)

multiplot(g1, g2)

p <- ggplot_build(g2)  # <---- INSTEAD OF `p <- print(m)`
extract<-p$data[[1]]

extract[extract$x > 0.003 & extract$x < 0.0035,]

# NOT RUNT - EXAMPLE OF SUMMING DENSITY PLOT

# Step 1 the data
# Species 1
TPDi1 <- rnorm(5, 10, 0.5)
# Species 2
TPDi2 <- rnorm(6, 9, 0.5)

# Step 2 density estimates for each species from the data
# Set from and to so that x values are the same
TPDs1 <- density(TPDi1, kernel = "gaussian", from = 0, to = 15)
TPDs2 <- density(TPDi2, kernel = "gaussian", from = 0, to = 15)
all.equal(TPDs1$x, TPDs2$x)

plot(TPDs1$x, TPDs1$y, xlim = c(4, 12), type = "l", col = "green")
lines(TPDs2$x, TPDs2$y, col = "blue")

TPDsum <- TPDs1$y + TPDs2$y
lines(TPDs1$x, TPDsum)

# If abundance weights are 0.4 and 0.6
TPDc1 <- 0.4 * TPDs1$y + 0.6 * TPDs2$y
# plot with same x from before
lines(TPDs1$x, TPDc1)

#*************************************************************************************

xlimits=seq(from=-0.05, to=0.07, length=1000)
output<-expand.grid(Taxa = unique(A012$Species), IND=factor(c(1:max(summary(A012$Species)))), xlimits=xlimits, DENSITY=NA)
output <- output[with(output, order(Taxa, IND, xlimits)), , drop = FALSE ]

for(j in unique(A012$Species)){
  IND<-A012[A012$Species %in% j, "Biomass"]
  for (i in 1:length(IND)){
    output[output$Taxa %in% j & output$IND %in% i, "DENSITY"]<-dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01))))/(sum(dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01)))))*length(IND))
  }
}

for(j in unique(A012$Species)){
  IND<-A012[A012$Species %in% j, "Biomass"]
  for (i in 1:length(IND)){
    output[output$Taxa %in% j & output$IND %in% i, "DENSITY"]<-dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01))))/(sum(dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01)))))*(0.05))
  }
}

ggplot(output[output$Taxa %in% 'Tellina fabula',], aes(x=xlimits, y=DENSITY, colour=IND)) +
  geom_line()+
  facet_wrap(~Taxa)

SPdens<-output[output$Taxa %in% 'Tellina fabula',] %>% 
  group_by(xlimits,Taxa) %>% 
  summarise(sumDENS=sum(DENSITY))

ggplot(SPdens, aes(x=xlimits, y=sumDENS)) +
  geom_line()+
  geom_line(data=output[output$Taxa %in% 'Tellina fabula',], 
            aes(x=xlimits, y=DENSITY, colour=IND))+
  scale_y_continuous(name = 'density')+
  scale_x_continuous(name = 'biomass')+
  theme_bw()+
  ggtitle('Tellina fabula')

## Deriving biovolume, respiration and filtering rate
# biovolume = Biomass*1.435e-04 - -9.894e-05
# Filtering rate (in ml/h) F = a*W exp b (Winter 1973)
# Where at 12 deg c (a=2410 & b = 0.73)
# Respiration 0.78–0.88 mg O2/g/h (Landes et al. 2015)
outputTfab<-output[output$Taxa %in% 'Tellina fabula',]

head(outputTfab)
# Can't have negative biomass - need to sort it later
# In the meantime, I'll just move the value across the range
outputTfab$biomass<-outputTfab$xlimits + 0.06
outputTfab$biovolume<-outputTfab$biomass*1.435e-04 - -9.894e-05
outputTfab$filtering<-2410*outputTfab$biomass*exp(0.73)
outputTfab$respiration<-0.83*outputTfab$biomass


outputTfabdf<-outputTfab %>%
  gather(physiology, value, biomass:respiration)

SPdens<-outputTfabdf %>% 
  group_by(physiology, value, Taxa) %>% 
  summarise(sumDENS=sum(DENSITY))


ggplot(SPdens, aes(x=value, y=sumDENS)) +
  geom_line()+
  scale_y_continuous(name = 'density')+
  scale_x_continuous(name = '')+
  theme_bw()+
  facet_wrap(~physiology,  scale = 'free')+
  ggtitle('Tellina fabula')


## Station A017 - Thracia gracilis
A017<-BACIA[BACIA$StationName %in% "A017",]
summary(A017)

xlimits=seq(from=-0.05, to=0.07, length=1000)
output<-expand.grid(Taxa = unique(A017$Species), IND=factor(c(1:max(summary(A017$Species)))), xlimits=xlimits, DENSITY=NA)
output <- output[with(output, order(Taxa, IND, xlimits)), , drop = FALSE ]

for(j in unique(A017$Species)){
  IND<-A017[A017$Species %in% j, "Biomass"]
  for (i in 1:length(IND)){
    output[output$Taxa %in% j & output$IND %in% i, "DENSITY"]<-dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01))))/(sum(dnorm(xlimits, mean=IND[i], sd=ifelse(length(IND)>1, hpi(IND), hpi(c(IND, IND+0.01)))))*length(IND))
  }
}

ggplot(output, aes(x=xlimits, y=DENSITY, colour=IND)) +
  geom_line()+
  facet_wrap(~Taxa)

SPdens2<-output[output$Taxa %in% 'Thracia gracilis',] %>% 
  group_by(xlimits, Taxa) %>% 
  summarise(sumDENS=sum(DENSITY, na.rm=T))

ggplot(SPdens2, aes(x=xlimits, y=sumDENS)) +
  geom_line()+
  geom_line(data=output[output$Taxa %in% 'Thracia gracilis',], 
            aes(x=xlimits, y=DENSITY, colour=IND))+
  scale_y_continuous(name = 'density')+
  scale_x_continuous(name = 'biomass')+
  theme_bw()+
  ggtitle('Thracia gracilis')


## Combining species
VIRTcom<-rbind(SPdens, SPdens2)
head(VIRTcom)

ggplot(VIRTcom, aes(x=xlimits, y=sumDENS, colour=Taxa)) +
  geom_line()

## Proportion of abundance
# 5 Tellina fabula (0.71) and 2 Thracia gracilis (0.29)
VIRTcom[VIRTcom$Taxa %in% "Tellina fabula",'sumDENSw']<-VIRTcom[VIRTcom$Taxa %in% "Tellina fabula", 'sumDENS']*0.71

VIRTcom[VIRTcom$Taxa %in% "Thracia gracilis",'sumDENSw']<-VIRTcom[VIRTcom$Taxa %in% "Thracia gracilis", 'sumDENS']*0.29

# sum of the weighted proportion
sumVIRTcom<-VIRTcom %>% 
  group_by(xlimits) %>% 
  summarise(DENSw = sum(sumDENSw))

ggplot(VIRTcom, aes(x=xlimits, y=sumDENSw, colour=Taxa, fill=Taxa)) +
  geom_line()

test<-cbind(sumVIRTcom, VIRTcom[VIRTcom$Taxa %in% "Tellina fabula", 'sumDENSw'], VIRTcom[VIRTcom$Taxa %in% "Thracia gracilis", 'sumDENSw'])
colnames(test)<-c('Biomass', 'DENSw', 'Tellina', 'Thracia')
head(test)

test$TellinaProp<-(test$Tellina/test$DENSw)*test$DENSw
test$ThraciaProp<-(test$Thracia/test$DENSw)*test$DENSw



ggplot(test, aes(x=Biomass, y=DENSw)) +
  geom_line()+
  scale_y_continuous(name = 'density')+
  scale_x_continuous(name = 'biomass')+
  geom_area(aes(x=Biomass, y=TellinaProp+ThraciaProp, fill='red'))+
  geom_area(aes(x=Biomass, y=TellinaProp, fill='blue'))+
  scale_fill_discrete(name="species",
                      breaks=c("red","blue"),
                      labels=c("Thracia", "Tellina"))
  
  
  
