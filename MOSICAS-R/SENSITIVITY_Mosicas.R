###################################################
########## mosicas sensitivity analysis ###########
###################################################
rm(list=ls())
library(sensitivity)

############################################################
#### Input and output folder have to be in the same place than "SENSITIVITY_Mosicas.R"
setIni <- setwd(c(dirname(rstudioapi::getActiveDocumentContext()$path)))
setInput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/input",sep='') 
setOutput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/output",sep='') 

# list of parameter to be used in the sensitivity analysis
MATPAR <- read.csv2(paste(setInput,'/Plant_parameter.csv',sep=''),header=T,stringsAsFactors = F,dec='.')
# all parameters
ParnameALL <- MATPAR$paramname[which(MATPAR$varcode=='R570')]
# specific parameters
Parname <- c("humaerdeb","humaerdec","humaerfin","humaertb", "humstemdeb","humstemdec", 
             "kcmax","ke","laicroi","laimax","laitb","laiwksen",
             "p01","p04","pracdeb","pracdec","pracfin", 
             "psilb", "pstrudec", "pstrufin", "pstrutb","pstrutcroi","pstemdeb","pstemdec",
             "pstemfin", "ruemax","ruetk", "ruetopt","sthydbio", "sthydcroi",  
             "taldebtt")
Parname <- ParnameALL
Nbpar <- length(Parname)

# initial value of parameter
Valpar <- MATPAR$plantcrop[which(MATPAR$varcode=='R570')] # plant crop or ratoon crop
Valpar <- Valpar[ParnameALL %in% Parname]

# limits, here +/-5%
Valpar[which(Valpar==0)] <- 0.001 
Bsup <- Valpar*1.05
Binf <- Valpar*0.95
bounds <- apply(cbind(Binf, Bsup), 1, function(x){list(min = x[1], max = x[2])})

# name of the output
ID_SENS = 'All_test2'

###################################################
######### sampling design #########################
### choice between: Moriss, Fast or Sobol methods
Method = 'Fast'  # alternative Moriss / Fast / Variance based analysis
rsensM=100 #default 100
rsensF=350 #default 450
rsensV=2000 #default 1000

if(Method=='Moriss'){
  sa<-morris(model=NULL,
            factors = Parname, 
            r = rsensM, #default 100 
            binf=Binf, 
            bsup=Bsup, 
            design = list(type = "oat", levels = 10, grid.jump = 5))

  # plot sampling
  plot(sa$X)
}

if(Method=='Fast'){
  sa <- fast99(model = NULL,
                factors = Parname,
                n = rsensF,
                q = rep("qunif", Nbpar),
                q.arg=bounds)
  plot(sa$X[,1:2])
}

if(Method=='Variance'){
  library(lhs)
  # we will test only 2nd order 
  Nblimit = factorial(Nbpar)/(factorial(2)*factorial(Nbpar-2))
  if(rsensV < 2*Nblimit){print('Warning, not enough sammpling combinations for interaction analysis')}
  sa <- list()
  #  sa$X <- improvedLHS(rsensV,Nbpar)
  sa$X <- randomLHS(rsensV,Nbpar)
  for(i in 1:Nbpar){ sa$X[,i] <- sa$X[,i]*(Bsup[i]-Binf[i]) + Binf[i]}
  sa$X <- as.data.frame(sa$X)
  names(sa$X) <- Parname
  plot(sa$X[,1:2])
}

###################################################
########## run simulations

Newpar = sa$X
# sensitivity function
source(paste(setIni,'R_MOSICAS/Mosicas_function_sensitivity.R',sep='/'))
# output
y <- MATOUTTEMP

# save
save.image(paste(setOutput,'/SENS_Mosicas_',Method,'_',ID_SENS,sep=''))



###################################################
########## sensitivity 

##### Variance based methods ######################

Method='Variance'
load(paste(setOutput,'/SENS_Mosicas_',Method,'_',ID_SENS,sep=''))

Varout=c('spari','agrfm','stamfm','stamsu')
MATIS <- list() ; CV <-c()

for(j in 1:length(Varout)){

  YSIM <- y[,Varout[j]]
  CV <- c(CV,sd(YSIM)/mean(YSIM))
  hist(YSIM)
  XPAR <- sa$X

  EXP <- paste('YSIM~(',Parname[1],sep='')
  for (i in 2:Nbpar){
    if(i<Nbpar){EXP <- paste(EXP,'+',Parname[i],sep='')}
    if(i==Nbpar){EXP <- paste(EXP,'+',Parname[i],')^2',sep='')}
  }

  fit <- lm(EXP,data=XPAR)
#  require(car)
 # ana <- Anova(fit,type='II')
  ana <- anova(fit)

  SUMSQ <- ana$`Sum Sq`
  FACT <- row.names(ana)

  ### principal sensitivity index
  ISP <- c(SUMSQ[1:31]/sum(SUMSQ),SUMSQ[length(SUMSQ)]/sum(SUMSQ))
  FACTIS <- c(Parname,'Residuals')

  ### total sensitivity index
  IST <- sapply(1:length(Parname), function(i) sum(SUMSQ[grep(Parname[i],FACT)])/sum(SUMSQ))
  IST <- c(IST,ISP[length(ISP)])

  MATIS[[j]] <- matrix(c(ISP,IST-ISP),nrow=2,byrow=T)

}

#x11()

par(mfrow=c(1,4),oma=c(4,6,1,1),mar=c(0,1,0,0))
barplot(MATIS[[1]]*100,names.arg = FACTIS,las=2,horiz=T,mgp=c(2,0.5,0),col=gray(c(0.8,0.1)),
        border=NA,xlim=c(0,50),width=1,space=0.5,ylim=c(1,47),tck=-0.02)
axis(2,at=(1:length(FACTIS))*1.5-0.5=rep('',length(FACTIS)),tck=-0.02)
box();text(25,47,labels=expression(Sigma~'iPAR'),cex=1.2)
barplot(MATIS[[2]]*100,names.arg = rep('',length(FACTIS)),las=2,horiz=T,xlim=c(0,80),
        col=gray(c(0.8,0.1)),border=NA,width=1,space=0.5,ylim=c(1,47),mgp=c(2,0.5,0),tck=-0.02)
axis(2,at=(1:length(FACTIS))*1.5-0.5,labels=rep('',length(FACTIS)),tck=-0.02)
box();text(40,47,labels=expression(ABV[FM]),cex=1.2)
barplot(MATIS[[3]]*100,names.arg = rep('',length(FACTIS)),las=2,horiz=T,xlim=c(0,30),
        col=gray(c(0.8,0.1)),border=NA,width=1,space=0.5,ylim=c(1,47),mgp=c(2,0.5,0),tck=-0.02)
axis(2,at=(1:length(FACTIS))*1.5-0.5,labels=rep('',length(FACTIS)),tck=-0.02)
box();text(15,47,labels=expression(Sta[FM]),cex=1.2)
barplot(MATIS[[4]]*100,names.arg = rep('',length(FACTIS)),las=2,horiz=T,xlim=c(0,30),
        col=gray(c(0.8,0.1)),border=NA,width=1,space=0.5,ylim=c(1,47),mgp=c(2,0.5,0),tck=-0.02)
axis(2,at=(1:length(FACTIS))*1.5-0.5,labels=rep('',length(FACTIS)),tck=-0.02)
box();text(15,47,labels=expression(Sug[DM]),cex=1.2)
mtext('% variance explained',side = 1,line = 2.5,cex=1.2,outer=T)




barplot(IST)
barplot(ISP,add=T,col='black')

barplot((SUMSQ/sum(SUMSQ))[1:62],names.arg = FACT[1:62],las=2)

Varout='agrfm'

tell(sa,y[,Varout])

fit <- lm(sa$y~sa$X[,1]*sa$X[,2])
aa<- anova(fit)
aa$`Sum Sq`/sum(aa$`Sum Sq`)

print(sa)
plot(sa,las=3)


# fast indice calculation
if(Method=='Fast'){
Ifast1 <- sa$D1/sa$V
Ifastt <- 1-sa$Dt/sa$V
}


