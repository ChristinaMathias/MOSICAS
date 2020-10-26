#installation des packages si necessaire
#install.packages('httr')
#install.packages('openssl')

#on charge les package

library(httr)
library(openssl)
library(readr)

# liste des varietes disponibles
VARIETES <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21)
VARIETES_DF <- data.frame(VARIETES)
rownames(VARIETES_DF) <-  c("Inconnue", "B 5992", "B 69379", "B 69566", "B 8008", "B 80689", "B 82139", 
														"B 47528", "B 51129", "CO 6415", "R 570", "R 579", "R580", "R581", "R582", "R583", 
														"NCO 376", "Mex 68 / 200", "R577", "R92 804")
VARIETES_SIMULEE <- "R 570"

# liste des stades
STADES <- c(1, 0)
STADES_DF <- data.frame(STADES)
rownames(STADES_DF) <- c("Repousse", "Vierge")
STADE_SIMULE <- "Repousse"

#initialisation du login et password
lgin <- "christina" #login smartis
pwd <- sha512("Math@01071988") #mot de passe entre les "

#creation de l'objet d'authentification
o_auth <- authenticate(lgin, pwd, 'basic') 


#fonction de fransformation web service -> dataframe
ws_transform <- function(x, dataframe = FALSE) {
	data <- content(x)
	data_format <- strsplit(data, '\n')
	x <- strsplit(data_format[[1]], ';')
	nb_col <- nrow(data.frame(x[1][1]))
	d <- matrix(ncol = nb_col, nrow = length(x) - 1)
	nb_ligne <- nrow(d)
	for(i in 2:length(x)) {
		for (j in 1:nb_col) {
			d[i - 1,j] <- x[[i]][j]
		}
	}
	colnames(d) <- x[[1]]
	if (dataframe) {
		df <- as.data.frame(d)
		d <- df
	}
	return(d)
}

call_web_service_ru <- function(long, lat) {
	line_str <- "https://smartis.re/api/WSRU?lat="
	line_str <- paste0(line_str, lat)
	line_str <- paste0(line_str, "&long=")
	line_str <- paste0(line_str, long)
	line_str <- paste0(line_str, "&format=")
	line_str <- paste0(line_str, FORMAT)
	req <- GET(line_str, o_auth)
	return(req)
}

call_web_service_mosiweb <- function(long, lat, ru) {
	line_str <- "https://smartis.re/api/WSMosicas?lat="
	line_str <- paste0(line_str, lat)
	line_str <- paste0(line_str,"&long=")
	line_str <- paste0(line_str, long)
	line_str <- paste0(line_str, "&startdate=")
	line_str <- paste0(line_str, DATE_DEB)
	line_str <- paste0(line_str, "&enddate=")
	line_str <- paste0(line_str, DATE_FIN)
	line_str <- paste0(line_str, "&remplissage=")
	line_str <- paste0(line_str, REMPLISSAGE)
	line_str <- paste0(line_str, "&p0=")
	line_str <- paste0(line_str, P0)
	line_str <- paste0(line_str, "&ru=")
	line_str <- paste0(line_str, ru)
	line_str <- paste0(line_str, "&profrac=")
	line_str <- paste0(line_str, PROF_RAC)
	line_str <- paste0(line_str, "&ressurf=")
	line_str <- paste0(line_str, RES_SURF)
	line_str <- paste0(line_str, "&variete=")
	line_str <- paste0(line_str, CODE_VARIETE_SIMULEE)
	line_str <- paste0(line_str, "&stade=")
	line_str <- paste0(line_str, CODE_STADE_SIMULE)
	line_str <- paste0(line_str, "&irdose=")
	line_str <- paste0(line_str, IR_DOSE)
	line_str <- paste0(line_str, "&freqdose=")
	line_str <- paste0(line_str, FREQ_DOSE)
	line_str <- paste0(line_str, "&freq=")
	line_str <- paste0(line_str, FREQ)
	line_str <- paste0(line_str, "&format=")
	line_str <- paste0(line_str, FORMAT)
	req <- GET(line_str, o_auth)
	return(req)
}

points_pour_mosicas_brut <- read.csv("COORD_SIMUL.csv", sep=";", dec = ".") #au format num_labo, longitude, latitude
points_pour_mosicas_sans_labo <- points_pour_mosicas_brut[,2:3]
point_labo <- points_pour_mosicas_brut[,1]

nb_ligne <- nrow(points_pour_mosicas_sans_labo)


MATOUT <- data.frame(matrix(rep(NA,nb_ligne*length(2002:2020)*73),ncol = 73))

#definition date debut et fin
#au format YYYY-MM-DD

for (YEAR in 2002:2020){
  
DATE_DEB <- paste(YEAR-1,"-09-15",sep='')
DATE_FIN <- paste(YEAR,"-09-15",sep='')
REMPLISSAGE <- 50
P0 <- 1#0.5 # 0.65 # http://www.fao.org/3/X0490e/x0490e0e.htm
PROF_RAC <- 1
RES_SURF <- 5
CODE_VARIETE_SIMULEE <- VARIETES_DF[VARIETES_SIMULEE,]
CODE_STADE_SIMULE <- STADES_DF[STADE_SIMULE,]
IR_DOSE <- 10# 10ll
FREQ_DOSE <- 3 # 3 jours
FREQ <- 1
FORMAT <- "raw"


for (lig in 1:nb_ligne) {
	print(lig)
	long <- points_pour_mosicas_sans_labo$LONG[lig] 
	lat <- points_pour_mosicas_sans_labo$LAT[lig]
	req_ru <- call_web_service_ru(long, lat)
	data_frame_ru <- as.numeric(as.vector(ws_transform(req_ru, TRUE)[[1]]))
	req_mosicas <- call_web_service_mosiweb(long, lat, data_frame_ru)
	data_frame_mosicas <- ws_transform(req_mosicas, TRUE)
	MATOUT[lig+(YEAR-2002)*nb_ligne,2:73] <- as.vector(unlist(c(data_frame_mosicas[dim(data_frame_mosicas)[1],])))
	MATOUT[lig+(YEAR-2002)*nb_ligne,1] <- YEAR
#	output_filename <- paste0("mosiweb_out_", point_labo[lig])
#	output_filename <- paste0(output_filename, ".csv")
#	write.csv(data_frame_mosicas, output_filename,quote=F)
#	write.csv(data_frame_mosicas[dim(data_frame_mosicas)[1],], output_filename,quote=F)
}}

names(MATOUT) <- c('Year',names(data_frame_mosicas))
head(MATOUT)

write.table(MATOUT, paste('SimulationNWS.csv',sep = ''),quote=F,row.names=F,col.names=T,dec='.',sep=';')


MATOUTWS5 <- read.table('SimulationWS5.csv',header=T,sep=';',dec='.',stringsAsFactors = F)
MATOUTWS3 <- read.table('SimulationWS3.csv',header=T,sep=';',dec='.',stringsAsFactors = F)
MATOUTNWS <- read.table('SimulationNWS.csv',header=T,sep=';',dec='.',stringsAsFactors = F)
head(MATOUTWS5)

points_pour_mosicas_brut <- read.csv("COORD_SIMUL.csv", sep=";", dec = ".") #au format num_labo, longitude, latitude
BAS <- rep(points_pour_mosicas_brut$LOC,length(2002:2020))


boxplot(MATOUTWS5$yldcan~BAS)
YEAR <- 2002:2020

MATMOYWS5 <- matrix(rep(NA,length(2002:2020)*5),nrow=5)
MATMOYWS3 <- matrix(rep(NA,length(2002:2020)*5),nrow=5)
MATMOYNWS <- matrix(rep(NA,length(2002:2020)*5),nrow=5)
for (i in 1:5){
  for (j in 1:19){
    YM <- mean(MATOUTWS5$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2006:2020))])
#    YM <- mean(MATOUTWS5$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2019))])
    YMY <- mean(MATOUTWS5$yldcan[(BAS %in% levels(BAS)[i])&(MATOUTWS5$Year ==YEAR[j])])
    MATMOYWS5[i,j] <- (YMY-YM)/YM*100
    YM <- mean(MATOUTWS3$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2006:2020))])
#    YM <- mean(MATOUTWS3$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2019))])
    YMY <- mean(MATOUTWS3$yldcan[(BAS %in% levels(BAS)[i])&(MATOUTWS5$Year ==YEAR[j])])
    MATMOYWS3[i,j] <- (YMY-YM)/YM*100
    YM <- mean(MATOUTNWS$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2006:2020))])
#    YM <- mean(MATOUTNWS$yldcan[BAS %in% levels(BAS)[i]&(MATOUTWS5$Year %in% (2019))])
    YMY <- mean(MATOUTNWS$yldcan[(BAS %in% levels(BAS)[i])&(MATOUTWS5$Year ==YEAR[j])])
    MATMOYNWS[i,j] <- (YMY-YM)/YM*100
  }
}

ORD <- c(2,1,4,5,3)

par(mar=c(4,4,1,1),mar=c(8,4,1,1))
barplot(MATMOYWS5[ORD,19],col=gray(0.2),border=NA,width=1,space=0.5,ylim=c(-25,20),
        ylab='Variation des rendements en 2020 (%) - (2006-2020)',cex.lab=1.2,mgp=c(2,0.5,0),tck=-0.02)
barplot(MATMOYNWS[ORD,19],add=T,col=gray(0.6),border=NA,width=1,space=0.5,axes=F)
legend('topleft',legend=c('Effet temp?rature','Effet stress hydrique'),pch=15,col=gray(c(0.5,0.2)),
       cex=1.2,bty='n')
Bas <- c('Beaufonds','Bois Rouge','Savanna','Grand Bois','Le Gol')
axis(1,at=(1:5)+0.5*(0:4),labels=Bas[ORD],las=2,cex.axis=1.2)

yr <- (2006:2020)

par(mfrow=c(5,1),mar=c(1,5,1,1),oma=c(4,0,0,0))

plot(2006:2020,MATMOYWS5[2,5:19],xlab='',ylab='',xaxt='n',mgp=c(2,0.5,0),tck=-0.02,bty='n',pch=16,cex=1.5,col=gray(0)) ; 
abline(lm(MATMOYWS5[2,5:19]~yr),lwd=2)
legend('topleft',legend=Bas[2],pch=NA,bty='n',cex=1.2)

plot(2006:2020,MATMOYWS5[1,5:19],xlab='',ylab='',xaxt='n',mgp=c(2,0.5,0),tck=-0.02,bty='n',pch=16,cex=1.5,col=gray(0.6)) ; 
abline(lm(MATMOYWS5[1,5:19]~yr),lwd=2)
legend('topleft',legend=Bas[1],pch=NA,bty='n',cex=1.2)

plot(2006:2020,MATMOYWS5[4,5:19],xlab='',ylab='',xaxt='n',mgp=c(2,0.5,0),tck=-0.02,bty='n',pch=16,cex=1.5,col=gray(0.1)) ; 
abline(lm(MATMOYWS5[4,5:19]~yr),lwd=2)
legend('topleft',legend=Bas[4],pch=NA,bty='n',cex=1.2)

plot(2006:2020,MATMOYWS5[5,5:19],xlab='',ylab='',xaxt='n',mgp=c(2,0.5,0),tck=-0.02,bty='n',pch=16,cex=1.5,col=gray(0.4)) ; 
abline(lm(MATMOYWS5[5,5:19]~yr),lwd=2)
legend('topleft',legend=Bas[5],pch=NA,bty='n',cex=1.2)

plot(2006:2020,MATMOYWS5[3,5:19],xlab='',ylab='',xaxt='n',mgp=c(2,0.5,0),tck=-0.02,bty='n',pch=16,cex=1.5,col=gray(0.3)) ; 
abline(lm(MATMOYWS5[3,5:19]~yr),lwd=2)
legend('topleft',legend=Bas[3],pch=NA,bty='n',cex=1.2)

axis(1,at=2006:2020, labels=2006:2020)
mtext('Variation des rendements annuels (%) - (2006-2020)',side=2,cex=1.2,outer=T,padj=2)




barplot(MATMOYWS5[,4:19],beside=T)

barplot(MATMOYWS5[,19],beside=T)
barplot(MATMOYNWS[,19],beside=T)
