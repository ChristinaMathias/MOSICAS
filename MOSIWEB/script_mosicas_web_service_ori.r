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

#definition date debut et fin
#au format YYYY-MM-DD

DATE_DEB <- "2017-12-04"
DATE_FIN <- "2018-12-04"
REMPLISSAGE <- 50
P0 <- 0.5
PROF_RAC <- 1
RES_SURF <- 5
CODE_VARIETE_SIMULEE <- VARIETES_DF[VARIETES_SIMULEE,]
CODE_STADE_SIMULE <- STADES_DF[STADE_SIMULE,]
IR_DOSE <- 15
FREQ_DOSE <- 2
FREQ <- 1
FORMAT <- "raw"


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

points_pour_mosicas_brut <- read.csv("D:/Mes Donnees/Admin/Sujet Stage/2020/Florent/mosicas/pt_savanna_gps.csv", sep=";", dec = ".", ) #au format num_labo, longitude, latitude
points_pour_mosicas_brut <- read.csv("pt_savanna_gps.csv", sep=";", dec = ".", ) #au format num_labo, longitude, latitude
points_pour_mosicas_sans_labo <- points_pour_mosicas_brut[,2:3]
point_labo <- points_pour_mosicas_brut[,1]

nb_ligne <- nrow(points_pour_mosicas_sans_labo)

for (lig in 1:2){#nb_ligne) {
	print(lig)
	long <- points_pour_mosicas_sans_labo[lig, 1]
	lat <- points_pour_mosicas_sans_labo[lig, 2]
	req_ru <- call_web_service_ru(long, lat)
	data_frame_ru <- ws_transform(req_ru, TRUE)[[1]]
	req_mosicas <- call_web_service_mosiweb(long, lat, data_frame_ru)
	data_frame_mosicas <- ws_transform(req_mosicas, TRUE)
	#output_filename <- paste0("~/Admin/Sujet Stage/Florent/mosicas/output_p0_0_5_ir_dose_15/mosiweb_out_", point_labo[lig])
	#output_filename <- paste0(output_filename, ".csv")
	#write.csv(data_frame_mosicas, output_filename)
}
