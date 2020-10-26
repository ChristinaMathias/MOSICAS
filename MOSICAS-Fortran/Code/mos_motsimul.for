c** MOTSIMUL
c** PART1A PRE-INITIALISATION
c** PART1B Reading files sim.txt and trait.txt
c          writing tables trt() irrigxim.txt and obsplantxim.txt      
c** PART1C Reading parcal.txt (calibration variable, parameters names and value)
c          writing par_nom,par_val,par_min,par_max,varipar,par_nbre,par_ecart,par_indice
c** PART1D Opening output files for plant,soil and weather inex. Writting first line

c** Test subroutines
c** PART1E Opening personal output file and writing 1st line   ******
c********* PART2 Begining of SENSITIVITY LOOP (test finparam, calcul nosim)
c** PART2A Update of finparam and paramèters index (CALL MAJ_INDICES_PARAM)
c********* PART3 Begining of TREATMENTS LOOP 
c** PART3A Reading following treatment in table trt()
c** PART3C Reading stamet.txt file
c** PART3D Reading soil general data in solgen.txt file
c** PART3E Reading soil layers data in solcouch.txt file
c** PART3F Reading initial soil layers data in solini.txt file
c** PART3G Calculation of soil waters capacities llw, satw, dulw
c** PART3H Reading irrifation data in irrigxim.txt file and filling of table irr()
c** PART3I Reading Nitrogen fertilisation in fertN.txt file filling of table fertN()
c** PART3J Reading soil observations data in obsol.txt file and filling of table obsol()	
c** PART3K Initialisation of soil parameters values
c** PART3L Reading of plant parameters in plante.txt (call Litplante)
c** PART3M Opening of plant observation file and search of first line of the treatment
c** PART3N variables initialisation
C** PART3O Open weather daily data file and search of first line (beginnng date)
C** PART3P Initialisation of soil water balance parameters 
C** PART3Q Update of parameters values for sensitivity (CALL MAJ_VAL_PARAM)
c** PART3R root initialisation
c********* PART4 Begining of DAILY LOOP
c** PART4A  Calculation of AGE, DATES and Test of daily loop ending
c** PART4B  Reading of daily xeather data
c** PART4C  Test of weather data
c** PART4D  Reading of irrigation data from table irr() and irrigation strategies calculations
c** PART4E  Reading of nitrogen fertilisation data from table fertN()
c** PART4F  Calls of growth and water balance modules
c** PART4G  FORCING with observations of first date
c** PART4H  FORCING with observations of all dates
c** PART4I  Writting observations in dnobs and Simulations results in simtemp 
c** PART4J  Writting outputs in ***.csv
c** PART4K  Writting in perso.csv
c** PART4L  Reading observation data of following day in observation table/file

!******** INITPLANTE, INITPARCAL, MAJ_VAL_PARAM




!***************** PROCEDURES
!******** appelmodul
		
        Subroutine motsimul
        INCLUDE 'cercane1.inc' !JFM000301
	  INCLUDE 'denopo.inc'   !JFM130116
		real valeur,av,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14
		integer n,n1,n2,n3,n4,n5,ios,nosim,ncouch,idlayr,inlayr,decade  !nvaripar,
		integer i,j,k, nforc1,nforct, irotation
		character*20 weather   !JFM000301
		character*20 nomsol,meteo(3),pluvio(3), meteosta !jfm0170524
		integer io1,ioplt,iotrt,iomet,ioplu,trtprec
		integer iom(3), iop(3)
		character fmeteoj1*50,fmeteoj2*50,foutgro*10, typetrt*3
		character c,s1*1000,s2*350,cos*3,mq*3,s3*20 
		character varobsim1*20, varobsimt*20
		character trtname*30,s_err*80,irrname*25    !JFM000301 jfm110920 irrname
	  character datedeb*10,datefin*10,dateact*10,datefich*10
	  integer  ddat1,mdat1,ydat1,ddeb,mdeb,ydeb,dfin,mfin,yfin
	  integer dday,mday,yday
	  character obptrt*25,obpdate*10, obsl1*300
		integer obpyr,obpjdate,obpmm,obpdd, forc1eff
	integer jdatemet(3), jdateplu(3)
	  character typcouchsol*4
	  integer l,nlayr,jdate
	  integer iyrbeg, jdatebeg, iyrend, jdateend  !JFM000301
		real ifuite,itransfert,ip0,iru,ill,idul,isat,ibd,iprofsol		
	  character*20 meteofile, bilanh  !JFM000301
		integer stirr1 !!options irrig
		integer nbmod
		real stirr2,stirr3,stirr4,stirr5,stirr6 !paramètres options irrig
		character meteostr*60  !chaine meteo
		integer metst11,metst21, metst22, nbstamet  !options météo 
		real metco(3),metpl(6) !pondérations météo et pluvio
		real metgr(3) !gradients temp, rayon et etp
	    character*20 metsta(20)
		character*10 metdeb(20),metfin(20)
		real mprecip(3),meeq(3),mtempmx(3),mtempmn(3),msolrad(3) !param stations meteo
		real mux(3),mun(3),mvt(3),mtempmm(3), pprecip(3)  ! param stations meteo
		real amprecip(3),ameeq(3),amtempmx(3),amtempmn(3),amsolrad(3) !param stations meteo
		real amux(3),amun(3),amvt(3), apprecip(3) ,amtempmm(3)  ! param stations meteo
		integer dmet(3),mmet(3),ymet(3),dplu(3),mplu(3),yplu(3)
	    character*3 finmet(3),finplu(3)
		real stalti(3),stlati(3),stlongi(3),sttav(3),stamp(3) ! param stations météo
	! déclarations	ghis 10/4/00
		character*8, allocatable::trait(:) !tableau des noms de traitements

c*************************************************************************************
c************** PART1 GENERALITES ****************************************************

c************** PART1A PRE-INITIALISATION ********************************************
	c=";"
	stirr1=1;stirr2=0;stirr3=0;stirr4=0;stirr5=0;stirr6=0
	perso='non'  !**** Put non ****
	testsub='non' !**** Put non ****
	stopmod=0  ! stop flag
	dstopmod=1000 
	mstopmod=1000 
	parac=0 
	nobs=0; dnobst(:)=-999.; simtempt(:)=-999.; nvobst(:)=""

	
c************************************************************************************
c************** PART1B READ files sim.txt and trait.txt		***************
	nbrecuit=nbrecuit+1
		if (nbrecuit==1) then
		 write (*,*) 'READING : SIM.TXT'
		 write(*,*) 'CHECKING : text files'
		s_err='Lecture fichier impossible:sim.txt'
		open(13,File='sim.txt',err=99) !
		READ(13,*,err=99) s1; read(13,FMT='(A100)') repmeteo ! jfm040304 old READ(13,*) repmeteo	    
		repmeteo=trim(repmeteo)//"\"; repdata=repmeteo
		do i=1,12
	     s2="trait;trstr;obsplant;obsol;irrig;fertn;plante;"
		 s2=trim(adjustl(s2))//"solcouch;solgen;solini;"
	     s2=trim(adjustl(s2))//"metdectemp;stamet"
		 call valcol(i,s2,';'); s2=trim(adjustl(s2))//".txt";iotrt=0
		 write(*,*) trim(repdata)//trim(adjustl(s2))
		 s_err='Erreur Ouverture : '//trim(repdata)//trim(adjustl(s2))	     
		 OPEN(3,File=trim(repdata)//trim(adjustl(s2)),err=99)
		 s_err='Lecture impossible fichier: '//trim(adjustl(s2))
		 READ(3,*,err=99) s1; 						
		 if (i<12) then
		 DO  !scrutation fichier txt
		  read(3,FMT='(A400)',IOSTAT=iotrt) s1;call integrcar(s1,s2)
		 	IF (iotrt<0) EXIT
		 END DO
		 end if
		 close(3)
		end do !do i=1,9
		READ(13,*) s1; read(13,FMT='(A100)') represult ! jfm040304 oldREAD(13,*) represult
		represult=trim(represult)//"\"
		READ(13,*) s1; read(13,FMT='(A3)') typetrt ! jfm150513 
		represult=trim(represult)//"\"		
		READ(13,*) s1
		READ(13,*) foutgro,irotation,bilanh,ifreq,iobs,iwatbal,inbal !jfm000315
		s3=trim(adjustl(bilanh))
		nbmod=0; call nbsep(nbmod,s3,"_");nbmod=nbmod+1
		READ(13,*) s1
	    READ(13,*) typeparam,typecalage !type de simul, type d'optimisation
		open(7,file=trim(represult)//trim(foutgro)//'.txt'
	1		,status='unknown',action='write',access='sequential')
	    write(7,*) 'RAPPORT DE SIMULATION'
		READ(13,*) s3		
	!*************** Reading forcing first date
		IF (adjustl(s3)/="==FORCAGE1") THEN
			write(*,*) "Erreur de lecture FORCAGE1 dans sim.txt"; 
			write(*,*) "Press any key to continue"; read(*,*); stop
		END IF
		READ(13,*) nforc1  !nbre forcing variables 
		if (nforc1>0) then
			do i=1,nforc1
				READ(13,*) varobsim1  !forcing variables
				call verifforc('1',varobsim1)
			end do
		end if
		READ(13,*) s3
	!*************** Reading forcing all dates
		IF (adjustl(s3)/="==FORCAGE2") THEN
			write(*,*) "Erreur de lecture FORCAGE2 dans sim.txt"; 
			write(*,*) "Press any key to continue"; read(*,*); stop
		END IF
		READ(13,*) nforct !nbre forcing variables
		if (nforct>0) then
			do i=1,nforct
				READ(13,*) varobsimt !forcing variables
				call verifforc('t',varobsimt)
			end do
		end if
		IF ((nforct==1) .OR. (nforc1==1)) iobs=1
		IF (ifreq<=0) ifreq=200000 !jfm000316
		call verifsimdat (irotation,iobs,
	1      	iwatbal,inbal,typeparam,nforc1,nforct) 
	!*************** Reading types of irrigation
		READ(13,*) s3
		IF (adjustl(s3)/="==IRRIGATION") THEN
			write(*,*) "Erreur de lecture IRRIGATION dans sim.txt"
			write(*,*) "Press any key to continue"; read(*,*); stop
		END IF
		READ(13,*) stirr1 !Type d'irrigation
	
		If ((stirr1<0).or.(stirr>5)) then
			write(*,*) "Erreur de lecture Type irrigation dans sim.txt"
			write(*,*) "Press any key to continue"; read(*,*); stop
		END IF
		Select case (stirr1)
			case (2,3,4)
			READ(13,*) stirr1, stirr2,stirr3,stirr4
			case (5)
			READ(13,*) stirr1, stirr2,stirr3,stirr4,stirr5
		end select

		write(7,*) 'version:    ',versionmos
		write(7,*) 'fichier     '//adjustr(foutgro)
		write(7,*) 'irotation:  ', irotation
		write(7,*) 'bilanh:     ', bilanh
		write(7,*) ' fréquence: ',ifreq
		write(7,*) 'iobs:       ', iobs
		write(7,*) 'iwatbal:    ',iwatbal
		write(7,*) 'typeparam:  ', typeparam
		write(7,*) 'typecalage: ', typecalage
		write(7,*) 'forcage 1:  ', nforc1
		write(7,*) 'forcage 2:  ', nforct
		write(7,*) 'irrigation: ', stirr1
	

	

	!*************** READ TRAIT.TXT **********************
	!PART1BB
		if (typeparam<2) then
		write (*,*) 'READING :  TRAITEMENTS'
		end if
		READ(13,*) s3
		IF (adjustl(s3)/="==TRAITEMENTS") THEN
			write(*,*) "Erreur de lecture TRAITEMENTS dans sim.txt"
			write(7,*) "Erreur sim.dat: ==TRAITEMENTS n'existe pas"
			write(*,*) "Press any key to continue"; read(*,*); stop
		END IF
		read(13,*) nbtrt
		if (nbtrt>trtmax) Then    !JFM131025
			write(*,*) "Erreur: Nbre trait > ",trtmax
			write(7,*) "Erreur: Nbre trait > ",trtmax
			write(*,*) "Press any key to continue"; read(*,*); stop
		end if	
	    if (typetrt=="obs") then !jfm150513
			OPEN(3,File=trim(repdata)//'trait.txt',err=99); iotrt=0 !jfm000315
			else
			OPEN(3,File=trim(repdata)//'trstr.txt',err=99); iotrt=0 !jfm000315
		end if		
		do i=1,nbtrt
			read(13,*)s1
			REWIND(3)
			READ(3,*) trtname
						
			DO  !scrutation of file trait.txt
				READ (3,*,IOSTAT=iotrt) nbbgi,trtname,
	1			ratoon,cultivar,rowspac,datedeb,datefin,
	2			lati,longi,alti,nomsol,vsolini,profrac,irrname,
     3			meteosta  !jfm20170524
				meteostr="0_1;0_100;0;0_100;0;0;0;0;0_0;0;0" !jfm20170524
				IF (adjustl(trtname)==adjustl(s1)) THEN
				    call datesplit(datedeb,mdeb,ddeb,ydeb) !JFM20100615
				    call datesplit(datefin,mfin,dfin,yfin) !JFM20100615
					trt(i).nbbgi=nbbgi;trt(i).name=trtname
					trt(i).ratoon=ratoon;trt(i).cultivar=cultivar
					trt(i).ecart=rowspac
					call dat_str(trt(i).datedeb,ddeb,mdeb,ydeb)			
					call dat_str(trt(i).datefin,dfin,mfin,yfin)
					trt(i).lati=lati;trt(i).longi=longi
					trt(i).alti=alti;trt(i).nomsol=nomsol
					trt(i).nomsol=nomsol; trt(i).profrac=profrac
					trt(i).vsolini=vsolini;trt(i).iname=irrname  !jfm110920,iname,irrname
					trt(i).metabdeb=0   !JFM120124
					trt(i).meteostr=trim(adjustl(meteostr)) !jfm20170524
					trt(i).meteosta=trim(adjustl(meteosta)) !jfm20170524
										EXIT
				END IF	
				IF (iotrt<0) THEN
				   write(*,*) 'Erreur:'
				   write(*,*) 'pas de traitement '// trim(adjustl(s1))
				   write(7,*) 'pas de traitement '// trim(adjustl(s1))
				  close(3);close(13)
				 write(*,*) "Press any key to continue";read(*,*);stop
				END IF
			END DO !scrutation of file trait.txt
		end do  !do i=1,nbtrt
		close (13)
		CLOSE(3)
	



	!*************** PART1BE  Write file irrigxim.txt ***************
	 if (stirr1==1) then
	 open(13,file=trim(repdata)//'irrigxim.txt',
	1         access='sequential',status='unknown',err=99)
		s1='nomirrig,typirrig,efficience,doseirr,jour,mois,an'  !jfm110920 nomirrig et " "=> ","
	  write(13,fmt="(A49)") s1 
		OPEN(3,File=trim(repdata)//'irrig.txt',err=99); iotrt=0 !jfm000315
		do i=1,nbtrt
		   irrname=trt(i).iname  !jfm110920
	       v1=1      
		   if (i>1) then 
			  do j=1,i-1 
	              if (irrname==trt(j).iname) v1=0	         
			  end do
		   end if
		   if (v1==1) then !jfm110930
			REWIND(3)
			READ(3,*) s1
			DO  !scrutation fichier irrig.txt				
			  READ (3,*,IOSTAT=iotrt) 
	1		  s1,kindirr,effirr,valeur,datefich   !avant 20100615 s1,kindirr,effirr,valeur,ddat1,mdat1,ydat1
			  IF (adjustl(irrname)==adjustl(s1)) THEN !jfm110920 irrname
	            j=len_trim(trtname)  !//" "//kindirr//" "//effirr//" "//valeur
			    call datesplit(datefich,mdat1,ddat1,ydat1)  !jfm20100615
				write(13,"(A20,A1,I1,A1,F4.2,A1,F5.1,2(A1,I2),A1,I4)") 
	1            adjustl(irrname),',',kindirr,',',effirr,',',valeur,
	2            ',',ddat1,',',mdat1,',',ydat1  !jfm110920 irrname
			  END IF
			  s1=''	
			  IF ((iotrt<0)) THEN
				EXIT
			  END IF
			END DO !scrutation fichier irrig.txt
		   end if  !v1==1
		end do  !do i=1,nbtrt
	 close(3)
	 close (13)
	 end if !if (stirr1==1) then


	!*************** PART1BF WRITE FILE obsplantxim.txt ************
	c=','
	!update jfm170312
	!** write first line obsplantxim.txt
	 open(13,file=trim(repdata)//'obsplantxim.txt',
	1         access='sequential',status='unknown',err=99)   !Last modif JFM150220
	s1='obptrt,obpdd,obpmm,obpyr,obnbtigv,obnbtigm,obnbtigfl,obhtvd,'
      s1=trim(s1)//'obnbt,obnbv,obltvdsurf,oblai,obei,obsucre,'
	s1=trim(s1)//'obrdcan,obmstigusi,obhumtig,obsucretig,obfibretig,'
	s1=trim(s1)//'obbrixjtig,puretejtig,obmslv,obmslt,obmfa,obmsa,'
	s1=trim(s1)//'obhumcan,obmsen,obsutil'	
	write(13,fmt="(A223)") trim(s1)  !309
	!** open then read obsplant.txt
	OPEN(3,File=trim(repdata)//'obsplant.txt',err=99); iotrt=0 !jfm000315
	do i=1,nbtrt
	 trtname=trt(i).name
	 REWIND(3)
	 READ(3,*) s1
	 DO  !reading and searching in obsplant.txt
	   READ (3,*,IOSTAT=iotrt) obptrt,datefich,obsutil,
	1   obnbtigv,obnbtigm,obnbtigfl,obhtvd,obnbt,obnbv,
     2   obltvdsurf,oblai,obei,obsucre,obrdcan,obmstigusi,obhumtig,
     3   obsucretig,obfibretig,obbrixjtig,puretejtig,obmslv,obmslt,
     4   obmfa,obmsa,obhumcan,obmsen
	!** Write in obsplantxim.txt good treatment and calage=1
	   IF ((adjustl(trtname)==adjustl(obptrt)).AND.(obsutil==1)) THEN  ! ecrit si meme trtname
		call datesplit(datefich,obpmm,obpdd,obpyr)
		write(13,"(A25,2(A1,I2),A1,I4,23(A1,F9.3),A1,I1)") 
	1	obptrt,c,obpdd,c,obpmm,c,obpyr,c,obnbtigv,c,obnbtigm,c,
	2    obnbtigfl,c,obhtvd,c,obnbt,c,obnbv,c,obltvdsurf,c,
     3    oblai,c,obei,c,obsucre,c,obrdcan,c,obmstigusi,c,obhumtig,c,
     4    obsucretig,c,obfibretig,c,obbrixjtig,c,puretejtig,c,obmslv,c,
     5    obmslt,c,obmfa,c,obmsa,c,obhumcan,c,obmsen,c,obsutil
 	   END IF
	   s1=''	
	   IF ((iotrt<0)) THEN
		 EXIT
	   END IF
	 END DO !scrutation fichier obsplant.txt
	end do  !do i=1,nbtrt
	close(3)
	close (13)
	c=';'
		

	    elaps1=timeF()

		end if  !if (nbrecuit==1) then	JFM061114

c********************* END PART1B reading files sim.txt et trait.txt	***************
c************************************************************************************
c************************************************************************************
c************************************************************************************

c************************************************************************************
c********** PART1C READ parcal.txt (parameters names and values)
	if (nbrecuit==1) then
		par_nom(1:12)="aucun"; par_val(1:12)=0; par_min(1:12)=0 
		par_max(1:12)=0; varipar(1:12)=0; par_nbre(1:12)=0
		varopt(1:12)="NA"; varoptpds(1:12)=0.0  !JFM130117
	end if		
	IF (typeparam>0) THEN  !JFM130116 130118 typeparam>0
		if (typeparam<2) then
			write (*,*) 'LECTURE  FICHIER PARAMETRES : PARCAL.TXT'
		end if
	    OPEN(9,file='PARCAL.txt',err=99)
		rewind(9)
		s_err='Ouverture fichier impossible: parcal.txt'
		IF (typeparam>1) THEN  !JFM140617 enlevé le JFM150227 car structure parcal.txt inadequate
		  read(9,*,err=99) nvaropt          !JFM130116		
	      do i=1,nvaropt  !JFM130116
	        read(9,*) varopt(i), varoptpds(i)    
		  end do
	      varparam=varopt(1) !JFM130116
		end if	    
		read(9,*) nbpar !JFM130116  enlevé le JFM150227 car structure parcal.txt inadequate		
		if (nbpar==0)then
			 s1="sensibilité ou calage: Nombre de parametres >0"
			 write(7,*) "Erreur: nbpar=0"
			 write(*,*) "Erreur: nbpar=0"
			 write(*,*) "Press any key to continue"; read(*,*); stop
		end if
		call verifparam1(typeparam,nbpar,varparam)		
		if (typeparam==1) then !jfm030929
			 do ipar=1,nbpar
				read(9,*)   par_nom(ipar), par_val(ipar), !new
	1						varipar(ipar), par_nbre(ipar)  !new
			 if (par_nom(ipar)=="profrac") parac=1  !JFM080701
			 !write(*,*) par_nom(ipar),par_val(ipar),varipar(ipar), par_nbre(ipar)
			 end do 
		end if !typeparam==1  <2
		if (typeparam>1) then !jfm090930
             if (nbrecuit>=1) then
			do ipar=1,nbpar !jfm030929
				read(9,*)   par_nom(ipar), par_val(ipar) !new
			end do
		   end if
		end if !jfm030929			
		close(9)
	END IF !IF (typeparam>0) THEN

      IF (typeparam==1) THEN
	    do ipar=1,nbpar
		par_ecart(ipar) = par_val(ipar)*varipar(ipar)/100
		par_min(ipar) = par_val(ipar) - par_ecart(ipar)*par_nbre(ipar)
		par_max(ipar) = par_val(ipar) + par_ecart(ipar)*par_nbre(ipar)
		par_indice(ipar)=1
	    par_nbre(ipar)=2*par_nbre(ipar)+1
	    !write(*,*) par_min(ipar),par_max(ipar),par_ecart(ipar)
		end do
			par_indice(1)=0
			!call stat_param("open")	
	  END IF
c********** END PART1C 
c************************************************************************************	




c************************************************************************************
c***** PART1D Opening output files for plant,soil and weather inex. Writting first line
	 c=";"; s1="nosimul"//c//"treatcode"//c//"lati"//c//"longi"
	 s1=trim(s1)//c//"alti"//c
	 open(6,file=trim(represult)//trim(foutgro)//'.csv',
	1         access='sequential',status='unknown',err=99)
	  write(6,"(A34)",ADVANCE="NO") adjustl(adjustr(trim(s1)))	  
	  IF (typeparam==1) then  !si sensibilité
			do i=1,nbpar
			 write(6,"(A16)",ADVANCE="NO") trim(par_nom(i))//c
			end do
	  END IF
	  !pheno
	  s1="dates"//c//"age"//c//"agetu"//c//"glai"//c//"ei"//c//"blagnb"
	  s1=trim(s1)//c//"blatnb"//c//"stagnb"//c//"statnb"//c//"stahtvd"
	  s1=trim(s1)//c//"blastvd"
	  !Stalk Biom+Qual
	  s1=trim(s1)//c//"stamfm"//c//"stamdm"//c//"stamwc"//c//"stamsu"
	  s1=trim(s1)//c//"stamsuc"//c//"stamfbc"
        s1=trim(s1)//c//"stampur"//c//"stambjc"
	  !Above ground Biom+Qual
        s1=trim(s1)//c//"agrfm"//c//"agrdm"//c//"agrwc"
        s1=trim(s1)//c//"bladm"//c//"blagdm"
	  !Observations !JFM150220
	  s1=trim(s1)//c//"ostagnb"//c//"ostadnb"//c//"ostafnb"
	  s1=trim(s1)//c//"ostahtvd"//c//"oblatnb"
	  s1=trim(s1)//c//"oblagnb"//c//"oblastvd"//c//"oglai"
	  s1=trim(s1)//c//"oei"//c//"ostamsu"//c//"ostamfm"//c//"ostamdm"
	  s1=trim(s1)//c//"ostamwc"//c//"ostamsuc"//c//"ostamfbc"
	  s1=trim(s1)//c//"ostambjc"//c//"ostampur"//c//"oblagdm"
	  s1=trim(s1)//c//"obladm"//c//"oagrfm"//c//"oagrdm"
	  s1=trim(s1)//c//"oagrwc"//c//"osendm"
	  !SOIL water-roots Simulated-Observed
	  s1=trim(s1)//c//"hv1"//c//"hv2"//c//"hv3"//c//"hv4"//c//"hv5"
        s1=trim(s1)//c//"hv6"//c//"stock"//c//"stockm"//c//"stockfr"
        s1=trim(s1)//c//"ohv1"//c//"ohv2"//c//"ohv3"//c//"ohv4"
	  s1=trim(s1)//c//"ohv5"//c//"ohv6"//c//"ostock"//c//"ostockfr"
	  s1=trim(s1)//c//"profrac"//c//"drac1"//c//"drac2"
	  s1=trim(s1)//c//"drac3"//c//"drac4"//c//"drac5"//c//"drac6"
	  s1=trim(s1)//c//"precip"//c//"tempmn"//c//"tempmx"//c//"tempmm"
	  s1=trim(s1)//c//"solrad"
	  s1=trim(s1)//c//"eeq"//c//"etm"//c//"etr"//c//"tmp"//c//"trp"
	  s1=trim(s1)//c//"swdef"//c//"swdf1"//c//"swdf2"//c//"swan"
	  s1=trim(s1)//c//"djbl"//c//"djgrow"//c//"pari"//c//"stmo"
	  s1=trim(s1)//c//"spar"//c//"spari"//c//"spr"
	  s1=trim(s1)//c//"sirr"//c//"setp"//c//"setm"//c//"setr"
	  s1=trim(s1)//c//"stmp"//c//"strp"//c//"sddbla"
	  s1=trim(s1)//c//"swdefm"//c//"swdf1m"//c//"swdf2m"//c//"sddgr"
	  s1=trim(s1)//c//"sdrain"//c//"srunoff"//c//"sddblat"//c//"sddgrt"
	  s1=trim(s1)//c//"sparit"//c//"swdeft"//c//"swdf1t"//c//"swdf2t"
	  s1=trim(s1)//c//"stmot"//c//"setmt"//c//"setrt"//c//"srunofft"
	  s1=trim(s1)//c//"sdraint"//c//"setpt"//c//"stmpt"//c//"strpt"
	  s1=trim(s1)//c//"durjour"//c//"srgx"
	  write(6,fmt="(A776)") s1  !762
c***** END PART1D

c**********TEST 1 or several routines ************************************************
	if (testsub=='oui') then
		open(1,file=trim(represult)//'testsub.csv'
	1					,access='sequential',status='unknown')
	end if !if (testsub=='oui')



c************************************************************************************
c***** PART1E PART1E Opening personal output file and writing 1st line	if ((typeparam<2) .AND. (perso=='oui')) then ! 
	if ((typeparam<2) .AND. (perso=='oui')) then ! 
		open(14,file=trim(represult)//'perso.csv'
	1					,access='sequential',status='unknown')
	    s1="codtrait;nomtrait;"
	    write(14,"(A18)",ADVANCE="NO") adjustl(adjustr(trim(s1)))	  
	    IF (typeparam==1) then  !si sensibilité
			do i=1,nbpar
			 write(14,"(A16)",ADVANCE="NO") trim(par_nom(i))//c
			end do
	    END IF
	    s1="dates;age;lai;profrac;stock;stockfr;valrsurf;valrde;valrfe"      !58
		s1=trim(s1)//";precip;irr;eos;es;eeq;etm;etr;tm;tr;drain;kce"        !46
		s1=trim(s1)//";kcp;cstr;swdef;swdf1;swdf2;lai;glait;glaibiom"                           !27
		write(14,FMT="(A150)") trim(adjustl(s1)) !
		av=11.111
	end if
c************************************************************************************



c***** PART1F   Open files Dnobs and Simtemp   ******************
	!IF (typeparam==2) then	 !JFM130118 ==2
	! 	s_err='Ouverture fichier impossible: '//'dnobs' 
	!	OPEN(10,File='dnobs.txt',err=99)                    !JFM130117
	! 	s_err='Ouverture fichier impossible: '//'simtemp' 
	!	OPEN(12,File='simtemp.txt',err=99)                  !JFM130117
	! 	s_err='Ouverture fichier impossible: '//'nvobs' 
	!	OPEN(11,File='nvobs.txt',err=99)                  !JFM130117
	!END IF
c************************************************************************************





c********************************************************************************
c********************************************************************************
c********** PART2 Begining of SENSITIVITY LOOP  *********************************
c********************************************************************************
        finparam = "NON"; nosim=0    
        Do While (finparam == "NON")
	  nosim=nosim+1

c********************************************************************************
c********** PART2A  updates finparam and indices of parameters ******************
				IF (typeparam==1) THEN	    
					CALL MAJ_INDICES_PARAM
				END IF
				IF ((typeparam==0) .OR. (typeparam==2)) THEN
					finparam="OUI"
				END IF
c********************************************************************************



c********************************************************************************
c********************************************************************************		
c********** PART3 Begining of TREATMENTS LOOP  **********************************

c********************************************************************************
c**********  PART3A Reading following treatment in table trt()
	DO notrait=1,nbtrt
	trtname=trt(notrait).name
	if ((nosim==1) .and. (typeparam<2)) then  !jfm030929
		write(7,*) ' '
		write(7,*) repeat('*',40)
		write(7,*) 'Traitement :   '//adjustl(trtname)
		write(7,*) repeat('*',40)
		write(*,*) notrait, nbtrt," ",  adjustl(trtname)
	end if
	meteofile=trt(notrait).meteo(1)
	call datesplit(trt(notrait).datedeb,mdeb,ddeb,ydeb); iyrbeg=ydeb
	dday=ddeb; 	mday=mdeb; yday=ydeb
	jdatebeg=julday(mdeb,ddeb,ydeb)
	call datesplit(trt(notrait).datefin,mfin,dfin,yfin); iyrend=yfin	
	jdateend=julday(mfin,dfin,yfin)
	njtot=jdateend-jdatebeg
	nbbgi=trt(notrait).nbbgi
	ratoon=trt(notrait).ratoon
	rowspac=trt(notrait).ecart
	cultivar=trt(notrait).cultivar
	lati=trt(notrait).lati;longi=trt(notrait).longi
	alti=trt(notrait).alti
	nomsol=trt(notrait).nomsol
	if (parac==0) profrac= trt(notrait).profrac !jfm080701 pas de calage,sensibilité de profrac
	vsolini=trt(notrait).vsolini; irrname=trt(notrait).iname   !jfm110920;iname;irrname
	meteosta=trt(notrait).meteosta !jfm20170524
	meteostr=trt(notrait).meteostr !jfm20170524
c********************************************************************************

	
c********************************************************************************
c*********** PART3B writing in output file **.txt
	if ((nosim==1) .and. (typeparam<2)) then  !jfm120124
		write(7,*) 'Bg plantes /mlin   ',nbbgi
		write(7,*) 'Cycle              ',ratoon
		write(7,*) 'Variété            ',cultivar
		write(7,*) 'Ecartement         ',rowspac
		write(7,*) 'Date début         ', trt(notrait).datedeb
		write(7,*) 'Date fin           ', trt(notrait).datefin
		write(7,*) 'Altitude           ' , alti
		write(7,*) 'Latitude           ',lati
		write(7,*) 'Longitude          ',longi
		write(7,*) 'Nom Sol            ',nomsol
		write(7,*) 'Sol Initial        ',vsolini
		write(7,*) 'Prof. racinaire    ',profrac
		write(7,*) 'Nom Irrigation     ',irrname
	    write(7,*) 'Options meteo      ', meteostr
		write(7,*) 'Meteo 1            ',meteosta  !jfm20170524
	end if
	s1= trt(notrait).datedeb // ' ' // trt(notrait).datedeb
c********************************************************************************



c********************************************************************************
c*********** PART3C READ stamet.txt *******
	IF (adjustl(meteostr)/="0_") THEN  !jfm20170524
	i=1
		stalti(1)=-999; stlati(1)=-999; stlongi(1)=-999
		sttav(1)=-999; stamp(1)=-999
	s_err='Ouverture fichier impossible:stamet.txt'	
		if (meteosta=="0") EXIT
		i=1
		OPEN(3,File=trim(repmeteo)//'stamet.txt',err=99); iotrt=0 !jfm000315
			READ(3,*) s1
			DO  !scrutation stamet.txt
			  READ (3,*,IOSTAT=iotrt) s1,v1,v2,v3,v4,v5 !v1=alti,v2=lati,v3=longi,v4=tav,v5=amp
			  IF (trim(adjustl(s1))==meteosta) THEN
			   stalti(1)=v1; stlati(1)=v2; stlongi(1)=v3
			   sttav(1)=v4; stamp(1)=v5;EXIT	
			  END IF
			  IF ((iotrt<0)) THEN
				  write(*,*) trtname,'  Erreur dans stamet.txt'
				  write(*,*) meteosta, ' inexistant'  
				  write (7,*) 'Erreur dans stamet, pas de ',meteosta !
				  close(3)
			   write(*,*) "Press any key to continue"; read(*,*); stop
			  END IF
			END DO
		CLOSE(3)
	IF (meteostr=="0_") EXIT
	END IF  !	IF (adjustl(meteostr)/="0_") THEN
c********************************************************************************


c********************************************************************************
c********** PART3D Reading soil general data in solgen.txt file *************
		OPEN(3,File=trim(repdata)//'solgen.txt',err=99); iotrt=0 !jfm000315
			REWIND(3)
			READ(3,*) s1
			DO  !scrutation fichier solprof.txt
				READ(3,*,IOSTAT=iotrt)s1,typsol,nlayr,ip0,iru,iprofsol !jfm040304
				IF (adjustl(nomsol)==adjustl(s1)) THEN
				  p0=ip0;				  
				  if (nlayr==1) then
					 if (vsolini==-1) then
					   vsolini=0.25
					   write(7,*)'vsolini=-1 abherant car monocouche'
					 end if
					 typcouchsol="mono"
				   else
					 typcouchsol="mult"
				  end if
				  EXIT
				END IF
!				s1=''	
				IF (iotrt<0) THEN
				  write(*,*) 'Erreur: pas de profil de sol '// s1(1:8)
				  !DEALLOCATE(obsol)
				  close(3);close(13)
				  write(*,*)"Press any key to continue";read(*,*);stop
				END IF
			END DO !scrutation of file solprof.txt
		CLOSE(3)
		!****Discrétisation in layers if soil is mono layer
		if (typcouchsol=="mono")  then 
		    nlayr=6
			DO i=1,6
	            if (i==1) then
				    dlayr(i)=5
				  else
					if ((iprofsol-5)<=75) then
						dlayr(i)=(iprofsol-5)/5
					  else
						if (i==2) dlayr(i)=15
						if (i>2)  dlayr(i)=(iprofsol-20)/4
					end if
				end if  !if (i==1) then
				ll(i)=0.3;dul(i)=ll(i)+iru/1000;sat(i)=dul(i)+0.05
				bd(i)=1	
				!write(*,*) i,dlayr(i),ll(i),dul(i)			
			END DO
			if (profrac>iprofsol) profrac=iprofsol   !jfm080701
		end if
c********************************************************************************



c********************************************************************************
c********** PART3E Reading soil layers data in solcouch.txt file *****************
		if (typcouchsol=="mult") then
		OPEN(3,File=trim(repdata)//'solcouch.txt',err=99); iotrt=0 !jfm000315
			REWIND(3)
			READ(3,*) s1
			i=0;j=0 !;nlayr=0
			DO  !scrutation file solcouch.txt
				READ (3,*,IOSTAT=iotrt) s1,ncouch,idlayr,
	1					ill,idul,isat,ibd
				IF (adjustl(nomsol)==adjustl(s1)) THEN
					i=i+1
					j=j+1 !; nlayr=nlayr+1
					dlayr(ncouch)=idlayr;ll(ncouch)=ill
					dul(ncouch)=idul;sat(ncouch)=isat;bd(ncouch)=ibd
				END IF	
				s1=''
				!<====== Vérification layers
				IF ((iotrt<0)) THEN
				 if (i==0) then
					write(*,*) trtname, nomsol
					write(*,*)'Erreur: pas de couche de sol '//s1(1:8)
					write (7,*) 'solcouch'
					write(7,*) trim(adjustl(nomsol)) // ': introuvable'
					close(3);close(13)
				  write(*,*)"Press any key to continue";read(*,*);stop
				 END IF
					EXIT
				END IF
			END DO !scrutation file solcouch.txt
			if (nosim==1) write(7,*) 'solcouch.txt lu'
		CLOSE(3)
		end if !nlayr>1
c********** End PART3E 
c********************************************************************************
	!write(*,*) "5 ",timeF()

c********************************************************************************
c********** PART3F Reading initial soil layers data in solini.txt file **********
		IF ((vsolini<0) .and. typcouchsol=="mult") THEN
		OPEN(3,File=trim(repdata)//'solini.txt',err=99); iotrt=0 !jfm000315
			REWIND(3)
			READ(3,*) s1
			i=0;j=0;inlayr=0
			DO  !scrutation file solcouch.txt
				READ (3,*,IOSTAT=iotrt) s1,ncouch,v1
				IF (adjustl(trtname)==adjustl(s1)) THEN
				  i=i+1
				  IF (ncouch>0) THEN
					 j=j+1; inlayr=inlayr+1
					 sw(ncouch)=v1
				   ELSE
					 write(*,*) trtname, nomsol
				     write(*,*)'Erreur solini : No couche = 0 '
					 close(3);close(13)
				     write(*,*)"Press any key to continue";read(*,*)
					 stop
				  END IF
				END IF
				s1=''	
				!<====== Vérifications Layers
				IF ((iotrt<0)) THEN
				 if (i==0) then
					write(*,*) trtname, nomsol
				write(*,*)'Erreur: pas de couches initiales '//s1(1:8)
					write (7,*) 'solini'
					write(7,*) trim(adjustl(nomsol)) // ': introuvable'
					!DEALLOCATE(obsol)
					close(3);close(13)
				  write(*,*)"Press any key to continue";read(*,*);stop
				 END IF
					EXIT
				END IF
			END DO !scrutation file solcouch.txt
				if (nosim==1) write(7,*) 'solini.txt lu'
		CLOSE(3)
		ELSE	!(vsolini>=0)  THEN
			if (vsolini==0.) vsolini=0.05
			DO i=1,nlayr
				sw(i)=ll(i)+(dul(i)-ll(i))*vsolini
			END DO
		END IF 	!(vsolini<0)  THEN
c********** End PART3F 
c********************************************************************************
	

c********************************************************************************
c********** PART3G Calculation of soil waters capacities llw, satw, dulw ********
				satw=0;dulw=0;llw=0;tsw=0.;tll=0.; stockfdeb=0
			do l=1,nlayr
				esw(l)=dul(l)-ll(l)
	!write(*,*) l, dul(l), ll(l)
			if (esw(l)<=0) then
			write(7,*) 'Erreur couche ',l
			write(7,*) 'dul <= ll'
			write(*,*) adjustl(trtname)," erreur calculs stocks d'eau"
			write(*,*)"Press any key to continue";read(*,*);stop
			end if
				tll=tll+ll(l)*dlayr(l)
				tsw=tsw+(sw(l)-ll(l))/esw(l)
				stockfdeb=stockfdeb+sw(l)*dlayr(l)*10
				llw=llw+ll(l)*dlayr(l)*10
				satw=satw+sat(l)*dlayr(l)*10
				dulw=dulw+dul(l)*dlayr(l)*10
			end do	
			!write(*,*) "toto",llw,dulw,iru	
		if (nosim==1) write(7,*) 'Calcul stocks'
c********** E,d PART3G
c********************************************************************************



c********************************************************************************
c********* PART3H Reading irrigation data and filling of table irr()
		!** détermination of irrigations number nirr
		nirr=0
		if (stirr1==1) then
		OPEN(3,File=trim(repdata)//'irrigxim.txt',err=99); iotrt=0 !JFM061118
			REWIND(3)
			READ(3,*) s1
			DO  !scrutation of file irrigxim.txt and détermination of nirr
				READ (3,*,IOSTAT=iotrt) s1
				IF (adjustl(irrname)==adjustl(s1)) THEN   !jfm110920 trtname=>irrname
					nirr=nirr+1
				END IF
				s1=''	
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation of file irrigxim.txt
		!** Filling table irr()
			IF (nirr>0) THEN
			REWIND(3);iotrt=0;l=0
			READ(3,*) s1
			DO  !scrutation fichier irrig.txt
			  READ (3,"(A20,A1,I1,A1,F4.2,A1,F5.1,2(A1,I2),A1,I4)"
	1		  ,IOSTAT=iotrt) 
	2		  s1,c,kindirr,c,effirr,c,valeur,c,ddat1,c,mdat1,c,ydat1  !JFM061118
				IF (adjustl(irrname)==adjustl(s1)) THEN    !jfm110920 trtname=>irrname
					l=l+1;irr(l).irrdose=valeur;irr(l).irryr=ydat1
					irr(l).irrday=julday(mdat1,ddat1,ydat1)
					irr(l).effirr=effirr
					irr(l).kindirr=kindirr !JFM060624
					!write (*,*) kindirr,effirr,valeur
				END IF
				s1=''	
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation fichier irrig.txt
			END IF !(nirr>0 THEN)
		CLOSE(3)
		if (nosim==1) write(7,*) 'irrigations ', nirr
	end if !if (stirr==1) then
!*********  END PART3H
c********************************************************************************


c********************************************************************************
!********* PART3I Reading Nitrogen fertilisations in fertN.txt and filling table fertN()
		nfertN=0
		OPEN(3,File=trim(repdata)//'fertN.txt',err=99); iotrt=0 !jfm000315
			REWIND(3)
			READ(3,*) s1
			DO  !scrutation file fertN.txt
				READ (3,*,IOSTAT=iotrt) s1
				IF (adjustl(trtname)==adjustl(s1)) THEN
					nfertN=nfertN+1
				END IF
				s1=''	
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation file fertN.txt
		!** filling table fertN()
			IF (nfertN>0) THEN
			REWIND(3);iotrt=0;l=0
			READ(3,*) s1
			DO  !scrutation file fertN.txt
				READ (3,*,IOSTAT=iotrt) 
	1			s1,idlayr,valeur,datefich   !avant 20100615 s1,idlayr,valeur,ddat1,mdat1,ydat1
				IF (adjustl(trtname)==adjustl(s1)) THEN
					l=l+1;fertN(l).fertNdose=valeur
					call datesplit(datefich,mdat1,ddat1,ydat1)  !New
					fertN(l).dlayN=idlayr;fertN(l).fertNyr=ydat1
					fertN(l).fertNday=julday(mdat1,ddat1,ydat1)
				END IF
				s1=''	
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation file fertN.txt
			END IF !(nfertN>0 THEN)
		CLOSE(3)
		if (nosim==1) write(7,*) 'Fertilisations N ', nfertN
c*********  END PART3I
c********************************************************************************

	!write(*,*) "9 ",timeF()
!
c********************************************************************************
c********* PART3J Reading soil observations in obsol.txt and filling table obsol()
		nobsol=0
		OPEN(3,File=trim(repdata)//'obsol.txt',err=99); iotrt=0 !jfm000315
			REWIND(3)
			READ(3,*) s1
			DO  !scrutation file obsol.txt
				READ (3,*,IOSTAT=iotrt) s1
				IF (adjustl(trtname)==adjustl(s1)) THEN
					nobsol=nobsol+1
				END IF
				s1=''
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation file obsol.txt
				!<====== Vérifications dates  AFAIRE
		!** filling table obsol()
			IF (nobsol>0) THEN
			REWIND(3);iotrt=0;l=0
			READ(3,*) s1
			DO  !scrutation file obsol.txt
				READ (3,*,IOSTAT=iotrt) 
	1			s1,datefich,v1,v2,v3,v4,v5,v6
				IF (adjustl(trtname)==adjustl(s1)) THEN
				call datesplit(datefich,mdat1,ddat1,ydat1) !JFM20100615
				l=l+1;obsol(l).swoyr=ydat1
				obsol(l).swoday=julday(mdat1,ddat1,ydat1)
				obsol(l).swo(1)=v1;obsol(l).swo(2)=v2
				obsol(l).swo(3)=v3;obsol(l).swo(4)=v4
				obsol(l).swo(5)=v5;obsol(l).swo(6)=v6
				END IF	
				s1=''
				IF ((iotrt<0)) THEN
					EXIT
				END IF
			END DO !scrutation file obsol.txt
			END IF !(nobsol>0 THEN)
		CLOSE(3)
		if (nosim==1) write(7,*) 'observ. sol : ', nobsol
!*********  END PART3J
c********************************************************************************
						
c********************************************************************************
c****** PART3K Initialisation of soil parameters values
C
      u=5.0
      alf=3.5                  
      icsdur=0  
      salb=0.12
      lat=45.
c********************************************************************************	



c********************************************************************************
c******* Reading of plant parameters in plante.txt (call Litplante)
	psilb=8.
	CALL LITPLANTE   !jfm
	if (nosim==1) write(7,*) 'Lecture plante.txt effectuee'
c********************************************************************************



c********************************************************************************
c******* PART3M Open plant observation file and search first line of the treatment
		s_err='Ouverture fichier impossible: '//'obsplant.prn'
		OPEN(4,File=trim(repdata)//'obsplantxim.txt',err=99); ioplt=0 !Last modif JFM150220
	  REWIND(4)
		READ(4,*) obsl1
		DO
		read(4,"(A25,2(A1,I2),A1,I4,23(A1,F9.3),A1,I1)",IOSTAT=ioplt) 
	1	obptrt,c,obpdd,c,obpmm,c,obpyr,c,obnbtigv,c,obnbtigm,c,
	2    obnbtigfl,c,obhtvd,c,obnbt,c,obnbv,c,obltvdsurf,c,
     3    oblai,c,obei,c,obsucre,c,obrdcan,c,obmstigusi,c,obhumtig,c,
     4    obsucretig,c,obfibretig,c,obbrixjtig,c,puretejtig,c,obmslv,c,
     5    obmslt,c,obmfa,c,obmsa,c,obhumcan,c,obmsen,c,obsutil
			IF (adjustl(obptrt)==adjustl(trtname))  THEN
			  !old obpjdate=julday(obpmm,obpdd,obpyr)-julday(1,1,obpyr)+1
			  obpjdate=julday(obpmm,obpdd,obpyr)
			  !write(*,*) obpdd, obpmm, obpyr		
			  EXIT
			END IF				
			IF (ioplt<0) EXIT
		END DO
c********************************************************************************	



c********************************************************************************      
c********* PART3N Initialisation calculated variables *********************
	  !1******. évaporation soil	
        sumes1=100.*(1-tsw) ! ; write(*,*) tsw, sumes1
        sumes2=0.
        tsw=0.
        swef=0.9-0.00038*(dlayr(1)-30)**2
	    !****** paramètres
	    Wver=0.3; K2ver=0.6  !paramètres des vertisols
	  !2******. innitialisation plant growth
		sirr=0; stmo=0;spar=0;spari=0;sdj=0;spr=0;sparim=0
		srgo=0;srgx=0; sdj1=0; htvdsdj=0
	    swdefsat=1
		sswdf1=0;sswdf2=0; sswdef=0; sdrain=0 ; srunoff=0; sevap=0
		setp=0;setm=0;setr=0; stmp=0; strp=0; sdraint=0; srunofft=0.
		setpt=0.;stmpt=0.;strpt=0.
		setmt=0.;setrt=0.
		aNitr=0;dNitr=0
		!*Stress
		cske=1;swdf1=1;swdf2=1;csrue=1;swdef=1
		!*blades
		blnbtot=0;blnbv=0;blsdj=0;bltvdsurf=0.;bltigsurf=0
		bltvdlong=0; bltvdlarg=0
	    !*stalks
	    stlev=0;sttal=0;talsdj=0.
		nbtigv=0.;nbtigt=0.;nbtigvl=0.;nbtigtl=0.;talsdj=0
		nbtigptl=0.;nbtigpvl=0.;nbtalvl=0.;nbtaltl=0
		iflor=0.; stflor=0; niflor=0.; astdurjourp=12.
		!*lai
		sdjlai=0.;lai=0.;ei=0.
		!*stalk height
		htvd=0
		!*biomasses
		dmtbl=0.;gdmbl=0.;rdcan=0.; dmbaert=0
		dmstru=0;mst=0;dmbaer=0.;dmst=0;dmstm=0;dmsug=0;watstm=0
		dmstst=0;dmtbl=0;dmgbl=0;dmsen=0;dmrac=0;yldcan=0
		dmndftu=0;dmndfae=0;dmceltu=0;dmcelae=0;dmhceltu=0
		dmhcelae=0;dmligntu=0;dmlignae=0
		ndftu=0;ndfae=0;holndftu=0;holndfae=0
		lignndftu=0;lignndfae=0
		ndef1=1;ndef2=1
		ind0=0;ind1=0;ind2=0;ind3=0;ind4=0
		ind5=0;ind6=0;ind7=0;ind8=0;ind9=0;ind10=0;ind11=0;ind12=0
		sdjblt=0.;sdjhtvdt=0.;sswdeft=0.;sswdf1t=0.;sswdf2t=0.
	    stmot=0.;stnt=0.;stxt=0.; spart=0.;sparit=0.
		stmofm=0.;sparifm=0.;sprfm=0.;swdeffm=0.;swdf1fm=0.;swdf2fm=0
		stnfm=0.;stxfm=0
		cum1=0;cum2=0;cum3=0;cum4=0;cum5=0  !variables vrtuelles.

        istage=4
	  !4****** divers,général
		deproo=0 !root depth (cm)
		nage=0 !initialisation days counter
		jat=0; forc1eff=0
!********* END PART3N
c********************************************************************************




c********************************************************************************
c********* PART3O Open weather daily data file and search first line (beginnng date)
	IF (adjustl(meteostr)/="0_") THEN  ! JFM120126	
		i=1 
	        finmet(1)='non'; finplu(1)='non'   
			if (meteosta/='0') then	!
			s_err='ouverture fichier meteo impossible : '//meteosta !
			open(15,file=trim(repmeteo)//trim(meteosta)//".txt",  !
	1					status='old',err=99)  		
			read(15,*) s1; iom(1)=0; n5=0
			DO  !scrutation meteo file
				read (15,*,IOSTAT=iom(i)) datefich,amprecip(i)  
	1			,ameeq(i),amtempmx(i),amtempmn(i), amsolrad(i)
 			    call datesplit(datefich,mdat1,ddat1,ydat1)  
				jdatemet(i)=julday(mdat1,ddat1,ydat1) 
				if (jdatemet(i)==jdatebeg)  EXIT
				if (jdatemet(i)>jdatebeg) THEN
					if (nosim==1) then
					write(*,*) 'Debut trait < debut meteo '//meteosta 
					write(7,*) 'Debut trait < debut meteo '//meteosta 
					end if
					exit
				end if
				IF ((iom(i)<0)) THEN !Fin fichier météo
					if (nosim==1) then
					write(*,*) 'Debut trait > fin meteo '//meteosta !jfm20170524
					write(7,*) 'Debut trait > fin meteo '//meteosta !jfm20170524
					end if
					finmet(i)='oui'
					exit
				END IF
			END DO !scrutation meteo file
			end if !if (meteosta/='0') then
			IF (meteostr=='0') EXIT
	END IF    !	IF (adjustl(meteostr)/="0_") THEN  ! JFM120126
c********* END PART3O
c********************************************************************************





c********************************************************************************
C ********* PART3P Initialisation of soil water balance parameters **************
c vitenr = growth rate of root depth (cm/degré-jour)
c pDSTMX = root length density (cm rac. / cm3 sol)
c pBASET = base temperature for root depth (C)
c rwumx = maximum root uptake mm / cm de racines
c  gam = gamma coefficient  (not in ceres original: See artcle gabriel)
c swcon = CH a saturation  (not in ceres original: See artcle gabriel)
c aaa = A coefficient for CH  (not in ceres original: See gabriel paper)
	pDSTMX=5.0; pBASET=0 !Ceres de Gabriel  (valeurs initiales)
	rwumx=0.035  ! Van antwerpen  (0.03 dans Ceres de gabriel) 
	gam=0.019;swcon=4.5;aaa=82.5  !Ceres de gabriel
	straw=2000; sdep=0;scn=75;root=300;rcn=75;cntot=9
	vitenr=0.08  !valeur dans Ceres
	if (adjustl(typsol)=="latosol") then !JFM040613
	vitenr=0.03; pDSTMX=0.5
	end if
	seuilruis=80;pcruis=80
C ********* END PART3P
c********************************************************************************


c********************************************************************************
C********* PART3Q update of parameters values for sensitivity (CALL MAJ_VAL_PARAM) 
		IF (typeparam==1) Then
			DO ipar = 1 , nbpar
					par_val(ipar) = par_min(ipar) + 
     s				par_ecart(ipar) * (par_indice(ipar) - 1)  
			END DO
		END IF
		IF (typeparam>0) Then
			CALL MAJ_VAL_PARAM
		END IF
c********************************************************************************


c********************************************************************************
C ********* PART3R root profile initialisation  ******************
	!** limitation of profrac (root depth) down to soil depth
	dep=0	
	do i=1,nlayr
	dep=dep+dlayr(i)
	end do
	if (profrac>dep) profrac=dep
	if (parac==1) then !JFM080701  si calage ou sensibilité de profrac
	DO ipar = 1 , nbpar
	if (par_nom(ipar)=="profrac")then !JFM080701 prfrac<dep (profsol)
	  if (par_val(ipar)>dep)  par_val(ipar)=dep
	end if  
	END DO !ipar = 1 , nbpar
	END IF !(parac==1)
	pRDPMX=profrac
	!** initialisation of rlv()
	do i=1,nlayr
		RLV(i)=0
	end do
	!** initialisation deproo and rlv() when ratoons
	IF (ratoon>0) then
	  dep=0
		do i=1,nlayr
			RLV(i)=0
		end do
		do i=1,nlayr
			if (profrac>=(dep+dlayr(i))) then
	          v1=384649.* (1/(dep+30))**(3.4163)
	          v2=384649.* (1/(dep+dlayr(i)+30))**(3.4163)
			  RLV(i)=(v1+v2)/2
			 else
	          v1=384649.* (1/(dep+30))**(3.4163)
	          v2=384649.* (1/(profrac+30))**(3.4163)
			  RLV(i)=(v1+v2)/2			  
			end if
			dep=dep+dlayr(i)
			deproo=dep  
			if (profrac<=dep) then
			deproo=profrac !jfm051007
			exit
			end if 	 
		end do
	END IF
C ********* END PART3R
c********************************************************************************



	stopmod=0
C**********************************************************************************
C ************ PART4 Begining of DAILY LOOP ***************************************
	CONTINUE  
	DO !Begining of DAILY LOOP

c********************************************************************************
	!******** PART4A  Calculation of AGE, DATES and Test of daily loop ending
	nage=nage+1  ! Age de la cultue (jours)
	call joursuivant(dday,mday,yday) ! jour, mois et année du jour actuel
	call CALC_DECADE (mday,dday,decade)! decade: N° de décade du jour actuel
	call dat_str(dateact,dday,mday,yday)! date (chaine) du jour actuel
	jcal=julday(mday,dday,yday)-julday(1,1,yday)+1 !jour calendaire JFM130513	
	IF (nage>jdateend-jdatebeg) EXIT ! Sortie de la boucle journal. si date de fin atteinte
	daytoend=jdateend-julday(mday,dday,yday)+1 !Nombre de jours avant la fin
      call Astrono1(nage,jcal,lati,longi,alti,astlatrad,astdr,astdecl,
	1	astomeg,astrgo,astrgx,astdurjour,astdurjourp)
	IF ((mday*100+dday)==(mstopmod*100+dstopmod)) stopmod=1
c************** END PART4A	


c********************************************************************************
c************** PART4B  Reading of daily weather data
			  i=1
			  if (jdatemet(i)<(jdatebeg+nage)) then 
				IF ((iom(i)>=0)) THEN
				read (14+i,*,IOSTAT=iom(i)) datefich,
	1			amprecip(i),ameeq(i),amtempmx(i),amtempmn(i),
     2			amsolrad(i),amux(i),amun(i),amvt(i)	,amtempmm(i)
				call datesplit(datefich,mdat1,ddat1,ydat1)  
				mmet(i)=mdat1; dmet(i)=ddat1; ymet(i)=ydat1  
				jdatemet(i)=julday(mmet(i),dmet(i),ymet(i))
				end if
			  end if !(jdatemet(i)<(jdatedeb+nage))
			  IF ((iom(i)<0)) THEN  
				s1= 'date fin traitement > date fin meteo '//meteosta !
				if (finmet(i)=='non') write(7,*) trim(s1)
				jdatemet(i)=jdatebeg+nage+10; finmet(i)='oui'
			  END IF
			  if (jdatemet(i)==(jdatebeg+nage)) then
				precip=amprecip(1); eeq=ameeq(1); tempmx=amtempmx(1); 
				tempmn=amtempmn(1); solrad=amsolrad(1)
				tempmm=amtempmm(1)  
			  end if
c************** END PART4B	

	spr=spr+precip
	setp=setp+eeq

c********************************************************************************
c************** PART4C  Test of weather data
	if ((precip*eeq*solrad<0).OR.((tempmn+tempmx)/2<0)) then
	 if (precip<0) then;precip=0.;write(7,*) "RR<0: ",dateact;end if
	 if (eeq<0) then;eeq=0.5;write(7,*) "ETP<0: ",dateact;end if
	 if (solrad<0) then;solrad=0.5;write(7,*) "RG<0: ",dateact;end if
	 if (tempmn<0) then;tempmn=0.;write(7,*) "TN<0: ",dateact;end if
	 if (tempmx<0) then;tempmx=1.;write(7,*) "TX<0: ",dateact;end if
	 if (tempmx<tempmn) then
	 tempmx=tempmn+5.;write(7,*) "TX<TN: ",dateact;end if						        
	end if
c************** END PART4C


c********************************************************************************
c******** PART4D  Reading of irrigation data from table irr() 
c******                 and irrigation strategies calculations
	 irdose=0.
	 SELECT CASE (stirr1)
	 CASE (0)   !Pluvial
	 irdose=0
	 CASE (1)   ! Irrigations reelles / lecture fichier
	 if (nirr>0) then
	 DO l=1,nirr	
	  IF (irr(l).irrday==(jdatebeg+nage)) THEN
	   kindirr=irr(l).kindirr !JFM060624
	   irdose=irr(l).irrdose * irr(l).effirr !jfm040317
	       sirr=sirr + irdose		
	  END IF
	 END DO
	 end if !	 if (nirr>0) then
	 CASE (2)
	  i=stirr2  
	 if (modulo(nage,i)==0) then ! Corrigé et modifié le 18/08/2007 JFM070828 
		if (stockfr<((dulwr-llwr)*stirr3/100.)) then
		      irdose=((dulwr-llwr)*stirr3/100.)-stockfr
	    end if
		if (irdose<0.) irdose=0.
		if ((jdateend-(jdatebeg+nage))<stirr4) irdose=0.
	       sirr=sirr + irdose
	 end if
	CASE (3)  
	  !stirr2 :   % seuil minimum de déclenchement de l'arrosage
	  !stirr3 :   % seuil maximum de remplissage 
	  !stirr4 :   sevrage jours avant récolte
		if (stockfr<((dulwr-llwr)*stirr2/100.)) then
		      irdose=((dulwr-llwr)*stirr3/100.)-stockfr
	    end if
		if (irdose<0.) irdose=0.
		if ((jdateend-(jdatebeg+nage))<stirr4) irdose=0.
	       sirr=sirr + irdose
	    	
	END SELECT
c************** END PART4D


c********************************************************************************
C******** PART4E  Reading of nitrogen fertilisation data from table fertN()	
	DO l=1,nfertN		
	  IF (fertN(l).fertNday==(jdatebeg+nage)) THEN
			aNitr=aNitr + fertN(l).fertNdose 
	       dNitr=dNitr + fertN(l).dlayN
	  ENDIF
	END DO
C******** END PART4E


c********************************************************************************
!******* Forcing before calls (for ei forcing by observations)
	isobei=0	
	    !******** 1ere date
	IF ((nforc1==1>0) .AND. (forc1eff==0)) THEN
	  ! Est ce que l'on arrive à une date d'observation
	  IF ((obpjdate==(jdatebeg+nage)).AND.
	1      (iobs==1).AND.(adjustl(obptrt) == adjustl(trtname))) THEN
		!valeur de la variable à forcer
		select case (trim(adjustl(varobsim1)))
			case ('ei')
			   if(obei/=-999) then 
			     isobei=1; forc1eff=1
			   end if
		end select
	  END IF
	END IF
		    !******** toutes dates
	IF (nforct>0) THEN
	  ! Est ce que l'on arrive à une date d'observation
	  IF ((obpjdate==(jdatebeg+nage)).AND.
	1      (iobs==1).AND.(adjustl(obptrt) == adjustl(trtname))) THEN
		select case (trim(adjustl(varobsimt)))
			case ('ei')
			   if(obei/=-999) then 
			   isobei=1
			   end if
		end select
	  END IF
	END IF
c************* END FORCING

 
c********************************************************************************
c******** PART4F  Calls of growth and water balance routines
	IF (stopmod==0) then !JFM061207 stop calls
		do l=1,nbmod 
			s3=trim(adjustl(bilanh))
			call valcol(l,s3,"_")	
	        call appelmodul(s3)
		end do 
	  !* transformation of var. plante mosicas in var; plante ncsoil
		dmlf=dmgbl;	dmtops=dmbaer-dmgbl; dmroot=dmrac
	  gdmtop=gdmbaer-gdmbl+sdmbl; gdmroo=gdmrac
	 !call nuptak
	END IF !IF (stopmod==0) then
c************** END PART4F


c********************************************************************************
c************** PART4G  FORCING with observations of first date
	! A t'on forcage et at'il déjà été effectué
	IF ((nforc1>0) .AND. (forc1eff==0)) THEN
	  ! Est ce que l'on arrive à une date d'observation
	  IF ((obpjdate==(jdatebeg+nage)).AND.
	1      (iobs==1).AND.(adjustl(obptrt) == adjustl(trtname))) THEN
		!valeur de la variable à forcer
		select case (trim(adjustl(varobsim1)))
	        case ('stahtvd')
			   if(obhtvd/=-999.) then
					htvd=obhtvd; forc1eff=1
			   end if
	        case ('hv') 				!A FAIRE
			case ('glai')
			   if(oblai/=-999) then
			            lai=oblai; forc1eff=1
	           end if 
			case ('ei')
			   if(obei/=-999) then 
			   ei=obei; lai=-log(1-ei)/ke
			   end if
	        case ('agrdm')
			   if(obmsa/=-999) then
					dmbaer=obmsa*100; forc1eff=1  !JFM141111
			   end if
			case ('stamdm')
			   if(obmstigusi/=-999) then
					dmstm=obmstigusi*100; forc1eff=1 !JFM141111
			   end if
	        case ('blatnb')
			   if(obnbt/=-999) then
					blnbtot=obnbt; forc1eff=1
			   end if
	        case ('blagnb')
			   if(obnbv/=-999) then
					blnbv=obnbv; forc1eff=1
			   end if
	        case ('blastvd')
			   if(obltvdsurf/=-999) then
					bltvdsurf=obltvdsurf; forc1eff=1
			   end if
	        case ('stagnb')
			   if(obnbtigv/=-999) then
					nbtigv=obnbtigv; forc1eff=1
			   end if
	        case ('stamfm')
			   if(obrdcan/=-999) then
					yldcan=obrdcan; forc1eff=1
			   end if
	        case ('stamsu')
			   if((obsucretig/=-999).AND.(obrdcan/=-999)) then
					sug_stm=obsucretig*obrdcan; forc1eff=1
	           end if
	        case ('stock')	    ! A FAIRE
		end select
	  END IF
	END IF
c************** END PART4G


c********************************************************************************
c************** PART4H  FORCING with observations of all dates
	! A t'on forcage et at'il déjà été effectué
	!write(*,*) nforct, obpjdate-jdatebeg-nage
	IF (nforct>0) THEN
	  ! Est ce que l'on arrive à une date d'observation
	  IF ((obpjdate==(jdatebeg+nage)).AND.
	1      (iobs==1).AND.(adjustl(obptrt) == adjustl(trtname))) THEN
		!valeur de la variable à forcer
		write(*,*) htvd, obhtvd
		select case (trim(adjustl(varobsimt)))
	        case ('stahtvd')
			   if(obhtvd/=-999) htvd=obhtvd
	        case ('hv') !A FAIRE
			case ('glai')
			   if(oblai/=-999) lai=oblai
			case ('ei')
			   if(obei/=-999) then 
			   ei=obei; lai=-log(1-ei)/ke
			   end if
	        case ('agrdm')
			   if(obmsa/=-999) dmbaer=obmsa*100 !JFM141111
	        case ('stamdm')
			   if(obmstigusi/=-999) dmstm=obmstigusi*100 !JFM141111
	        case ('blatnb')
			   if(obnbt/=-999) blnbtot=obnbt
	        case ('blagnb')
			   if(obnbv/=-999) blnbv=obnbv
	        case ('blastvd')
			   if(obltvdsurf/=-999) bltvdsurf=obltvdsurf
	        case ('stagnb')
			   if(obnbtigv/=-999) nbtigv=obnbtigv
	        case ('stamfm')
			   if(obrdcan/=-999) yldcan=obrdcan
	        !case ('rich')
			!   if(obsucrec/=-999) htvd=obsucrec
	        case ('stamsu')
			   if((obsucretig/=-999).AND.(obrdcan/=-999))	then
			     sug_stm=obsucretig*obrdcan
	           end if
	        case ('stock')	    ! A FAIRE
			   if(obstock/=-999) stock=obstock

		end select
		write(*,*) htvd, obhtvd

	  END IF
	END IF
c************** END PART4H




c********************************************************************************
c******** PART4I  Writting observations in dnobs and Simulations results in simtemp 
	if (nobsol>0) then !calcul des stocks obs et sim pour calage et sortie foutgro
		ios=0
		DO l=1,nobsol
				IF (obsol(l).swoday==(jdatebeg+nage)) THEN
	        ios=l; exit
				end if
		end do
		if (ios>0) then
		  if (obsol(ios).swo(1)==-999.) then
			obstock=-999; obstockfr=-999
	      else
		    obstock=0.;v1=0
			do l=1,nlayr
				obstock=obstock+dlayr(l)*obsol(ios).swo(l)*10
		    end do
		    obstock=obstock-llw
		    dep=0; obstockfr=0
			do l=1,nlayr
			 if (profrac>dep+dlayr(l)) then
				obstockfr=obstockfr+dlayr(l)*obsol(ios).swo(l)*10
			  else
				obstockfr=obstockfr+dlayr(l)*obsol(ios).swo(l)*10
				obstockfr=obstockfr*(profrac-dep)/dlayr(l)
		     end if
			 dep=dep+dlayr(l)
			 if (profrac<=dep) exit 
		    end do

		  end if !if (obsol(ios)%swo(1)==-999.) then
		end if !if (ios>0) then
	end if
	do i=1,nvaropt  !JFM130117
	varparam=varopt(i)  !JFM130117
	If (typeparam>1) then   !JFM130118 typeparam>1
	 if((trim(varparam)=="stock").OR.(trim(varparam)=="stockfr")) then
			if (ios>0) then
			  if(trim(varparam)=="stock") then
				write(10,FMT='(F13.6)') obstock  ! écriture valeur dans dnobs
				write(12,FMT='(F13.6)') stock  ! écriture valeur dans Simtemp
	          end if
			  if(trim(varparam)=="stockfr") then
				write(10,FMT='(F13.6)') obstockfr  ! écriture valeur dans dnobs
				write(12,FMT='(F13.6)') stockfr  ! écriture valeur dans Simtemp
	          end if
			end if
		else
	 		IF ((obpjdate==(jdatebeg+nage)).AND.
	1          (adjustl(obptrt) == adjustl(trtname))) THEN
				call VALVARCAL(varparam,robs,rcalc)
				if (robs>-999) then
						!!calage pour certaines valeurs de la variable à caler robs<=v ou robs>=... 
	                    !if (robs<=6.2) then    ! exemple nage<220, lai:robs<=5
				  !write(10,FMT='(F13.6)') robs  ! écriture valeur dans dnobs
				  !write(12,FMT='(F13.6)') rcalc  ! écriture valeur dans dnobs
				  !write(11,FMT='(A15)') adjustl(varparam)  ! écriture valeur dans nvobs JFM130117
				  nobs=nobs+1; dnobst(nobs)=robs						!JFM130117
				  simtempt(nobs)=rcalc;nvobst(nobs)= adjustl(varparam)	!JFM130117			
						!write(*,*) robs, lai, rcalc
						!end if
				end if
			end if
	 end if	!(trim(varparam)=="stock")
	end if  !(typeparam==2)
	end do   !JFM130117
c******** END PART4I




c********************************************************************************	
c************** PART4J  Writting outputs in file***.csv
	!iobs: entrée de sim.txt 1 (sortie aux dates d'observ) 0 (pas de sortie
	!ios: =1 si observation sol existante ce jour et iobs=1 ou sortie frequentiel ce jour sinon 0
	c=';'
	valeur=setr/setm;v9=-999.0
		If (typeparam<2) then  !Si normal ou sensibilité
			ios=0
			IF ((nobsol>0).AND.
	1             ((iobs==1).OR.(modulo(nage,ifreq)==0))) then ! Il existe des observations sol pour ce traitement
				DO l=1,nobsol ! recherche dans le tableau obsol()
				IF (obsol(l).swoday==(jdatebeg+nage)) THEN  ! Une date d'observation existe
				ios=l; exit
				end if
				end do
			END if
	 	IF (((iobs==1).AND.((ios>0).OR.((obpjdate==(jdatebeg+nage))
	1      .AND.(adjustl(obptrt) == adjustl(trtname)))))  .OR.
     2      (modulo(nage,ifreq)==0) .OR. (nage==jdateend-jdatebeg)) THEN
		  !call wr_trtname(notrait,trtname,lati,longi,alti,c)
		  i=len(trim(trtname));  Write( s1, '(i20)' ) i   ! Converti i en chaine
	      s1="(I3,A1,A"//trim(adjustl(s1))//",A1,F12.6,A1,F12.6,A1,F6.1)"
	      write(6,FMT=trim(adjustl(s1)),ADVANCE="NO") notrait,c,
	1             trtname,c,lati,c,longi,c,alti
		  IF (typeparam==1) then  !si sensibilité
			do i=1,nbpar
			 write(6,FMT="(A1,F10.5)",ADVANCE="NO") c,par_val(i)
			end do
		  END IF
	      ! ***** Sorties variables plante simulées dans tous les cas !JFM141111 modif pour T/Ha ../100
	write(6,FMT="(A1,A10,2(A1,I4),2(A1,F9.6),20(A1,F9.3))",ADVANCE="NO") 
     1		c,dateact,c,nage,c,jat,c,lai,c,ei,c,blnbv,c,blnbtot,
	2		c,nbtigv,c,nbtigt,c,htvd,c,bltvdsurf,c,yldcan,c,dmstm/100,
	3		c,hum_stm,c,dmsug/100,c,sug_stm,c,v9,c,v9,c,v9,
	4		c,dmbaer/(100-hum_aer),c,dmbaer/100,c,hum_aer,
     5		c,v9,c,v9
		  ! ***** output observed plant data   !JFM150220  16775 - 193.87 * Humidité%
	 	  v1=obnbv+obnbm; if (v1==-1998) v1=-999
		  IF ((obpjdate==(jdatebeg+nage)) !observation plante ce jour
	1          .AND.(adjustl(obptrt) == adjustl(trtname))) THEN
			write(6,FMT="(23(A1,F9.3))",ADVANCE="NO")
	1	    c,obnbtigv,c,obnbtigm,c,obnbtigfl,c,obhtvd,c,
     2		obnbt,c,obnbv,c,obltvdsurf,c,oblai,c,obei,c,obsucre,c,
     3        obrdcan,c,obmstigusi,c,obhumtig,c,obsucretig,c,
     4        obfibretig,c,obbrixjtig,c,puretejtig,c,obmslv,c,obmslt,c,
     5        obmfa,c,obmsa,c,obhumcan,c,obmsen
			ELSE ! if no plant observations this day
			write(6,FMT="(23(A1,F9.3))",ADVANCE="NO") 
     1			c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,
	2			c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,
     3			c,v9,c,v9,c,v9
		  END IF 
		  ! ***** output observed soil data
	      IF (ios>0) then !observation sol ce jour (ios>=1)
			write(6,FMT="(17(A1,F9.3))",ADVANCE="NO") 
	1			c,sw(1),c,sw(2),c,sw(3),c,sw(4),c,sw(5),c,sw(6),
	2			c,stock,c,stocksurface,c,stockfr,c,obsol(ios).swo(1),
     3			c,obsol(ios).swo(2),c,obsol(ios).swo(3),
	4			c,obsol(ios).swo(4),c,obsol(ios).swo(5),	    
     5			c,obsol(ios).swo(6),c,obstock,c,obstockfr   	 
			ELSE ! pas d'observation sol de jour (ios=0)
			write(6,FMT="(17(A1,F9.3))",ADVANCE="NO") 
	1			c,sw(1),c,sw(2),c,sw(3),c,sw(4),c,sw(5),c,sw(6),
	2			c,stock,c,stocksurface,c,stockfr,
	3			c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9,c,v9	 
		  end if  
		  ! ***** output observed root data
		  ! ***** Sortie variables d'état racinaires simulées Sol
		  write(6,FMT="(7(A1,F9.3))",ADVANCE="NO")
	1	  c,deproo,c,rlv(1),c,rlv(2),c,rlv(3),c,rlv(4),c,rlv(5),c,rlv(6)
		  ! ***** output bio-weather index
		  write(6,FMT="(51(A1,F9.3))")
     1	  c,precip,c,tempmn,c,tempmx,c,tempmm,c,solrad,c,eeq,c,etm,
	2	  c,etr,c,tmp,c,trp,c,swdef,c,swdf1,c,swdf2,c,irdose,c,bldj,
	3	  c,htvddj,c,pari,c,stmo,c,spar,c,spari,c,spr,c,sirr,c,setp,
	4	  c,setm,c,setr,c,stmp,c,strp,c,blsdj,c,sswdef/nage,
	5	  c,sswdf1/nage,c,sswdf2/nage,c,htvdsdj,c,sdrain,c,srunoff,
	6	  c,sdjblt,c,sdjhtvdt,c,sparit,c,sswdeft,c,sswdf1t,c,sswdf2t,
	7	  c,stmot,c,setmt,c,setrt,c,srunofft,c,sdraint,c,setpt,c,stmpt,
	8	  c,talsdj,c,astdurjour,c,srgx  !strpt
		ENDIF  
	end if  !(typeparam==2)
c************** END PART4J





c********************************************************************************					
c******** PART4L  Reading observation data of following day in observation table/file
		IF ((ioplt<0) .OR. (adjustl(obptrt) /= adjustl(trtname))) THEN !Last modif JFM150220
	    ELSE
		  IF (obpjdate==(jdatebeg+nage)) THEN
		   DO	
		read(4,"(A25,2(A1,I2),A1,I4,23(A1,F9.3),A1,I1)",IOSTAT=ioplt) 
	1	obptrt,c,obpdd,c,obpmm,c,obpyr,c,obnbtigv,c,onnbtigm,c,
	2    obnbtigfl,c,obhtvd,c,obnbt,c,obnbv,c,obltvdsurf,c,
     3    oblai,c,obei,c,obsucre,c,obrdcan,c,obmstigusi,c,obhumtig,c,
     4    obsucretig,c,obfibretig,c,obbrixjtig,c,puretejtig,c,obmslv,c,
     5    obmslt,c,obmfa,c,obmsa,c,obhumcan,c,obmsen,c,obsutil
			 IF (ioplt<0) EXIT
			 IF (adjustl(obptrt) == adjustl(trtname)) THEN
     				obpjdate=julday(obpmm,obpdd,obpyr)
				EXIT	        
			 END IF
			END DO
		  END IF	
		END IF
c************** END PART4L


	END DO !début boucle journalière	
C ************ END DAILY LOOP *******************************************
C**********************************************************************************


	if (meteosta/='0') close(15)
	ii=julday(1,1,iyrend)+jdateend-julday(1,1,iyrbeg)-jdatebeg
	write(7,*)'nage=',nage,ii
	close (4) ! Fichier observ plante
	!stop

	
	END DO  !do notrait=1,nbtrt
c ************* END TREATMENT LOOP *********************************************
C**********************************************************************************


      !**Verification e,d paramétrage and output
	  CALL VER_FIN_PARAM

	close(7)

	  
      END DO  !Do While (finparam == "NON")
c ************* END SENSITIVITY LOOP *********************************************
C**********************************************************************************



	  !fichier simulation sim.dat
	 if (testsub=='oui') close(30)
	 close(6)
	 close(9);	!close(10);	close(12); close(11)
	if (perso=="oui") close(14)

	return
99	write(7,*) s_err; write(*,*) s_err
	stop
	
      end Subroutine motsimul
C**** end simulating engine
C******************************************************************************
C******************************************************************************
C******************************************************************************	

	SUBROUTINE LITPLANTE
	! Lit le fichier plante.txt
        INCLUDE 'cercane1.inc' !JFM000301
	    character*80 s1
	    integer io
	    
		io=0
		OPEN(3,File=trim(repdata)//'plante.txt')  !jfm000315
			REWIND(3); READ(3,*) s1
	!       !write(*,*) "n",cultivar
			DO  !scrutation fichier plante.txt, remplissage avec R570
				READ (3,*,IOSTAT=io) var,param,v1,v2	
				IF (adjustl(var)=='R570') THEN 
					val=v2; IF  (ratoon==0) val=v1
					call initplante
				END IF
				var=''	
				IF ((io<0)) EXIT
			END DO !scrutation fichier plante.txt, remplissage avec cultivar exact
			io=0;REWIND(3); READ(3,*) s1
			DO  !scrutation fichier plante.txt
				READ (3,*,IOSTAT=io) var,param,v1,v2
	            
				IF (adjustl(trim(var))==adjustl(trim(cultivar))) THEN
					val=v2; IF  (ratoon==0) val=v1
	!write(*,*) var,param,v1
					call initplante
				END IF
				var=''	
				IF ((io<0)) EXIT
			END DO !scrutation fichier plante.txt
		CLOSE(3)
		!write(*,*) blsurfk,blsurfmax
	END SUBROUTINE LITPLANTE

      SUBROUTINE INITPLANTE
          INCLUDE 'cercane1.inc' !JFM000301
           SELECT CASE (adjustl(param))             
                Case ("aa1")
                    aa1 = val !;Print *,"aa1",aa1
                Case ("aa2")
                    aa2 = val !;Print *,"aa2",aa2
                Case ("aa3")
                    aa3 = val !;Print *,"aa3",aa3
                Case ("aa4")
                    aa4 = val !;Print *,"aa4",aa4
                Case ("aa5")
                    aa5 = val !;Print *,"aa4",aa4
                Case ("blk")
                    blk = val !Print *,"bhumtig", bhumtig
                Case ("blkdec")
                    blkdec = val !Print *,"bhumtig", bhumtig
                Case ("blsurfk")
                    blsurfk = val !Print *,"bhumtig", bhumtig
                Case ("blsurfmax")
                    blsurfmax = val !Print *,"bhumtig", bhumtig 
                Case ("bltb")
                    bltb = val !Print *,"bhumtig", bhumtig
                Case ("blvmax")
                    blvmax = val !Print *,"bhumtig", bhumtig
			  Case ("con_va")
                    conv_a = val !;Print *,"conv_a",conv_a
                Case ("conv_b")
                    conv_b = val !;Print *, "conv_b",conv_b
                Case ("deb_tigus")
                    deb_tigus = val !;Print *,"deb_tigus",deb_tigus
                Case ("debei")
                    debei = val !;Print *,"debei", debei
                Case ("debmsa")
                    debmsa = val !;Print *,"debmsa", debmsa
                Case ("dsenftig")
                    dsenftig = val !;Print *,"dsenftig", dsenftig
                Case ("ebmax")
                    ebmax = val !;Print *,"ebmax", ebmax
                Case ("ec")
                    ec = val !;Print *,"ec", ec
                Case ("eideb")
                    eideb = val !;Print *,"eideb", eideb
                Case ("eimax")
                    eimax = val !;Print *,"eimax", eimax
                Case ("epb")
                    epb = val !;Print *,"epb", epb
			  Case ("htvdcan")
                    htvdcan = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdgro1")
                    htvdgro1 = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdgro2")
                    htvdgro2 = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdtb")
                    htvdtb = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdtm")
                    htvdtm = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdto")
                    htvdto = val !;Print *,"htvdrdeb", htvdrdeb
			  Case ("humtigk1")
                    humtigk1 = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("humtigk2")
                    humtigk2 = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerdeb")
                    humaerdeb = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerfin")
                    humaerfin = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaertb")
                    humaertb = val !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerdec")
                    humaerdec = val !;Print *,"htvdrdeb", htvdrdeb
                !Case ("kintercplai")
                !    kintercplai = val !;Print *,"kintercplai", kintercplai
                Case ("kcmax")
                    kcmax = val !;Print *,"kcmax", kcmax
                Case ("ke")
                    ke = val !;Print *,"ke", ke
                Case ("keidj")
                    keidj = val !;Print *,"keidj", keidj
                Case ("kexrn")
                    kexrn = val !;Print *,"kexrn", kexrn
                Case ("klaidj")
                    klaidj = val !;Print *,"klaidj", klaidj
                Case ("ksenftig")
                    ksenftig = val !;Print *,"ksenftig", ksenftig
                Case ("ksenlai")
                    ksenlai = val !;Print *,"ksenlai", ksenlai
                Case ("ksurf")
                    ksurf = val !;Print *,"ksurf", ksurf 
                 Case ("laitb")
                    laitb = val !;Print *,"laideb", laideb
                Case ("laiwksen")
                    laiwksen = val !;Print *,"laideb", laideb
                Case ("laicroi")
                    laicroi = val !;Print *,"laideb", laideb
                Case ("maxsurf")
                    maxsurf = val !;Print *,"maxsurf", maxsurf
                Case ("msadeb")
                    msadeb = val !;Print *,"msadeb", msadeb
			  Case ("p01")
                    p01 = val !;Print *,"ph", ph
			  Case ("p02")
                    p02 = val !;Print *,"ph", ph
			  Case ("p03")
                    p03 = val !;Print *,"ph", ph
			  Case ("p04")
                    p04 = val !;Print *,"ph", ph
			  Case ("p05")
                    p05 = val !;Print *,"ph", ph
			  Case ("p06")
                    p06 = val !;Print *,"ph", ph
			  Case ("p07")
                    p07 = val !;Print *,"ph", ph
			  Case ("p08")
                    p08 = val !;Print *,"ph", ph
			  Case ("p09")
                    p09 = val !;Print *,"ph", ph
			  Case ("p10")
                    p10 = val !;Print *,"ph", ph
			  Case ("pbaset")
                    pbaset = val !;Print *,"ph", ph
			  Case ("ph")
                    ph = val !;Print *,"ph", ph
                Case ("pblfin")
                    pblfin = val !;Print *,"ph", ph
                Case ("pbldif")
                    pbldif = val !;Print *,"ph", ph
                Case ("pblkdec")
                    pblkdec = val !;Print *,"ph", ph
                Case ("pdstmx")
                    pdstmx = val !;Print *,"ph", ph
			  Case ("pflor")
                    pflor = val !;Print *,"ph", ph
                Case ("rwumx")
                    rwumx = val !;Print *,"ph", ph
                Case ("slafin")
                    slafin = val !;Print *,"ph", ph
                Case ("sladif")
                    sladif = val !;Print *,"ph", ph
                Case ("slakdec")
                    slakdec = val !;Print *,"ph", ph
                Case ("phyll")
                    phyll = val !;Print *,"phyll", phyll
                Case ("phyltbase")
                    phyltbase = val !;Print*,"phyltbase",phyltbase
                Case ("plai")
                    plai = val !;Print *,"plai", plai
                Case ("profrac_var")
                    profracvar = val*100 !;Print*,"profracvar",profracvar
                Case ("prue")
                    prue = val !;Print *,"prue ", prue
                Case ("psilb")
                    psilb = val !;Print *,"prue ", prue
                Case ("pstrutb")
                    pstrutb = val  !;Print *,"prue ", prue
                Case ("pstrutcroi")
                    pstrutcroi = val  !;Print *,"prue ", prue
                Case ("pstrufin")
                    pstrufin = val  !;Print *,"prue ", prue
                Case ("pstrudec")
                    pstrudec = val  !;Print *,"prue ", prue
                Case ("ptig")
                    ptig = val !;Print *,"kracdeb", kracdeb
                Case ("ptigdeb")
                    ptigdeb = val !;Print *,"kracdeb", kracdeb
                Case ("ptigdec")
                    ptigdec = val !;Print *,"kracdec", kracdec
                Case ("ptigfin")
                    ptigfin = val !;Print *,"kracfin", kracfin
                Case ("ruemax")
                    ruemax = val !;Print *,"ruemax", ruemax
                Case ("ruetk")
                    ruetk = val !;Print *,"ktconv", ktconv
                Case ("ruetopt")
                    ruetopt = val !;Print *,"ruemax", ruemax
                Case ("sthydcroi")
                    sthydcroi = val !;Print *,"ruemax", ruemax
                Case ("sthydbio")
                    sthydbio = val !;Print *,"ruemax", ruemax
			  Case ("taldebtt")
                    taldebtt = val !;Print *,"talfin", talfin
                Case ("talfinval")
                    talfinval = val !;Print *,"talrdeb", talrdeb
                Case ("talpeaktt")
                    talpeaktt=val !;Print*,"talrdebsen",talrdebsen
                Case ("talpeakval")
                    talpeakval = val !;Print *,"talrk", talrk
                Case ("taltb")
                    taltb = val !;Print *,"talrksen", talrksen
!			  Case ("tal2fin")
!                   tal2fin = val !;Print *,"talfin", talfin
!              Case ("tal2deb")
!                 tal2deb = val !;Print *,"talrdeb", talrdeb
!            Case ("tal2debsen")
!     !               tal2debsen=val !;Print*,"talrdebsen",talrdebsen
!                Case ("tal2kapp")
!                    tal2kapp = val !;Print *,"talrk", talrk
!                Case ("tal2ksen")
!                    tal2ksen = val !;Print *,"talrksen", talrksen
!  			  Case ("talfin")
!                    talfin = val !;Print *,"talfin", talfin
!                Case ("taldeb")
!                    taldeb = val !;Print *,"talrdeb", talrdeb
!                Case ("taldebsen")
!                    taldebsen=val !;Print*,"talrdebsen",talrdebsen
!                Case ("talkapp")
!                    talkapp = val !;Print *,"talrk", talrk
!                Case ("talksen")
!                    talksen = val !;Print *,"talrksen", talrksen
!                Case ("taltb")
!                    taltb = val !;Print *,"talrdeb", talrdeb  !initplant
!                Case ("taltm")
!                    taltm = val !;Print *,"talrdeb", talrdeb
!                Case ("talto1")
!                    talto1 = val !;Print *,"talrdeb", talrdeb
!                Case ("talto2")
!                    talto2 = val !;Print *,"talrdeb", talrdeb
!                Case ("taltp")
!                    taltp = val !;Print *,"talrdeb", talrdeb
                Case ("tbase")
                    tbase = val !;Print *,"tbase", tbase
                Case ("tbasehtvd")
                    tbasehtvd = val !;Print*,"tbasehtvd",tbasehtvd
                Case ("tbaselai")
                    tbaselai = val !;Print *,"tbaselai", tbaselai
                Case ("trescro")
                    trescro = val !;Print *,"trescro", trescro
                Case ("tresent")
                    tresent = val !;Print *,"tresent", tresent
                Case ("vitenr")
                    vitenr = val  !;Print *,"tresent", tresent
			!** Paramètres CLC Damien
			  Case ("tuaera")
                    tuaera = val !;Print *,"tuaera", tuaera
			  Case ("tuaerb")
                    tuaerb = val !;Print *,"tuaerb", tuaerb
			  Case ("tuaerk")
                    tuaerk = val !;Print *,"tuaerk", tuaerk
			  Case ("ndftua")
                    ndftua = val !;Print *,"ndftua", ndftua
			  Case ("ndftub")
                    ndftub = val !;Print *,"ndftub", ndftub
			  Case ("ndftuk")
                    ndftuk = val !;Print *,"ndftuk", ndftuk
			  Case ("celtua")
                    celtua = val !;Print *,"celtua", celtua
			  Case ("celtub")
                    celtub = val !;Print *,"celtub", celtub
			  Case ("celtuk")
                    celtuk = val !;Print *,"celtuk", celtuk
			  Case ("hceltua")
                    hceltua = val !;Print *,"hceltua", hceltua
			  Case ("hceltub")
                    hceltub = val !;Print *,"hceltub", hceltub
			  Case ("hceltuk")
                    hceltuk = val !;Print *,"hceltuk", hceltuk
			  Case ("ligntua")
                    ligntua = val !;Print *,"ligntua", ligntua
			  Case ("ligntub")
                    ligntub = val !;Print *,"ligntub", ligntub
			  Case ("ligntuk")
                    ligntuk = val !;Print *,"ligntuk", ligntuk
	      END SELECT
      END SUBROUTINE INITPLANTE






	SUBROUTINE datesplit(xdate,mm,id,iyyy)
	! converti un date chaine(*10) en 3 entiers id(jou),   jfm19990215
        CHARACTER xdate*10
	  INTEGER mm,id,iyyy
        mm=(ichar(xdate(4:4))-48)*10+(ichar(xdate(5:5))-48)
	  id=(ichar(xdate(1:1))-48)*10+(ichar(xdate(2:2))-48)	  
	  iyyy=(ichar(xdate(7:7))-48)*1000+(ichar(xdate(8:8))-48)*100
	  iyyy=iyyy+(ichar(xdate(9:9))-48)*10+(ichar(xdate(10:10))-48)	
	END SUBROUTINE datesplit

	SUBROUTINE dat_str(xdate,xjour,xmois,xan) !jfm19990215
	! fonction renvoie sous forme de texte (DATE_STR*10) une date spécfiée
	! par le jour(xjour), le mois (xmois) et l'année (xan) qui sont des entiers
	  INTEGER xjour,xmois,xan
        CHARACTER xdate*10
		xdate(3:3)='/'; xdate(6:6)='/'
	  xdate(1:1)=achar(xjour/10+48)
        xdate(2:2)=achar(xjour-(xjour/10)*10+48)
	  xdate(4:4)=achar(xmois/10+48)
        xdate(5:5)=achar(xmois-(xmois/10)*10+48)
	  xdate(7:7)=achar(xan/1000+48)
	  xdate(8:8)=achar((xan-(xan/1000)*1000)/100+48)
	  xdate(9:9)=achar((xan-(xan/100)*100)/10+48)
	  xdate(10:10)=achar((xan-(xan/10)*10)+48)
	END SUBROUTINE dat_str


	FUNCTION DATE_STR(xjour,xmois,xan) !jfm19990215
	! fonction renvoie sous forme de texte (DATE_STR*10) une date spécfiée
	! par le jour(xjour), le mois (xmois) et l'année (xan) qui sont des entiers
	  INTEGER xjour,xmois,xan
        CHARACTER DATE_STR*10
		DATE_STR(3:3)='/'; DATE_STR(6:6)='/'
	  DATE_STR(1:1)=achar(xjour/10+48)
        DATE_STR(2:2)=achar(xjour-(xjour/10)*10+48)
	  DATE_STR(4:4)=achar(xmois/10+48)
        DATE_STR(5:5)=achar(xmois-(xmois/10)*10+48)
	  DATE_STR(7:7)=achar(xan/1000+48)
	  DATE_STR(8:8)=achar((xan-(xan/1000)*1000)/100+48)
	  DATE_STR(9:9)=achar((xan-(xan/100)*100)/10+48)
	  DATE_STR(10:10)=achar((xan-(xan/10)*10)+48)
	END FUNCTION DATE_STR

	!********************************************** DEBUT CALC_DECADE **********
      ! convertit une date ( chaine xdate*10) en la decade correspondate
	!attentin ne fonctionne pas
	SUBROUTINE CALC_DECADE(mois,jour,dec) !jfm19990215
	  INTEGER dec,mois,jour,jj
        jj = Abs(Int((jour - 1) / 10)) + 1
          If (jj > 3) Then 
             jj = 3
	    END IF
          dec = (mois - 1) * 3 + jj
      END SUBROUTINE CALC_DECADE



	FUNCTION julday(mm,id,iyyy) !jfm19990215
	! calcule le jour julien du jour id du mois im de l'année iyyy
        INTEGER julday,id,iyyy,mm,IGREG
        PARAMETER (IGREG=15+31*(10+12*1582))
        INTEGER ja,jm,jy
        jy=iyyy
        if (jy.eq.0) pause 'julday: there is no year zero'
        if (jy.lt.0) jy=jy+1
        if (mm.gt.2) then
            jm=mm+1
          else
            jy=jy-1
            jm=mm+13
        endif
        julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
        if (id+31*(mm+12*iyyy).ge.IGREG) then
          ja=int(0.01*jy)
          julday=julday+2-ja+int(0.25*ja)
        endif
        return
      END  FUNCTION julday

	SUBROUTINE datcal(xdate,xyr,xjcal)
	! calcule le jour calendaire xjcal et l'année xyr pour une date xdate (__/__/____)
		character xdate*8
		integer xyr,xmm,xdd,xjcal
		call datesplit(xdate,xyr,xmm,xdd)
		xjcal=julday(xmm,xdd,xyr)-julday(1,1,xyr)+1
	END SUBROUTINE datcal
	
	
	SUBROUTINE INITPARCAL(nompar,valpar)
	INCLUDE 'cercane1.inc'
	character nompar*15
	real valpar
		select case (adjustl(nompar))
			  Case ("aa1")
                    aa1 = valpar  !Print *,"aa1",aa1
                Case ("aa2")
                    aa2 = valpar  !Print *,"aa2",aa2
                Case ("aa3")
                    aa3 = valpar  !Print *,"aa3",aa3
                Case ("aa4")
                    aa4 = valpar  !Print *,"aa4",aa4
                Case ("aa5")
                    aa5 = valpar  !Print *,"aa4",aa4
                Case ("blk")
                    blk = valpar  !Print *,"aa4",aa4
                Case ("blkdec")
                    blkdec = valpar  !Print *,"aa4",aa4
                Case ("blsurfk")
                    blsurfk = valpar  !Print *,"aa4",aa4
                Case ("blsurfmax")
                    blsurfmax = valpar  !Print *,"aa4",aa4
                Case ("bltb") 
                    bltb = valpar  !Print *,"aa4",aa4
                Case ("blvmax")
                    blvmax = valpar !Print *,"bhumtig", bhumtig
 			  Case ("ke")
                    ke = valpar  !;Print *,"klaidj", klaidj
                Case ("kcmax")
                    kcmax = valpar  !;Print *,"klaidj", klaidj
                Case ("kexrn")
                    kexrn = valpar  !;Print *,"klaidj", klaidj
			  Case ("htvdcan")
                    htvdcan = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdgro1")
                    htvdgro1 = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdgro2")
                    htvdgro2 = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdtb")
                    htvdtb = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdtm")
                    htvdtm = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("htvdto")
                    htvdto = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humtigk1")
                    humtigk1 = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humtigk2")
                    humtigk2 = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerdeb")
                    humaerdeb = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerfin")
                    humaerfin = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaertb")
                    humaertb = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("humaerdec")
                    humaerdec = valpar !;Print *,"htvdrdeb", htvdrdeb
                Case ("laitb")
                    laitb = valpar !;Print *,"laideb", laideb
                Case ("laiwksen")
                    laiwksen = valpar !;Print *,"laideb", laideb
                Case ("laicroi")
                    laicroi = valpar !;Print *,"laideb", laideb
                Case ("p0")
                    p0 = valpar !;Print *,"p0", p0
                Case ("p01")
                    p01 = valpar !;Print *,"p0", p0
                Case ("p02")
                    p02 = valpar !;Print *,"p0", p0
                Case ("p03")
                    p03 = valpar !;Print *,"p0", p0
                Case ("p04")
                    p04 = valpar !;Print *,"p0", p0
                Case ("p05")
                    p05 = valpar !;Print *,"p0", p0
                Case ("p06")
                    p06 = valpar !;Print *,"p0", p0
                Case ("p07")
                    p07 = valpar !;Print *,"p0", p0
                Case ("p08")
                    p08 = valpar !;Print *,"p0", p0
                Case ("p09")
                    p09 = valpar !;Print *,"p0", p0
                Case ("p10")
                    p10 = valpar !;Print *,"p0", p0
                Case ("pbaset")
                    pbaset = valpar !;Print *,"laideb", laideb
                Case ("pblfin")
                    pblfin = valpar !;Print *,"ph", ph
                Case ("pbldif")
                    pbldif = valpar !;Print *,"ph", ph
                Case ("pblkdec")
                    pblkdec = valpar !;Print *,"ph", ph
                Case ("pdstmx")
                    pdstmx = valpar !;Print *,"ph", ph
                Case ("pflor")
                    pflor = valpar !;Print *,"ph", ph
                Case ("profrac")
                    profrac = valpar !;Print *,"profrac", profrac
                Case ("profracvar")
                    profracvar = valpar !;Print *,"profrac", profrac

                Case ("psilb")
                    psilb = valpar !;Print *,"ph", ph
                Case ("slafin")
                    slafin = valpar !;Print *,"ph", ph
                Case ("sladif")
                    sladif = valpar !;Print *,"ph", ph
                Case ("slakdec")
                    slakdec = valpar !;Print *,"ph", ph
                 Case ("prue")
                    prue = valpar  !;Print *,"prue ", prue
                Case ("pstrutb")
                    pstrutb = valpar  !;Print *,"prue ", prue
                Case ("pstrutcroi")
                    pstrutcroi = valpar  !;Print *,"prue ", prue
                Case ("pstrufin")
                    pstrufin = valpar  !;Print *,"prue ", prue
                Case ("pstrudec")
                    pstrudec = valpar  !;Print *,"prue ", prue
			  Case ("ptig")
			      ptig = valpar !;Print *,"kracdeb", kracdeb                
			  Case ("ptigdeb")
                    ptigdeb = valpar !;Print *,"kracdeb", kracdeb
                Case ("ptigdec")
                    ptigdec = valpar !;Print *,"kracdec", kracdec
                Case ("ptigfin")
                    ptigfin = valpar !;Print *,"kracfin", kracfin
                Case ("ruemax")
                    ruemax = valpar  !;Print *,"ruemax", ruemax
                Case ("ruetk")
                    ruetk = valpar  !;Print *,"ruemax", ruemax
                Case ("ruetopt")
                    ruetopt = valpar  !;Print *,"ruemax", ruemax
                Case ("rwumx")
                    rwumx = valpar  !;Print *,"ruemax", ruemax
                Case ("sthydcroi")
                    sthydcroi = valpar !;Print *,"ruemax", ruemax
                Case ("sthydbio")
                    sthydbio = valpar !;Print *,"ruemax", ruemax

			  Case ("taldebtt")
                    taldebtt = valpar !;Print *,"talfin", talfin
                Case ("talfinval")
                    talfinval = valpar !;Print *,"talrdeb", talrdeb
                Case ("talpeaktt")
                    talpeaktt=valpar !;Print*,"talrdebsen",talrdebsen
                Case ("talpeakval")
                    talpeakval = valpar !;Print *,"talrk", talrk
                Case ("taltb")
                    taltb = valpar !;Print *,"talrksen", talrksen


                Case ("tbase")
                    tbase = valpar  !;Print *,"tbase", tbase
                Case ("vitenr")
                    vitenr = valpar  !;Print *,"tbase", tbase

			  !******** CLC
                Case ("tuaera")
                    tuaera = valpar  !;Print *,"tbase", tbase
                Case ("tuaerb")
                    tuaerb = valpar  !;Print *,"tbase", tbase
                Case ("tuaerk")
                    tuaerk = valpar  !;Print *,"tbase", tbase
                Case ("ndftua")
                    ndftua = valpar  !;Print *,"tbase", tbase
                Case ("ndftub")
                    ndftub = valpar  !;Print *,"tbase", tbase
                Case ("ndftuk")
                    ndftuk = valpar  !;Print *,"tbase", tbase
                Case ("celtua")
                    celtua = valpar  !;Print *,"tbase", tbase
                Case ("celtub")
                    celtub = valpar  !;Print *,"tbase", tbase                
			  Case ("celtuk")
                    celtuk = valpar  !;Print *,"tbase", tbase
                Case ("hceltua")
                    hceltua = valpar  !;Print *,"tbase", tbase
                Case ("hceltub")
                    hceltub = valpar  !;Print *,"tbase", tbase                
			  Case ("hceltuk")
                    hceltuk = valpar  !;Print *,"tbase", tbase
			  Case ("ligntua")
                    ligntua = valpar  !;Print *,"tbase", tbase
			  Case ("ligntub")
                    ligntub = valpar  !;Print *,"tbase", tbase
			  Case ("ligntuk")
                    ligntuk = valpar  !;Print *,"tbase", tbase
!	real tuaera,tuaerb,tuaerk, ndftua,ndftub,ndftuk
!	real celtua,celtub,celtuk, hceltua,hceltub,hceltuk
!	real ligntua,ligntub,ligntuk
		end select
	 END SUBROUTINE INITPARCAL

	SUBROUTINE VALVARCAL(nomvar,valobs,valcal)
	! MAJ des valeurs observées et simulées de la variable de calage
		INCLUDE 'cercane1.inc'
		character nomvar*12
		real valobs,valcal
	    valobs=-999
		select case (adjustl(nomvar))
				Case ('ei')
                    valobs = obei  !Print *,"aa1",aa1
                    valcal = ei  !Print *,"aa1",aa1
				Case ('glai')
                    valobs = oblai  !Print *,"aa1",aa1
                    valcal = lai  !Print *,"aa1",aa1
				Case ('stagnb')
                    valobs = obnbtigv  !Print *,"aa1",aa1
                    valcal = nbtigv  !Print *,"aa1",aa1
				Case ('blatnb')
	              !if(obnbt/=-999)	then
                    valobs = obnbt  !write (*,*) valobs
                    valcal = blnbtot  !Print *,"aa1",aa1	
				  !end if
				Case ('blagnb')
                    valobs = obnbv  !write (*,*) valobs
                    valcal = blnbv  !Print *,"aa1",aa1	
				Case ('blastvd')
                    valobs = obltvdsurf  !write (*,*) valobs
                    valcal = bltvdsurf  !Print *,"aa1",aa1	
				Case ("agrdm")
                    valobs = obmsa  !Print *,"aa1",aa1 !JFM141111
				  valcal = dmbaer/100
				Case ("agrwc")
                    valobs = obhumcan  !Print *,"aa1",aa1 !JFM141111
				  valcal = hum_aer
				Case ("blagdm")
                    valobs = obmslv  !Print *,"aa1",aa1
								valcal=dmgbl/100
				!Case ("dmstru")
                  !  valobs = (obmsa-obsucre)  !Print *,"aa1",aa1 !JFM141111
				!				valcal=dmstru/100
				  !if ((obmsa<=-0.9).OR.(obsucre<=-0.9)) valobs=-1.0   !avant JFM141111
			!if ((obmsa*100<=-0.9).OR.(obsucre*100<=-0.9)) valobs=-1.0   !JFM141111
				Case ("stamsu")
                    valobs = obsucre  !Print *,"aa1",aa1 !JFM141111
				  valcal=dmsug/100				
				Case ("stamfm")
                    valobs = obrdcan  !Print *,"aa1",aa1
								valcal = yldcan
				Case ("stahtvd")
                    valobs = obhtvd  !Print *,"aa1",aa1
								valcal = htvd
				Case ("stamdm")
                    valobs = obmstigusi  !Print *,"aa1",aa1  !JFM141111
								valcal = dmstm/100
				Case ("sucrecms")
	              if((obsucretig/=-999).AND.(obrdcan/=-999))	then
					if (obmstigusi>0.AND.(dmstm>0)) then
						valobs = obsucretig*obrdcan/obmstigusi  !Print *,"aa1",aa1
						valcal=dmsug/dmstm
					end if
				  end if						
				Case ("stamsuc")
                    valobs = obsucretig  !Print *,"aa1",aa1
				  valcal= sug_stm				
								!valcal=0
								!if (yldcan>0) valcal = dmsug/yldcan
				!Case ("watstm")
				!				valobs = obrdcan*100-obmstigusi
				!				valcal = watstm
				!				if ((obrdcan<=-0.9).OR.(obmstigusi<=-0.9)) 
	1			!									valobs=-1.0
     				Case ("stock")
								valobs = obstock
								valcal = stockf
		end select
	End subroutine



	SUBROUTINE VERIFSIMDAT(irotation,iobs,
	1			iwatbal,inbal,typeparam,nforc1,nforct)
	! routine de verification des données de sim.dat
	!character bilanh*20, foutgro*10
	!integer irotation,typeparam,typecalage,nforc1,nforct
	!if((adjustl(bilanh)=="ceres").OR.(adjustl(bilanh)=="pfactor")) THEN
	!else
	!write(*,*) "Nom de Bilan Hydrique inconnu"; stop
	!end if
	integer typeparam
	if ((irotation >1).OR.(irotation <0)) then
	write (*,*) "erreur: irotation different de 0 ou 1"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((iobs >1).OR.(iobs <0)) then
	write (*,*) "erreur: iobs different de 0 ou 1"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((iwatbal >1).OR.(iwatbal <0)) then
	write (*,*) "erreur: iwatbal different de 0 ou 1"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((inbal >1).OR.(inbal <0)) then
	write (*,*) "erreur: inbal different de 0 ou 1"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((typeparam>2).AND.(typeparam <0)) then
	write (*,*) "erreur: typeparam different de 0,1 ou 2"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((nforc1 <0)) then
	write (*,*) "erreur: nforc1 inferieur = 0"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if ((nforct <0)) then
	write (*,*) "erreur: nforct inferieur à 0"
	write(*,*)"Press any key to continue";read(*,*);stop 
	end if

	END SUBROUTINE

	SUBROUTINE verifparam1(typeparam,nbpar,varparam)
	!Verifie si lors de la lecture de sim.txt si la variable de calage est adequate
	character varparam*15
	integer nbpar, typeparam
	if ((nbpar>12).OR.(nbpar <0)) then
	    write(*,*) "erreur: nombre de paramètres"
		write(*,*)"Press any key to continue";read(*,*);stop 
	end if
	if (typeparam==2) THEN
		select case (adjustl(varparam))		
				Case ('stagnb','stahtvd',"blatnb","blagnb",'blagdm')
				Case ('blastvd','glai','ei','stamsu','stamsuc')
				Case ("stamfm",'stamdm','agrdm','agrwc','stock')							
				case default
					 write(*,*) "erreur: variable de calage"
					 write(*,*)"Press any key to continue";read(*,*)
					 stop				
		end select
	END IF !if (typeparam==2) THEN
	END SUBROUTINE
	
	SUBROUTINE verifforc(cc,s3)
	! vérifie, lors de la lecture de sim.txt si la variable de forcage existe
		character cc,s3*20
	    !write(*,*) s3
		If (cc=='1') Then
			select case (adjustl(s3))
					Case ('stagnb','stahtvd',"blatnb","blagnb")
					Case ('blastvd','glai','ei','stamsu')
					Case ("stamfm",'stamdm','agrdm','stock')					
				Case default
				  write(*,*) "erreur de variable de forcage date 1"
				  write(*,*)"Press any key to continue";read(*,*);stop
			end select
		end if
		If (cc=='t') Then
			select case (adjustl(s3))
					Case ('stagnb','stahtvd',"blatnb","blagnb")
					Case ('blastvd','glai','ei','stamsu')
					Case ("stamfm",'stamdm','agrdm','stock')					
				Case default
				  write(*,*) "erreur de variable de forcage dates"
				  write(*,*)"Press any key to continue";read(*,*);stop
			end select
		end if
	END SUBROUTINE





	SUBROUTINE valcol(k,chain,sep)
	!JFM020820  Verification: ok
	! Calcule la valeur de la colonne k de la chaine chain avec separateur sep
	! Attention la chaine chain est modifiée et retournée
	! Attention la valeur est fauuse si k > nbre total de colonnes
	!call calc_valcol(k,z,sep)
	CHARACTER(*):: chain
	CHARACTER sep	
	integer i,j,k
	do i=1,k
		chain=adjustl(chain)
		j=scan(chain,sep)
		if (i<k) chain=chain(j+1:)
	end do
	if (j>0) then
			chain=chain(:j-1)
	end if	
	END SUBROUTINE


	SUBROUTINE nbsep(nb,chain,sep)
	!Calcule le nombre de caractère sep dans une chaine chain
	!JFM020820  Verification: ok
	CHARACTER(*):: chain
	CHARACTER sep
	integer i,j,nb
	i=0
	do
		i=i+1
		chain=adjustl(chain)
		j=scan(chain,sep)
		if (j==0) exit
		chain=chain(j+1:)
	end do
	nb=i-1
	END SUBROUTINE

	SUBROUTINE integrcar(s1,s2)
	! vérifie que la chaine ne contient pas d'espace , / \ et ,
	CHARACTER(*):: s1
	CHARACTER(*):: s2
	integer i,j,k,l
	i=0; j=0; k=0; l=0
	i=scan(trim(s1)," ")
	j=scan(trim(s1),"/")
	k=scan(trim(s1),"\")
	l=scan(trim(s1),",")
	if ((i+j+k+l)>0) then 
	 write(*,*) trim(adjustl(s2))," incoherent: espaces ou , / ou \"; 
	 write(*,*) trim(adjustl(s1))
	 write(*,*) "Press any key to continue"; read(*,*); stop
	end if
	END SUBROUTINE


	

	subroutine joursuivant(idd,imm,iyr)
	! calcul le jour suivant (idd,imm,iyr) du jour (idd,imm,iyr)
	! attention les idd, imm et iyr changent
	! JFM verifie le 10/4/03
		integer idd, imm, iyr, id,im,iy,feb
		feb =28
	id=idd; im=imm; iy=iyr
		if (mod(iyr-1900,4)==0) feb=29
		Select case (imm)
			case (1,3,5,7,8,10,12)
				if (id<31) then
					id=id+1
					else
					id=1	
					if (im<12) then
						im=im+1
						else
						im=1; iy=iy+1
					end if
				end if
			case (4,6,9,11)
				if (id<30) then
					id=id+1
					else
					id=1; im=im+1
				end if
			case (2)
				if (id<feb) then
					id=id+1
					else
					id=1; im=3
				end if
		END SELECT
		idd=id 
		imm=im 
		iyr=iy
	end subroutine



