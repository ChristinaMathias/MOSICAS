c****************************************************
c	Subroutines	Definition
c****************************************************	
c	bladesnb0	Appearance, senescence and growth of green blades
c	convert1	Conversion into total dry mass
c	elong2		Stalks: growth of Tvd height
c	flor0		Flowering
c	humcan1		Above ground biomass water content
c	indic1		Indicators
c	indic2		Indicators
c	indic3		Indicators
c	intercep2	Interception efficiency
c	iwatnbal	water and nirogen stress index according to simulation options
c	laiglob0	Green leaf area index
c	partmst1	Partition of total dry mass into above ground dry mass
c	partsucre1	Partition of millable stalk dry mass into sucrose and structure
c	parttige1	Partition of above ground dry mass into millable stalk dry mass
c	rdtsucre2	Millable stalks: yield, sucrose and water contents
c	strhydrgen	water stress index
c	talcas1		Stalks: Appearance and senescence
c	thermo2		Degree days Calculation
c****************************************************


!********************* thermo2 **********************************/
		!thermo2: Degree days (xdj) Calculation
		!mode: calculation type of degree days xdj
	    !tempmn,tempmx: observed daily températures mini and maxi
		!tmo: calculated average température
		!ptmo: parameter to calculate average temperature (standard temperate = 0.5)
		!xdj: degree days of the day
	    !tb,topt1,topt2,tm: parameters: base, optimum (topt1=topt2) et maximum temperatures
		!p02,p03: additional parameters (not used)
		!tpent: coefficient calculatiing the evolution between tb and topt1. 1:linear
		!hourly types  x1: linéaire, x2:ceres x3: sigmoïde x4:Liu(non implémenté)
!**********************************************************************/
	SUBROUTINE thermo2 (mode,tempmn,tempmx,tmo,xdj,
	1                    tb,topt1,topt2,tm,tpent,ptmo,p02,p03)
		INTEGER mode, hour
		REAL tempmn,tempmx,tb,topt1,topt2,tm,tpent
		REAL v0,v1,v2,v3,ptmo,p02,p03
		INTEGER n
		REAL tmfac,tg,ta
		REAL tr
		!OUT
		REAL xdj
		xdj=0; hour=1
		tmo=ptmo*tempmx+(1-ptmo)*tempmn  !standard: ptmo=0.5
		select case (mode)
		 !********** ALGO using average temperature tmo
		 case (0)  !T moyenne algo simple
		 xdj = tmo - tb        
		 IF (xdj < 0) xdj=0
	     case (1) ! T moyenne Bornage de Tmin à tbase et Tmax à topt1 (Ecotrop)
		 xdj=(max(tempmn,tb)+min(tempmx,topt1))/2 - tb
		 IF (xdj < 0) xdj=0
		 case (2) !T moyenne   Ecotrop: module EvalDegresJourVitMoy 
		  !Hypothèse linéaire entre tn et tx, linéaire entre tb et top1,topt2 et tm, plant entre topt1 et topt2
		  v1=((max(tempmn,tb)+min(tempmx,topt1))/2-tb)/(topt1-tb)
		  v2=(tm-max(tempmx,topt2))/(tm-topt2)
		  v3=v1*(min(tempmx,topt1)-tempmn)
		  v3=v3+min(topt2,max(topt1,tempmx))-topt1
		  v3=v3+v2*(max(topt2,tempmx)-top2)
		  v3=v3/(tempmx-tempmn)
		  xdj=v3*(topt1-tb)
		  IF (xdj < 0) xdj=0
		 case (5) !T moyenne Tmoday-BetaFunction Beta fonction sur température moyenne Réponse non linéaire 
			xdj=0
	        if ((tmo>tb).and.(tmo<tm)) then
			v1=((tm-tmo)/(tm-topt1)); v2=((tmo-tb)/(topt1-tb))
			v3=((topt1-tb)/(tm-topt1))
			xdj=v1*v2**v3
			end if
			!write(*,*) tmo,xdj
		 !********** ALGO calculating and using hourly temperatures tr
		 case (11,12,13,21,22,23,41,42,43,51)    !T horaires
		  hour=24  ! nombre de scrutations en 24 heures
		  do n=1,hour
		  select case (mode)   !tr calculation
			case (41,11,21,51) !Thor linear  
	        tr=tempmn+(tempmx-tempmn)/hour*(n-0.5)
			case (42,12,22) !Thor Ceres  (hour=8, every 3 hours)
			tmfac=0.931+0.114*n-0.0703*(n**2)+0.0053*(n**3)	
			tr= tempmn + (tempmx-tempmn)*tmfac
			case (43,13,23)  !Thor  Sigmoïd from colimaçons (hour=8, every 3 hours
			tr= tempmn + (tempmx-tempmn)*tmfac
		  end select
		  select case (mode) ! xdj calculation from tr
			case (41,42,43)	! Thor Linéaire
			  if ((tr<=tb).or.(tr>=tm)) v0=0
			  if ((tr>tb).and.(tr<topt1)) then 
				v0 = (((tr-tb)/(topt1-tb))**tpent)*(topt1-tb)
			  end if
			  if ((tr>=topt1).and.(tr<=topt2)) v0=topt1-tb
			  if ((tr>topt2).and.(tr<tm)) then
				v0=(((tm-tr)/(tm-topt2))**1)*(topt1-tb)
			  end if
			xdj=xdj+v0/hour
			case (11,12,13)	!Thor-tblin
			if (tr<=tb) xdj=xdj+0
			if (tr>tb) xdj=xdj+(tr-tb)/hour
			!write(*,*) tempmn,tempmx,tr,xdj
			case (21,22,23)	!Thor-tbolin
			if (tr<=tb) xdj=xdj+0
			if ((tr>tb).and.(tr<=topt1)) xdj=xdj+(tr-tb)/hour
			if (tr>topt1) xdj=xdj+(topt1-tb)/hour
			case (51)	!JFM101110 ThorLin-BetaFunction
			v0=0
	        if ((tr>tb).and.(tr<tm)) then
			v1=((tm-tr)/(tm-topt1)); v2=((tr-tb)/(topt1-tb))
			v3=((topt1-tb)/(tm-topt1));
			v0=v1*v2**v3
			end if
			xdj=xdj+v0/hour
		  end select
		 end do
		end select
	END SUBROUTINE !thermo2




!********************* blades1 **********************************/
	!Partition of above ground dry mass into green blades dry mass
	!Sla calculation
	!gdmbaer: accroissement journalier de ms aérienne hors-sol (g/m2/d)
	!dmbaer: matière sèche aérienne hors-sol (g/m2)
	!gdmbl: accroissement journalier de ms de limbes titaux formés (g/m2/d)
	!dmtbl: matière sèche de limbes totaux formés (g/m2)
	!slam: surface massique des limbes totaux formés (m2/g)
	!pblfin: parameter, taux de partition final aux limbes
	!pbldif: parameter, différence de partition entre le début et la fin
	!pblkdec: parameter, affectant le taux de partition aux limbes
	!slafin: parameter, surface massique finale
	!sladif: parameter, différence surface massique entre le début et la fin
	!slakdec: parameter,  affectant l'evolution de la surface massique des limbes
!**********************************************************************/
	SUBROUTINE blades1(pblfin,pbldif,pblkdec,slafin,sladif,slakdec,
	1					dmbaer,gdmbaer,dmtbl,gdmbl,slam)
		  real pblfin,pbldif,pblkdec,slafin,sladif,slakdec
		  real dmbaer,gdmbaer
		  real dmtbl,gdmbl,slam
			if (gdmbl>0) then
			  gdmbl=(pblfin+pbldif*exp(-pblkdec/1000*dmbaer))*gdmbaer !30/05/02 0.17  0.394  0.00085
			  dmtbl=dmtbl+gdmbl !30/05/02
			  slam=(slafin+sladif*exp(-slakdec/1000*dmtbl))/10000  !92.3  47.2  0.035
			end if
	END SUBROUTINE blades1


!********************* flor0 *************************************/
	!Flowering calculation
	!nage: age (days)
	!pflor,iflor: % of flowering
	!stflor: flowering stage (0-1)
	!stmo: cumul of temperatures (Tx+Tn)/2 since p,ating/harvest
	!astdurjour,asdurjourp: day duration of current and previous days
!*****************************************************************/
	subroutine flor0 (nage,pflor,stmo,astdurjour,astdurjourp,stflor,
	1		iflor,niflor)
	real pflor
	real stmo, astdurjour,asdurjourp
	real iflor
	integer stflor, nage, i, jcal, niflor, niflorprec
	niflorprec=niflor

	if (stmo>2500) then
	if ((astdurjour<=12.7).and.(astdurjour>=12.3)) then
	if (astdurjour<astdurjourp) then
		stflor=1; niflor = niflor + 1.
	end if
	end if
	else
	niflor=0; stflor=0
	end if
	if ((niflor>1).and. (niflor>niflorprec))  then 
	iflor=pflor*stflor 
	end if
	end subroutine flor0
!*****************************************************************/


!********************* laiglob0 **********************************/
	!Green leaf area index calculation
	! Calculate glai without using mass and population  JFM051128
	!mode: type of calculation of degree days xdj (not used here)
	!iflor: % flowering 
	!tempmn,tempmx,laitb,laicroi,laiwksen,iflor,swdf2,sdjlai,lai,nbtigv
!**********************************************************************/
	subroutine laiglob0(tempmn,tempmx,laitb,laicroi,laiwksen,
	1		   iflor,swdf2,sdjlai,lai,nbtigv)
	real laitb,laicroi,laiwksen,laimax !Input parameters
	real tempmn,tempmx,swdf2,iflor,nbtigv  !Input variables
	real glai, slai, slainh2o,tmo,xdj,xdjlai,pari  !temporay variables
	real sdjlai,lai    !Input/Output state variables
		glai = 0.; slai = 0.;slainh2o=0. ;xdjlai=0.; laimax=7.     
		if (nbtigv>0.) then
	  	 if (lai==0.) then
		 lai=0.05 !0.1
		 end if
		 xdjlai = (tempmx+tempmn)/2 - laitb; xdjlai=max(0.,xdjlai)
		 glai=(laicroi)*lai*(1.-lai/laimax)*xdjlai*swdf2 !previous version (gompertz derivative)
		 glai=glai*(100-iflor)/100  !flowering effect
		 slainh2o=laiwksen*(1-swdf2)*lai 
	  	 slai=max(slai,slainh2o)
      	 lai=lai+glai-slai
		 if (lai<0.) lai=0.		
		end if
	end subroutine laiglob0
!**********************************************************************/


!********************* intercep2 **********************************/
	!calculation of interception efficiency (ei) 
	!calculation of intercepted photosynthetical active radiation (pari)
	!ke: parameter, extenction coefficient
	!solrad, par: incident global radiation and intercepted photosynthetical radiation
	!lai: green leaf area index
!**********************************************************************/
	SUBROUTINE intercep2(isobei,obei,ke,lai,solrad,ei,par,pari)  !jfm121219
	!in
	REAL ke,lai,ec,solrad
	!temp
	REAL cske
	!out
	REAL ei,par,pari
		cske=1
		ei = 1. - Exp(-(cske * ke * lai))
		if (isobei==1) then
			ei=obei; lai=-log(1-ei)/ke
		end if		
		par=0.5*solrad; pari=0.5*solrad*ei

	END SUBROUTINE
!**********************************************************************/

!********************* convert1 **********************************/
! Conversion of daily intercepted photosynthetic radiation (pari) into 
!    daily total dry mass accumulation (gmst)
!astrgx : Global extraterrestrial  radiation (MJ/m2/d)
!pari: daily intercepted photosynthetic radiation (MJ/m2/d) 
!tempmn,tempmx : daily min and max temperatures	
!ruemax: conversion coefficient of intercepted photosynthetic radiation into total dry mass (gr/MJ)
!rueopt: Optimum temperature for conversion
!ruetk: Effect of temperature on conversion
!solrad: Global Incident radiation (MJ/m2/d)
!swdf1: water stress index for dry mass accumulation (0-1)
!gmst: daily total dry mass accumulation (g/m2/d)
!mst: total dry mass (g/m2)
!p01: coefficient for maintenance effect on conversion
!p04: coefficient for diffused radiation effect on conversion
!**********************************************************************/
	SUBROUTINE convert1(ruemax,ruetk,ruetopt,tempmn,
	1	tempmx,pari,solrad,astrgx,swdf1,gmst,mst,p01,p04)
	INTEGER mainttype !IN  mode de calcul de l'effet age
	REAL tempmn,tempmx,pari,swdf1,p01,p04  !IN
	REAL tdiu, ktemp,rm, ruex,rap   !TEMP
	REAL gmst,mst,kage  !OUT
		gmst=0; kage=1; ktemp=1
		! Temperature effect ktemp (Kiniry)
		 tdiu=0.25*tempmn+0.75*tempmx ! tdiu: température moyenne diurne	
		 ktemp = 1.-ruetk*(abs(tdiu-ruetopt))**0.5 !ktconv,tconvopt   0.0015, 30 jfm160209 corrigé (abs)
		 ktemp=max(0.,min(1.,ktemp))		
		! Diffuse radiation effect
		 rap=min(1.0, solrad/astrgx)
		 ruex=ruemax
		 pari=pari*(1-p04*(rap-0.75)) !jfm160301
		! Maintenance effect (kage)				
		 kage=1-(mst*p01/10000)**3;
		 kage=max(0.,min(1.,kage))
		! Conversion
		   gmst=ruex*kage*pari*swdf1*ktemp
		mst=mst+gmst
	END SUBROUTINE
!**********************************************************************/



!********************* partmst1 **********************************/
	!Partition of total dry mass into above ground biomass
	!ratoon: crop cycle (0: plant crop, 1: first ratoon, ..)
	!gmst: daily total dry mass accumulation (g/m2/d)
	!mst: total dry mass (g/m2)
	!sdj: Thermal time since planting or previous harvesting (base temperature=0)
	!gdmrac: daily root dry mass accumulation (g/m2/d) (g/m2/d)
	!dmrac: root dry mass (g/m2)
	!gdmbaer: daily above ground dry mass accumulation (g/m2/d) (g/m2/d)
	!dmbaer: above ground dry mass (g/m2)
!**********************************************************************/
	SUBROUTINE partmst1(ratoon,gmst,sdj, mst,
	1			gdmrac,dmrac,gdmbaer,dmbaer)
	INTEGER ratoon
	REAL pracdeb,pracfin,pracdec,gmst,sdj, mst
	REAL krac
	REAL gdmrac,dmrac,gdmbaer,dmbaer
	gdmrac=0; gdmbaer=0
	pracdeb=0.7;  ! initial allocation to roots
	pracfin=0.1;  ! final allocation to roots
	pracdec=0.001 ! attenuation coefficient for allocation to roots
		if (ratoon>=1) then
			krac=pracfin
		  else
			krac=(pracdeb-pracdec*sdj)
			if (krac<pracfin) krac=pracfin
		end if
		gdmrac=krac*gmst
		dmrac=dmrac+gdmrac			
		IF ((mst-dmrac)>500) Then ! allocation to above ground mass
			gdmbaer=(gmst-gdmrac) 
		ELSE 
			gdmbaer=(0.7+(mst-dmrac)/500*0.3)*(gmst-gdmrac)
		END IF
		dmbaer=dmbaer+gdmbaer
	END SUBROUTINE
!**********************************************************************/


!********************* parttige1 **********************************/
!Partition of above ground dry mass into millable stalk dry mass
!gdmbaer: daily above ground dry mass accumulation (g/m2/d) (g/m2/d)
!dmbaer: above ground dry mass (g/m2)
!jat: age when millable stalk appearance
!ptigdeb:Beginning of millable stalk dry mass appearance (Plant parameter (g/m2)
!ptigdec: Extinction coefficient of  daily fraction of aboveground dry mass 
!            allocated to millable stalk  (Plant parameter)
!ptigfin	Final  daily fraction of aboveground dry mass allocated to 
!             millable stalk (Plant parameter)
!gdmstm: daily increase in millable stalk dry mass (g/m2/d)
!dmstm: millable stalk dry mass (g/m2)
!**********************************************************************/
	SUBROUTINE	parttige1 (ptigdeb,ptigfin,ptigdec,dmbaer,gdmbaer,
	1						jat,gdmstm,dmstm)
	REAL ptigdeb,ptigfin,ptigdec,dmbaer,gdmbaer !IN
	REAL gdmstm,dmstm !OUT
	INTEGER jat
	  if (dmbaer>ptigdeb) then !debmstu=500 le 28/04/02
	    jat=jat+1
		gdmstm=ptigfin*(1-exp(-ptigdec/1000*(dmbaer-ptigdeb)))*gdmbaer
		else
		jat=0;gdmstm=0
	  end if  
		dmstm=dmstm+gdmstm
	END SUBROUTINE
!**********************************************************************/


!********************* partsucre1 **********************************/
!Partition of millable stalk dry mass into sucrose and structure
!gdmstm: daily increase in millable stalk dry mass (g/m2/d)
!dmstm: millable stalk dry mass (g/m2)
!tempmn,tempmx : daily min and max temperatures	
!pstrudec	Extinction coefficient of  daily fraction of millable stalk dry mass 
!            allocated to structures   Plant parameter
!pstrufin	Final daily fraction of millable stalk dry mass allocated to structures 
!            Plant parameter  g/g	
!pstrutb	Temperature treshold from which fraction of millable stalks dry mass 
!            allocated to structures is decreasing.  Plant parameter
!pstrutcroi	Temperature effect on daily fraction of millable stalks dry mass 
!            allocated to structures  /°C.	  Plant parameter
!gdmsug: daily increase in succrose mass (g/m2/d)
!dmsug: succrose mass (g/m2)
!swdf1: water stress index for dry mass accumulation (0-1)
!swdf2: water stress index for growth (0-1)
!**********************************************************************/
	SUBROUTINE	partsucre1(tempmn, tempmx,gdmstm,dmstm,swdf1,swdf2,
	1                pstrutb,pstrutcroi,pstrufin,pstrudec,gdmsug,dmsug)
	REAL tempmn, tempmx,gdmstm, tmo, dmstm, swdf1, swdf2 !IN
	REAL pstrutb,pstrutcroi,pstrufin,pstrudec
	REAL gdmstst, swdfpart1 !TEMP
	REAL gdmsug,dmsug !OUT
		gdmstst=gdmstm	   
		swdfpart1=1
		tmo=(tempmn+tempmx)/2
		if (tmo>pstrutb) swdfpart1=1+pstrutcroi*(tmo-pstrutb)  !0.03
	    gdmstst=(1-pstrufin*(1-exp(-pstrudec/1000*dmstm)))*gdmstm !0.59  6
		gdmstst=gdmstst*swdfpart1*(1-(swdf1-swdf2))
		if (gdmstst>gdmstm) gdmstst=gdmstm
		gdmsug=gdmstm-gdmstst
		dmsug=dmsug+gdmsug
	END SUBROUTINE
!**********************************************************************/


!********************* rdtsucre2 **************************************/
!Calculation of Millable stalks:fresh mass (yield), sucrose and water contents
!dmstm: millable stalk dry mass (g/m2)
!hum_stm millable stalk water content (%)
!dmsug: succrose mass (g/m2)
!yldcan, Yield or millable stalk fresh mass (T/Ha)
!sug_stm, millable stalk fresh mass sucrose content (%)
!jat: age when millable stalk appearance
!**********************************************************************/
	SUBROUTINE rdtsucre2(dmstm,hum_stm,dmsug,yldcan,sug_stm,jat)
		REAL dmstm,dmsug,watstm !IN
		integer jat !IN
	    REAL hum_stm
		REAL yldcan, sug_stm !OUT
	if (dmstm>0) then
		hum_stm = 88 - 0.070374 * jat  !rmse=8.14 Hum%tige=85.77-0.603*agetu (R570, Rep, CalReu) 121107o  (jfm121109)
		hum_stm=min(hum_stm,90.); hum_stm=max(hum_stm,68.)
		yldcan=dmstm/(100.-hum_stm)
		if (yldcan>0.1) sug_stm=dmsug/yldcan
	end if
	END SUBROUTINE
!**********************************************************************/


!********************* humcan1 **********************************/
!Above ground biomass water content	
!humaerdeb: initial ground biomass water content % Plant parameter
!humaerfin: begining in decrease of ground biomass water content  Plant parameter
!humaertb: Base temperature for ground biomass water content  Plant parameter
!humaerdec: Rate of decrease in ground biomass water content  Plant parameter
!hum_aer: above ground biomass water content
!tempmn,tempmx : daily min and max temperatures	
!humaersdd: Thermal time for ground biomass water content
!dmbaer: above ground dry mass (g/m2)
!**********************************************************************/
	SUBROUTINE humcan1(humaerdeb, humaerfin, humaertb, humaerdec,
	1			hum_aer,tempmx,tempmn,humaersdd,dmbaer) 
		REAL humaerdeb, humaerfin, humaertb, humaerdec !IN parametres
		REAL tempmx,tempmn,dmbaer    !IN
		REAL  dd !TEMPO
		REAL hum_aer                 !OUT
		REAL humaersdd               !IN-OUT
		if (dmbaer<=0.) then
			humaersdd=0.; hum_aer=humaerdeb
		  else
			dd=max(0.,(tempmx+tempmn)/2-humaertb)
			humaersdd=humaersdd + dd
			if (humaersdd<humaerfin) then
				hum_aer=humaerdeb 
			  else
				hum_aer=hum_aer - dd * humaerdec				
			end if
			hum_aer=max(50.,hum_aer)
		end if		
	END SUBROUTINE
!**********************************************************************/

!********************* iwatnbal **************************************/
!Updates of water and nirogen stress index according to simulation options
!iwatbal Simulation option for water stress (0 water stress is not taken into account
!     or 1 water stress is taken into account)
!inbal   Simulation option for nitrogen stress (0 nitrogen stress is not taken into account
!     or 1 nitrogen stress is taken into account)
!swdf1: water stress index for dry mass accumulation (0-1)
!swdf2: water stress index for growth (0-1)
!swdef  water stress index of soil (0-1)
!ndef1  nitrogen stress index for dry mass accumulation (0-1)
!ndef2  nitrogen stress index for growth (0-1)
!**********************************************************************/
	SUBROUTINE iwatnbal(iwatbal,inbal,swdf1,swdf2,swdef,ndef1,ndef2)

		INTEGER iwatbal, inbal !IN
		REAL swdf1,swdf2,swdef,ndef1,ndef2
		IF (iwatbal==0) THEN
				swdf1=1; swdf2=1; swdef=1
		END IF
		IF (inbal==0) THEN
				ndef1=1; ndef2=1
		END IF
	END SUBROUTINE
!**********************************************************************/




!********************* indicl **************************************/
	SUBROUTINE indic1 (swdef,sswdef,swdf1,sswdf1,swdf2,sswdf2,
	1	pari,spari,astrgo,astrgx,srgo,srgx)		
		REAL swdef,swdf1,swdf2,pari,astrgo,astrgx !IN		
		REAL sswdef, sswdf1, sswdf2,spari,srgo,srgx !OUT
		sswdef=sswdef+swdef
		sswdf1=sswdf1+swdf1
		sswdf2=sswdf2+swdf2
		spari=spari+pari
	    srgo=srgo+astrgo; srgx=srgx+astrgx
	END SUBROUTINE
!**********************************************************************/


!********************* indic2 **************************************/
	SUBROUTINE indic2 (jat,bldj,sdjblt,htvddj,sdjhtvdt,swdef,sswdeft,
	1	swdf1,sswdf1t,swdf2,sswdf2t,tmo,stmot,par,spart,pari,sparit,
     2    etm,setmt,etr,setrt)		
		integer jat
		REAL bldj,htvddj,swdef,swdf1,swdf2,tmo,par,pari !IN		
		REAL sdjblt,sdjhtvdt,sswdeft,sswdf1t,sswdf2t,stmot,spart,sparit !OUT
		if (jat>=1) then
			sdjblt=sdjblt+bldj
			sdjhtvdt=sdjhtvdt+htvddj
			sswdeft=sswdeft+swdef
			sswdf1t=sswdf1t+swdf1
			sswdf2t=sswdf2t+swdf2
			stmot=stmot+tmo
			spart=spart+par
			sparit=sparit+pari
			setmt=setmt+etm 
			setrt=setrt+etr
		end if !if (jat>=1) then
	END SUBROUTINE
!**********************************************************************/



!********************* indic3 **************************************/
	SUBROUTINE indic3 (jat,etm,setm,etr,setr,tmp,stmp,trp,strp,runoff,
	1 srunoff,es,sevap,drain,sdrain,sdraint,srunofft,
     2 eeq,setpt,stmpt,strpt)
	real etm,setm,etr,setr,tmp,stmp,trp,strp,eeq,setpt,stmpt,strpt
	real runoff,srunoff,es,sevap,drain,sdrain
	setm=setm+etm 
	setr=setr+etr
	stmp=stmp+tmp
	strp=strp+trp 
	srunoff=srunoff+runoff 
	sevap=sevap+ES 
	sdrain=sdrain+drain
	if (jat>0) then
	sdraint=sdraint+drain
	srunofft=srunofft+runoff
	setpt=setpt+eeq
	stmpt=stmpt+tmp
	strpt=strpt+trp
	else
	sdraint=0.; srunofft=0.; setpt=0.; stmpt=0.; strpt=0.
	end if
	END SUBROUTINE
!**********************************************************************/



!********************* BLADESnb0 **********************************/
!Blades appearance, senescence and tvd surface on primary stalks
!tempmn, tempmx: daily mini and maxi temperatures
!bltb: base temperatures for blades appearance (plant parameter)
!blk   rate of blade appearance (plant parameter)
!blkdec: attenuation coefficient for blade appearance (plant parameter)
!blvmax: maximum green blades number (plant parameter)
!blsurfk area growth rate of tvd blades according to blades number (plant parameter)
!blsurfmax: maximun area of tvd blades (plant parameter)
!xdj,blsdj: degree days and thermal time for blades appearance
!blnbtot  number of total blades on primary stalks
!blnbv:   number of green blades on primary stalks
!bltvdsurf: area of tvd blade on primary stalks
!**********************************************************************/
	SUBROUTINE	bladesnb0(tempmn,tempmx,bltb,blk,blkdec,
	1			blvmax,blsurfk,blsurfmax,nbtigv,blnbtot,blnbv,
     2            blsdj,bltvdsurf)  !JFM170418)
	implicit none
	real tempmn,tempmx,nbtigv      !input variables	
	real bltb,blk,blkdec           !input param apparition
	real blvmax,blsurfk,blsurfmax  !input param tvdsurf
	real blnbtot,blnbv,blsdj,bltvdsurf  !input/output variables
	real xdj, gblnbt  !tempo
	gblnbt=0;
	if (nbtigv>0.) then  !jfm170408
	  !Thermal times
	    xdj=max((tempmn+tempmx)/2-bltb,0.); blsdj=blsdj+xdj;
	  !blades appearance on primary stalk (param: blk, blkdec)  
	     blnbtot= blk*(blsdj**blkdec)  !power law function: y=0.1*(x**0.7) Bonnet
	  !senescence on primary stalk (param: blvmax)
	     blnbv=min(blnbtot,blvmax)
	  !surface of tvd blade (param= blsurfk, blsurfmax)
	     bltvdsurf=min(blsurfk*blnbtot,blsurfmax)
	else
	blnbtot=0; blnbv=0; blsdj=0
	end if
	END SUBROUTINE bladesnb0
!**********************************************************************/




!********************* TALCAS1 **********************************/
!talcas1 alive stalks number calculation 
!Tempmn et tempmx: températures minimum et maximum
!taltb: base temperature for stalk appearance  (plant parameter)
!taldebtt: TT (termal time) to reach appearance of first stalks (plant parameter)
!talpeaktt: TT to reach tillering peak (plant parameter)
!talpeakval: Value of tillering peak (plant parameter)
!talfinval: Final value of alive stalks (/m2) (plant parameter)
!talsdj: TT of stalk appearance
!xdj: degree days of the day
!nbtigv: Alive stalks per m2
!**********************************************************************/
	SUBROUTINE talcas1(tempmn,tempmx,talsdj,nbtigv,ratoon, taltb, 
	1           taldebtt,talpeaktt,talpeakval,talfinval) 
	real tempmn,tempmx  !input variables
	real taltb,taldebtt,talpeaktt,talpeakval,talfinval  !input parameters
	real nbtigv,talsdj,xdj,tmo  !input/output variables
	xdj=max(0.,(tempmn+tempmx)/2-taltb)
	talsdj=talsdj+xdj	
	if (talsdj>=taldebtt) then  !après levee
		if (talsdj<talpeaktt) then
			nbtigv=nbtigv+talpeakval/(talpeaktt-taldebtt)*xdj
		endif
		if (talsdj>=talpeaktt) then
		if(talsdj<=talpeaktt*1.5) then
			nbtigv=nbtigv-(talpeakval-talfinval)/(0.5*talpeaktt)*xdj
		end if
		end if
	 else
	  nbtigv=0.
	end if
	END SUBROUTINE talcas1
!**********************************************************************/


!********************* ELONG2 **********************************/
!/elong2 Calculation of stalks TVD height  (primary stalks)
!tempmn, tempmx: daily mini and maxi temperatures
!dhtvd daily growth rate of tvd height (cm/d)
!htvd; tvd height of primary stalks (cm)
!htvdtb,htvdto,htvdtm: base, optimum and maximum temperatures (plant parameters)
!htvdcan: TVD height stage after which elongation rate is greater (plant parameter)
!htvdgro1: first rate of elongation (cm/°J)  (plant parameter)
!htvdgro2: second rate of elongation after htvdcan stage (cm/°J)  (plant parameter)
!nbtigv: alive stalks /m2 
!xdj: degree days
!dhtvd,htvd: daily elongation and tvd height  (cm)
!swdf2: water stress index for growth (0-1)
!modexdj: type of calculation of degreedays (xdj)
!**********************************************************************/
	SUBROUTINE elong2(tempmn,tempmx,htvdtb,htvdto,htvdtm,htvdcan,
	1		htvdgro1,htvdgro2, swdf2, nbtigv,htvd,dhtvd,modexdj)
		Integer modexdj  !input
		REAL tempmn,tempmx,htvdtb,htvdto,htvdtm !IN
		REAL htvdcan,htvdgro1,htvdgro2, swdf2 !IN
		REAL xdj,tmo,htvd,dhtvd, nbtigv !TEMP
		call thermo2 (modexdj,tempmn,tempmx,tmo,xdj,
	1                   htvdtb,htvdto,htvdto,htvdtm,1.,0.5)
	 dhtvd=0.
	 if (nbtigv>0.) then  !après levee	   		
	  if (htvd<htvdcan) then
		dhtvd=xdj*htvdgro1*swdf2
	   else
		dhtvd=xdj*htvdgro2*swdf2
	  end if
	  htvd=htvd+dhtvd
	 end if
	END SUBROUTINE elong2
!**********************************************************************/




!********************* strhydrgen **********************************/
!Calculation of water stress index swdf1 and swdf2 according to several types (mode)
!mode: type of calculation of swdf1 and swdf2
!swdef  water stress index of soil (0-1)
!swdf1: water stress index for dry mass accumulation (0-1)
!swdf2: water stress index for growth (0-1)
!tmp; Maximum transpiration of the day according to actual canopy (mm)
!trp; actual transpiration (mm)
!**********************************************************************/
	SUBROUTINE strhydrgen(mode,tmp,trp,swdef,swdf1,swdf2,psilb)
		INTEGER mode
		REAL tmp,trp, swdef,sthydcroi,sthydbio
		REAL swdf1, swdf2, psilb
		if (swdef==0) swdef=0.001
		!swdf1=1.; swdf2=1.   !JFM160921  mis en commentaire 
		select case (mode)
		  case (1) ! Brisson 1992 and  Slabber 1982  8.5  (Ozier 12)
			if (tmp==0) then
			srue=0.  ! JFM050211 : ajout si tmp=0
			else
			srue=(0.94-0.26*psilb/tmp)  !psilb=8 for R570 0.94- car normalement psilb<0
			end if
			if (swdef>(1.5*srue)) then   !JFM160921  correction <= (1.5*srue)
				swdf2=1
	            else
				swdf2=swdef/(1.5*srue)
			end if
			if (swdef>srue) then 
				swdf1=1
	            else
				swdf1=swdef/srue
			end if
		  case (2) ! swdf1 et swdf2 donnés par Ceres
			swdf1=swdf1
			swdf2=swdf2
		  case (3) !utilisation de p0
		  case (4) !utilisation de trp/tmp ou cstr
			if (tmp>0) swdf1=trp/tmp
			swdf2=swdf1 
	      case (5) !Eigelmann
		end select
		if (swdf1<0) swdf1=0
		if (swdf1>1) swdf1=1		
		if (swdf2<0) swdf2=0
		if (swdf2>1) swdf2=1
	END SUBROUTINE strhydrgen
c********************************************************************************




	 


