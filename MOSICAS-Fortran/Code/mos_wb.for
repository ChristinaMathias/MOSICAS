c****************************************************
c     WATER BALANCE MODULES
c****************************************************
c astrono1 : astronomic calculations
c WATBAL           : CERES water balance
c wbinit1          : FAO56 WB - variables initialisation
c wbkce            : FAO56 WB - calculation of evaporation index
c wbkcp            : FAO56 WB - calculation of transpiration index
c wbevapotfao      : FAO56 WB - calculation of maximum avapotranspiration
c wbevapfao        : FAO56 WB - calculation of actual avapotranspiration
c wbcstranspfao    : FAO56 WB - calculation of maximum transpiration
c wbtranspfao      : FAO56 WB - calculation of actual transpiration
c wbtransfert1fao  : FAO56 WB - Calculation of water flows between layers
c wbswdef          : FAO56 WB - calculation of water stress index
c****************************************************


C****   Astronomic calculations JFM130513 *****
	SUBROUTINE astrono1(nage,jcal,lati,longi,alti,astlatrad,astdr,
	1	astdecl,astomeg,astrgo,astrgx,astdurjour,astdurjourp)
		integer nage,jcal
		real lati,longi,alti,astlatrad,astdr,astdecl
		astlatrad= lati * 3.14159265358979 / 180				
	    astdr=1+0.033*COS(0.0172*jcal)							
		astdecl=0.409*SIN(0.0172*jcal-1.39)						
	    astomeg=ACOS(-TAN(astlatrad)*TAN(astdecl))				
		astrgo=cos(astdecl)*cos(astlatrad)*sin(astomeg)			
		astrgo= 37.6*astdr*(astomeg*sin(astdecl)*sin(astlatrad)+astrgo)
		astdurjourp=astdurjour
		astdurjour=7.64*astomeg									
	    astrgx=(0.75+0.00002*alti)*astrgo						
	END SUBROUTINE !&strono1



      SUBROUTINE WATBAL
        INCLUDE 'cercane1.inc'  !JFM000301
		if (adjustl(typsol)=="latosol") then
		do l=1,nlayr
			  if ((sw(l)-dul(l))>0.2*(sat(l)-dul(l))) then
			  sw(l)=dul(l)+0.2*(sat(l)-dul(l))
			  end if
		end do
		end if
	icsdur=icsdur+1 
      DRAIN=0.
      PINF=0.
      RUNOFF=0.
      stocki=0.
      do 109 l=1,nlayr
         stocki=stocki+dlayr(l)*sw(l)*10
         flux(l)=0.
	!write(*,*) l,ll(l),dul(l),sat(l),sw(l)
 109      continue
c assuming zero runoff
      WINF=0.
c precip = daily rainfall (mm)
      IF (PRECIP.EQ.0.) GO TO 110
  110 PINF=PRECIP+irdose-RUNOFF
	
C********** CALCULATE DRAINAGE AND SOIL WATER REDISTRIBUTION *******
      WINF=PINF
      FLUX(1)=PINF*0.1
	lai=3
	write(*,*) swdf1,swdf2
c flux(1) is water flux infiltrating in the top layer (cm/day)
	drl1=0.
c drl1 is drainage flux from first layer to 2nd layer (cm/d)
       sw(1)=sw(1)+flux(1)/dlayr(1)
       flux(1)=0.
	if (sw(1).gt.sat(1)) then
		drl1=(sw(1)-sat(1))*dlayr(1)
		sw(1)=sat(1)
       	l=1 
	j=1
        do while ((drl1.gt.0.) .AND. ((l+j).lt.nlayr))
                sw(l+j)=sw(l+j)+drl1/dlayr(l+j)
                if ((sw(l+j)).gt.(sat(l+j))) then
                        drl1=(sw(l+j)-sat(l+j))*dlayr(l+j)
                        sw(l+j)=sat(l+j)
               	else 
                   drl1=0. 
		endif
c drl1 is drainage flux from layer j to below layer j+1 (cm/d)
                j=j+1
        end do
	flux(nlayr)=flux(nlayr)+drl1
c flux(nlayr) =  drainage from the bottom of the profile
	drl1=0.
	endif
c       assuming all incoming water can infiltrate
c	excess of SAT is by-pass flow
c 
        DO 240 L=1,NLAYR
c 
        HOLD=(SAT(L)-SW(L))*DLAYR(L)
c 
c    computing drain using Gardner's approximation 
c     here drain occurs under gravity forces only
c       and each layer loses water above DUL
c         independently (within the one day time step
c           from what may drip from the upper layer
c
        if (sw(l).lt.dul(l)) then
                flux(l)=0.
                go to 240
                endif 
        dbar=SWCON*exp(aaa*(SW(L)-SAT(L)))
 	if (dbar.gt.SWCON) dbar=SWCON   
        DRAIN=dbar
        if (drain.gt.((sw(l)-dul(l))*dlayr(l))) 
     1 drain=((sw(l)-dul(l))*dlayr(l))
      	if (l.lt.nlayr) then 
	 if (drain.gt.((sat(l+1)-sw(l+1))*dlayr(l+1))) then 

        drain=(sat(l+1)-sw(l+1))*dlayr(l+1)
	endif
	endif

        sw(l)=sw(l)-drain/dlayr(l)      
        flux(l)=flux(l)+drain
        GO TO 240
  240 CONTINUE
      IF (L.GE.NLAYR) L=NLAYR
      DRAIN=FLUX(L)*10.0
      do 241 l=2,nlayr
        sw(l)=sw(l)+(flux(l-1)/dlayr(l))
  241 continue 
      stockf=0.
      do 249 l=1,nlayr
         stockf=stockf+dlayr(l)*sw(l)*10
 249      continue
  250 CONTINUE	


C********** CALCULATE POTENTIAL EVAPORATION ************************
        eo=eeq  
  302 EOS=EO*(1.-0.43*LAI)   !Ceres Original
      IF (LAI.GT.1.) EOS=EO/1.1*EXP(-0.4*LAI)   !Ceres original
C********** CALCULATE SOIL EVAPORATION *****************************
	
      IF (SUMES1.GE.U.AND.WINF.GE.SUMES2) GO TO 400
      IF (SUMES1.GE.U.AND.WINF.LT.SUMES2) GO TO 410
      IF (WINF.GE.SUMES1) GO TO 440
      SUMES1=SUMES1-WINF
      GO TO 450
  400 IF (WINF.LT.SUMES2) GO TO 410
      WINF=WINF-SUMES2
      SUMES1=U-WINF
      T=0.
      IF (WINF.GT.U) GO TO 440
      GO TO 450
  410 T=T+1.
      ES=alf*T**0.5-SUMES2
      IF (WINF.GT.0.) GO TO 420
      IF (ES.GT.EOS) ES=EOS
      GO TO 430
  420 ESX=0.8*WINF
      IF (ESX.LE.ES) ESX=ES+WINF
      IF (ESX.GT.EOS) ESX=EOS
      ES=ESX
  430 CONTINUE
      SUMES2=SUMES2+ES-WINF
      T=(SUMES2/alf)**2
      GO TO 470
  440 SUMES1=0.
  450 SUMES1=SUMES1+EOS
      IF (SUMES1.GT.U) GO TO 460
      ES=EOS
      GO TO 470
  460 ES=EOS-0.4*(SUMES1-U)
      SUMES2=0.6*(SUMES1-U)
      T=(SUMES2/alf)**2
  470 SW(1)=SW(1)-ES*.1/DLAYR(1)
      IF (SW(1).GE.ll(1)*SWEF) GO TO 480
      ES1=(ll(1)*SWEF-SW(1))*DLAYR(1)*10.
      SW(1)=ll(1)*SWEF
      ES=ES-ES1


C********** CALCULATE UNSATURATED FLOW BELOW DRAINED UPPER LIMIT ***
c
C********   k-theta relation
C********** and van Keulen and Wolf (eds, 1986), for the
C********** h-theta relation       (dbar in cm/d)
C**********       ....in order to apply a coherent Darcy-Richards'law...
C
  480 NIND=NLAYR-1
	!if (nage<50) write(*,*) nage, 10*(sw(1)-ll(1))*dlayr(1),eos,es
      DO 490 L=1,NLAYR
        FLOW(L)=0.0
        SWX(L)=SW(L)
  490 CONTINUE
      IST=1
        IF (DLAYR(1).EQ.5.0) IST=2
      DO 500 L=IST,NIND
        if (flux(l).gt.0) go to 500     
        MU=L+1
        if (sw(l).ge.sat(l)) sw(l)=sat(l)!-0.0001      !jfm170514  -0.0001
        if (sw(mu).ge.sat(mu)) sw(mu)=sat(mu)!-0.0001  !jfm170514  -0.0001   
        THET1=SW(L)-sat(l)
c        IF (THET1.LT.0.) THET1=0.
        THET2=SW(MU)-sat(mu)
	  
        DBAR=swcon*EXP(aaa*(THET1+THET2)*0.5)
C******* correlation coefficients for Yolo loam  (gam,aaa)       
        IF (DBAR.GT.swcon) DBAR=swcon
        avth=(DLAYR(L)+DLAYR(MU))*0.5
        psi(l)=exp(sqrt(log(sat(l)/abs(sw(l)))/gam))
	  psi(mu)=exp(sqrt(log(sat(mu)/abs(sw(mu)))/gam))  
        if ((psi(l)-psi(mu)-avth).gt.0.) then
        FLOW(L)=DBAR*(psi(l)-psi(mu)-avth)/avth
        endif
        WAT1=DUL(1)-SW(1)
        IF (FLOW(1).GT.WAT1) FLOW(1)=WAT1
        IF (WAT1.LT.0.0) FLOW(1)=0.0
c
c       computing capillary rise according to tables by Rijtema (1969)
c       (mettre tout a zero dans la boucle 509 si pas de remontees de nappe
c       
        SWX(L)=SWX(L)+FLOW(L)/DLAYR(L)
        SWX(MU)=SWX(MU)-FLOW(L)/DLAYR(MU)
c       
c       same iterative computation as in drain
c         the flows fill the upper layer (index L)
c           but empty the lower layer at the end of
c             procedure only
c  
  500 CONTINUE
        mu=nlayr        
       go to 509       
C by-passe les remontées capillaires
        if (flux(nlayr).gt.0) go to 509         
        if (psi(mu).lt.250) then 
                flow(mu)=0.0
                else if (psi(mu).lt.500) then
                        flow(mu)=0.03 
                        else if (psi(mu).lt.1000) then 
                        flow(mu)=0.05
                        else if (psi(mu).lt.2500) then
                        flow(mu)=0.06
                        else if (psi(mu).lt.5000) then
                        flow(mu)=0.08
                        else 
                        flow(mu)=0.09
                        endif   
        swx(mu)=swx(mu)+flow(mu)/dlayr(mu)
  509  continue
      DO 510 L=1,MU
        SW(L)=SWX(L)
  510 CONTINUE
      CES=CES+ES

C********** CALCULATE TRANSPIRATION ********************************
      EP=0.  !; write(*,*)istage
      IF (ISTAGE.GT.4) GO TO 520

	EP=EO*kcmax*(min(1.,lai/3))**0.5  !JFM051121
	IF ((EP+ES)>(kcmax*EO)) EP=kcmax*EO-ES
	
	if (adjustl(typsol)=="latosol") then !JFM060413
	!write(*,*) kcmax
	EP=EO*(kcmax*1.03)*(1-EXP(-kexrn*LAI)) ! jfm050902
	IF ((EP+ES)>(kcmax*EO)) EP=kcmax*EO-ES
	end if

	tmp=EP !JFM050905  Transpiration maxi mm
	etm=EP+ES !JFM050905 EvapoTranspiration maxi mm

      GO TO 600
  520 ET=ES
      CET=CET+ET
      TSW=0.
      DO 530 L=1,NLAYR
	TSW=TSW+SW(L)*DLAYR(L)
        flux(l)=0.		
  530 CONTINUE
      PESW=TSW-TLL
      RETURN
C********** CALCULATE ROOT GROWTH AND DEPTH ************************
  600 IF (ratoon>0) GO TO 641 !jfmodif
c---- calculate root depth
	RTDEP=DEPROO
        tempm=(tempmx+tempmn)/2
	if (tempm.gt.pBASET) deproo=deproo+tempm*vitenr
	if (deproo.gt.pRDPMX) deproo=pRDPMX

	dep=0.
        DO 640 L=1,NLAYR
	if (dep.gt.deproo) go to 700 
c---- calculate root length density. RLV(L)  ----
		DRLV(L)=1.6e-2*pDSTMX*(DEPROO-RTDEP) !JFM<051013 Equation initiale de watbal
	if (adjustl(typsol)=="latosol") then  !JFM060413
		DRLV(L)=0.8e-2*pDSTMX*(DEPROO-RTDEP) !JFM051013   Calage sur Données Patricia Laclau
	end if
c---- effect of water stress on growth of root lenth density 
        SWDF=1.
         ESW(L)=SAT(L)-LL(L)
         IF (SW(L)-LL(L).LT.0.25*ESW(L)) SWDF=4.*(SW(L)-LL(L))/
     1    ESW(L)
        IF (SWDF.LT.0.) SWDF=0.
        IF (SWDF.GT.1.) SWDF=1.
        DRLV(L)=SWDF*DRLV(L)
	RLV(L)=RLV(L)+DRLV(L)
	if (rlv(l).gt.pDSTMX) rlv(l)=pDSTMX
	dep=dep+dlayr(l)
  640   CONTINUE
  641   continue

C********** CALCULATE WATER UPTAKE AND SOIL DEFICIT FACTORS ********
  700  IF (EP.EQ.0.) GO TO 750
      EP1=EP*0.1
      TRWU=0.

      DO 710 L=1,NLAYR
	IF (RLV(L).EQ.0.0) GO TO 720
	RWU(L)=2.67E-3*EXP(62.*(SW(L)-LL(L)))/(6.68-ALOG(RLV(L)))	
c RWUMX (max transpiration by unit of root) set to 0.03 here    
	 IF (RWU(L).GT.RWUMX) RWU(L)=RWUMX 
 
       IF (SW(L).LT.LL(L)) RWU(L)=0.
       RWU(L)=RWU(L)*DLAYR(L)*RLV(L)
       if (rwu(l).gt.(dlayr(l)*(sw(l)-ll(l)))) 
     1  rwu(l)=(dlayr(l)*(sw(l)-ll(l)))
	 IF (RWU(L)<0.) RWU(L)=0  !JFM030527
       TRWU=TRWU+RWU(L)
  710 CONTINUE
  720 WUF=1
	IF (EP1.Lt.TRWU) WUF=EP1/TRWU
      TSW=0.
      DO 730 L=1,NLAYR
	RWU(L)=RWU(L)*WUF
	SW(L)=SW(L)-RWU(L)/DLAYR(L)
	TSW=TSW+SW(L)*DLAYR(L)	
  730 CONTINUE
      PESW=TSW-TLL
      SWDF2=1.
      IF (TRWU/EP1.LT.1.5) SWDF2=0.67*TRWU/EP1
  740 SWDF1=1.
      IF (EP1.LT.TRWU) GO TO 750
      SWDF1=TRWU/EP1
      EP=TRWU*10.
	
  750 ET=ES+EP  	
	etr=et 
	trp=ep 
C ** water balance checking
      stockf=0.
      do 809 l=1,nlayr
         stockf=stockf+dlayr(l)*sw(l)*10
 809      continue
          error=stockf-stocki+(et+drain-pinf)
		if (abs(error)>0.1) then  
		drain=drain-error
          error=stockf-stocki+(et+drain-pinf)
	end if
C ** water stress and capacities calculations
		dep=0;stockfr=0;llwr=0;dulwr=0;satwr=0; stock=0 !jfm0319
		do l=1,nlayr
			  stock=stock+dlayr(l)*sw(l)*10          
		end do	    
		do l=1,nlayr
			if (deproo>dep+dlayr(l)) then
			  stockfr=stockfr+dlayr(l)*sw(l)*10
			  llwr=llwr+ll(l)*dlayr(l)*10
			  satwr=satwr+sat(l)*dlayr(l)*10
			  dulwr=dulwr+dul(l)*dlayr(l)*10
			 else
			  stockfr=stockfr+dlayr(l)*sw(l)*10*(deproo-dep)/dlayr(l)
			  llwr=llwr+ll(l)*dlayr(l)*10*(deproo-dep)/dlayr(l)
			  satwr=satwr+sat(l)*dlayr(l)*10*(deproo-dep)/dlayr(l)
			  dulwr=dulwr+dul(l)*dlayr(l)*10*(deproo-dep)/dlayr(l)
			end if
			dep=dep+dlayr(l)
			if (deproo<=dep) exit 
		end do
		stock=stock-llw; stockfr=stockfr-llwr
		if (stock<0.) stock=0.      !Total available water capacity
		if (stockfr<0.) stockfr=0.  !Root available water capacity
		swdef=1. 
		if(stockfr<(dulwr-llwr)) then   
			swdef=stockfr/(dulwr-llwr)
		end if
		if (swdef>1) swdef=1; if (swdef<0) swdef=0
      END






c******************* FIN WATBAL ************************************************


c**************  Initiation sol 1er jour
	SUBROUTINE wbinit1  !JFM051216
		INCLUDE 'cercane1.inc'  !JFM000301
	real rutot,pevap,Stockini1, deprooini,tempm
	if (nage==1) then !initialisation
		stock=0.; depsol=0; deproo=0; rtdep=0.; dulw=0.; llw=0	 
		do l=1,nlayr  !calcul des stocks et profondeurs totales
			stock=stock+dlayr(l)*sw(l)*10.
			llw=llw+dlayr(l)*ll(l)*10.
			dulw=dulw+dlayr(l)*dul(l)*10.				
			depsol=depsol+dlayr(l)				
		!write(*,*) sw(l)
		end do	
		llw=ceiling(1000*llw)/1000
		dulw=ceiling(1000*dulw)/1000		
		stock=ceiling(1000*(stock-llw))/1000   !Le stock total devient stock utile en mm
		
		if (pRDPMX>depsol) pRDPMX=depsol
		profsurf=dlayr(1)
		pevap=1.-p0		
		CapaRevap=0.5*dlayr(1)*ll(1)*10. !JFM051121
		CapaRuSurf=ceiling(1000.*profsurf*(dul(1)-ll(1))*10.)/1000.  !JFM051121		
		CapaRfe=pevap*CapaRuSurf !JFM051216 different de Ecotrop Caparevap+caparusurf!!!!!!!
		CapaRde=CapaRuSurf-CapaRfe		
		Stockini1=(sw(1)-0.5*ll(1))*profsurf*10		
		ValRSurf=min(Stockini1,CapaRevap+CapaRde)  !JFM051216		
		ValRde=max(0.,ValRSurf-CapaRevap)
		!ValRde=ceiling(1000.*ValRde)/1000.		
		ValRfe=max(0.,Stockini1-(CapaRevap+CapaRde))
		!ValRfe=ceiling(1000.*ValRfe)/1000.  !JFM051121
		StockSurface=ValRde+ValRfe
	    stockfr=0.; stockft=0.
		if (ratoon>0) deproo=pRDPMX
		if (ratoon==0) deproo=0.7   !JFM130603
		stockfr=stock*deproo/depsol
		llwr=llw*deproo/depsol; dulwr=dulw*deproo/depsol
		satwr=satw*deproo/depsol
		Caparac=dulwr-llwr; Capatot=dulw-llw
		!write(*,*) "res", nage,valrsurf,valrfe,valrde,stockfr
		
		Ltr=1
	end if
	if ((ratoon==0).and.(nage>1)) then   !JFM130603
		deprooini=deproo; tempm=(tempmx+tempmn)/2
 	    if (tempm.gt.pBASET) deproo=deproo+tempm*vitenr  !type ceres JFM130604
		if (deproo.gt.pRDPMX) deproo=pRDPMX
		if (deproo<=profsurf) then 
			stockfr=(valrfe+valrde)*deproo/profsurf
	     else
	        stockfr= stockfr+stock*(deproo-deprooini)/depsol
		end if
		llwr=llw*deproo/depsol; dulwr=dulw*deproo/depsol
		satwr=satw*deproo/depsol
		Caparac=dulwr-llwr; Capatot=dulw-llw
	end if
	END SUBROUTINE wbinit1
c=========================================================================


c=========================================================================
	SUBROUTINE wbkce(mode,kce,lai)  
	real lai            !IN
	character mode*10   !IN
	real kce            !OUT
      SELECT CASE (trim(adjustl(mode)))  
		case('ceres1----')
			kce=1.-0.43*lai
			IF (lai>1.) kce=EXP(-0.4*lai)/1.1
	end select
	END SUBROUTINE wbkce
c=========================================================================

c=========================================================================
	SUBROUTINE wbkcp(mode,kcp,kcmax,lai)  
	real lai,kcmax       !IN
	character mode*10    !IN
	real kcp             !OUT
      SELECT CASE (trim(adjustl(mode)))
		case('FAOLAI1---')  !Ecotrop
			kcp=kcmax*(min(1.,lai/3))**0.5  !
	end select
	END SUBROUTINE wbkcp
c=========================================================================

c=========================================================================
	SUBROUTINE wbevapotfao(eeq,kce,eos)
	real eeq,kce       !INVAR
	real eos           !OUTVAR
	eos=kce*eeq 
	END SUBROUTINE wbevapotfao
c=========================================================================

c=========================================================================
	SUBROUTINE wbevapfao(eos,valrsurf,valrde,valrfe,stocksurface,
	1           stockfr,stock,caparevap,caparde,es,deproo,profsurf)
	!***** FAO/Ecotrop EVAPORATION !JFM051217
	real eos,caparevap,caparde                   !IN
	real stockfr,stock,valrsurf,valrde,valrfe    !INOUT
	real es,stocksurface                         !OUT
	real evap1,evap2,evaprde,valrdei,valrfei,defr     !TEMP
	valrfei=valrfe; valrdei=valrde; defr=0
	evaprde=0; es=0
	if (valrfe>=eos) then !JFM051121
		evap1=eos; evap2=0. 
		else				
		evap1=valrfe		
		evap2=(eos-evap1)*(valrsurf)/(caparevap+caparde)  
	    evap2=max(0.,evap2)
	end if  
	es=evap1+evap2  
	valrfe=valrfe-evap1  
	valrsurf=valrsurf-evap2
	valrde=max(0.,valrde-evap2)
	stocksurface=valrfe+valrde	 
	defr=(valrfe+valrde)-(valrfei+valrdei)*min(1.,(deproo/profsurf))  
	stockfr=stockfr+defr	   
	stock=stock+(valrfe+valrde)-(valrfei+valrdei)	
	END SUBROUTINE wbevapfao
c=========================================================================

c=========================================================================
	SUBROUTINE wbcstranspfao(eeq,es,kcp,stockfr,caparac,p0,cstr)
!	!***** FAO/ECOTROP  STRESS TRANSPIRATION   CSTR  
	real eeq,es,kcp,stockfr,caparac,p0  !IN
	real kctot,swdef,pfact                 !TEMP
	real cstr                              !OUT
	kctot=es/eeq+kcp !kce*kr+kcp
	swdef=stockfr/caparac
	pfact=(1-p0)+0.04*(5.-kctot*eeq)  !pfact=parpfact+0.04*(5.-kctot*et0)
	pfact=max(0.1,pfact)
	pfact=min(0.8,pfact)
	cstr=max(0.,min(swdef/(1-pfact),1.)) !cstr=max(0.,min(ftsw/(1-pfact),1.))
	END SUBROUTINE wbcstranspfao
c=========================================================================

c=========================================================================
	SUBROUTINE wbtranspfao(cstr,eeq,es,kcp,deproo,profsurf,stockfr,
	1     stock,stocksurface,valrde,valrfe,valrsurf,tmp,trp,etr,etm)
	!***** FAO/ECOTROP  TRANSPIRATION				!JFM051217
	real  cstr,eeq,es,kcp,deproo,profsurf			!IN
	real stockfr,stock,valrde,valrfe,valrsurf		!INOUT
	real stocksurface,tmp,trp,etr,etm				!OUT
	real trsurf										!TEMP
	tmp=kcp*eeq ; trp=tmp*cstr
	etr=es+trp ; etm=es+tmp	
	if (deproo>profsurf) then
	  trsurf=trp*stocksurface/stockfr
	  else
	  trsurf=trp
	end if
	if (trsurf>=valrfe) then
	    valrde=valrde-(trsurf-valrfe)
		valrsurf=valrsurf-(trsurf-valrfe); valrfe=0.
		else
		valrfe=valrfe-trsurf
	end if
	stockfr=stockfr-trp ; stock=stock-trp		

	stocksurface=valrde+valrfe
	END SUBROUTINE wbtranspfao
c=========================================================================

c=========================================================================
	SUBROUTINE wbtransfert1(precip,irdose,seuilruis,pcruis,caparevap,
	1	valrsurf,caparfe,valrfe,caparde,valrde,deproo,profsurf,runoff,
     2    stocksurface,caparac,stockfr,capatot,stock,drain)
	!***** FAO/ECOTROP  Remplissage réservoirs     
	real precip,irdose,seuilruis,pcruis,caparevap,caparfe,caparde  !IN
	real deproo,profsurf,caparac,capatot                           !IN
	real  runoff,drain,stocksurface                                !OUT
	real  valrsurf,valrfe,valrde,stockfr,stock                     !INOUT
	real  eaudispo,eautemp,eautransp,stockrprof,caparprof          !TEMP
  
	drain=0.; runoff=0.; eaudispo=0.; eautransp=0.; eautemp=0
	!** Streaming
		if (precip>seuilruis) runoff=(precip-seuilruis)*pcruis/100  
		eaudispo=max(0.,precip+irdose-runoff)
	!** fill revap if valrsurf<caparevap	
		eautemp=min(eaudispo,max(0.,caparevap-valrsurf))
		eaudispo=max(0.,eaudispo-eautemp)  
		valrsurf=valrsurf+eautemp
		eautransp=eaudispo
	!** fill rfe
		eautemp=min(eaudispo,caparfe-valrfe)
		eaudispo=max(0.,eaudispo-eautemp)   
		valrfe=valrfe+eautemp
	!** fill rde
		eautemp=min(eaudispo,caparde-valrde)
		eaudispo=max(0.,eaudispo-eautemp)  
		valrde=valrde+eautemp
	    valrsurf=valrsurf+eautemp
		stocksurface=valrfe+valrde
	!** update stockfr
		if (deproo<=profsurf) then
			stockfr=stocksurface*deproo/profsurf
		else
			eautemp=min(eautransp,caparac-stockfr)
			eaudispo=max(0.,eautransp-eautemp)   
			stockfr=stockfr+eautemp
		end if
	!** update stock
	    drain=0.
		eautemp=min(eautransp,capatot-stock)
		drain=max(0.,eautransp-eautemp)  
		stock=stock+eautemp
	END SUBROUTINE wbtransfert1
c=========================================================================

c=========================================================================
	SUBROUTINE wbswdef(swdef,swdefsat,stockfr,satwr,dulwr,llwr)
	!***** calculate water stress on soil
	real stockfr,satwr,dulwr,llwr !IN
	real swdef,swdefsat           !OUT
		swdef=1.;swdefsat=1.
		if(stockfr<(dulwr-llwr)) then !fillig ratio on available water capacity
			swdef=(stockfr)/(dulwr-llwr)
		end if
		if(stockfr>(dulwr-llwr)) then !fillig ratio on saturated zone
			swdefsat=(stockfr-dulwr)/(satwr-dulwr)
		end if
	END SUBROUTINE wbswdef


