
c*************************************************
c        CALL SUBROUTINES (from motsimul)
c**************************************************
c**************************************************



	SUBROUTINE appelmodul(s3)
	INCLUDE 'cercane1.inc'
	CHARACTER s3*20,s4*20
	real xx
	character*20 mo1,mo2,mo3,mo4,mo5,mo6	
	s4=s3



	SELECT CASE (adjustl(s3))
	 CASE ('MOSICASCERES','MOSICASFAO56','MOSICASFAO56b',
	1       'MOSICASCERESb')


c********************** WATER BALANCE MODULE  *******************************

	  !*** Call of water balances
	  	SELECT CASE (adjustl(s4))
	      CASE ('MOSICASCERES','MOSICASCERESb')	! ceres water balance module  		
		   kcmax=1.1  !JFM070529
		   call watbal
	      CASE ("MOSICASFAO56","MOSICASFAO56b")	! FAO56 water balance module
		   kcmax=1.0  !JFM070529
		   call wbinit1
		   call wbtransfert1(precip,irdose,seuilruis,pcruis,caparevap,
	1	   valrsurf,caparfe,valrfe,caparde,valrde,deproo,profsurf,
     2       runoff,stocksurface,caparac,stockfr,capatot,stock,drain)
		   call wbkce('ceres1----',kce,lai)
		   call wbkcp('FAOLAI1---',kcp,kcmax,lai)
		   call wbevapotfao(eeq,kce,eos)
		   call wbevapfao(eos,valrsurf,valrde,valrfe,stocksurface,
	1           stockfr,stock,caparevap,caparde,es,deproo,profsurf)
     		   call wbcstranspfao(eeq,es,kcp,stockfr,caparac,p0,cstr)
		   call wbtranspfao(cstr,eeq,es,kcp,deproo,profsurf,stockfr,
	1       stock,stocksurface,valrde,valrfe,valrsurf,tmp,trp,etr,etm)
		   call wbswdef(swdef,swdefsat,stockfr,satwr,dulwr,llwr)
		END SELECT
	  !*** End of water balances calls

	  !*** Calls of Water stress calculations
		call iwatnbal(iwatbal,inbal,swdf1,swdf2,swdef,ndef1,ndef2)
		SELECT CASE (adjustl(s4))
		CASE ('MOSICASCERESb','MOSICASFAO56b')
		call strhydrgen(1,tmp,trp,swdef,swdf1,swdf2,psilb) !type ceres with slabers
		END SELECT
		swdf1=swdf1**sthydbio
		swdf2=swdf2**sthydcroi  !jfm170308 sthydbio avant
	  !*** End of Water stress calculations


c********************** GROWTH MODULE  ***************************************
		
c********** Global temperature and radiation indicators
			tmo=(tempmn+tempmx)*0.5; stmo = stmo + tmo
			par=0.5*solrad;  spar=spar+par;	

c********** Thermal time from beginning (plantation or previous harvest)
		sdj=sdj+max(0.,(tempmn+tempmx)/2-bltb)    !jfm20170225

c********** Flowering
		call flor0 (nage,pflor,stmo,astdurjour,astdurjourp,stflor,iflor,
	1		niflor)

c********** Stalks appearance and senescence
		call talcas1(tempmn,tempmx,talsdj, nbtigv,ratoon, taltb, !jfm20150515
	1           taldebtt,talpeaktt,talpeakval,talfinval) 

c********** Blades appearance, senescence and area growth  (primary stalks)
		call bladesnb0(tempmn,tempmx,bltb,blk,blkdec,
	1			blvmax,blsurfk,blsurfmax,nbtigv,blnbtot,blnbv,
     2            blsdj,bltvdsurf)  !JFM170418)

c********** Stalks elongation (primary stalks)
		call elong2(tempmn,tempmx,htvdtb,htvdto,htvdtm,htvdcan,
	1		htvdgro1,htvdgro2, swdf2, nbtigv,htvd,dhtvd,0)

c********** Calculation of Green leaf area index and interception
		call laiglob0(tempmn,tempmx,laitb,laicroi,
	1		   laiwksen,iflor,swdf2,sdjlai,lai,nbtigv)
		call intercep2(isobei,obei,ke,lai,solrad,ei,par,pari)  !jfm121219

c********** Convertion into above ground dry mass (convertion and partition)
		call convert1(ruemax,ruetk,ruetopt,tempmn,tempmx,pari,
	1				solrad,astrgx,swdf1,gmst,mst,p01,p04)
		call partmst1(ratoon,gmst,stmo, mst,
	1			gdmrac,dmrac,gdmbaer,dmbaer)

c********** Partition of above ground into green blades dry mass and sla calculation
!		call blades1(pblfin,pbldif,pblkdec,slafin,sladif,slakdec,
!	1					dmbaer,gdmbaer,dmtbl,gdmbl,slam)

c********** Partition of above ground into millable stalks dry mass
			call parttige1 (ptigdeb,ptigfin,ptigdec,dmbaer,gdmbaer,
	1						jat,gdmstm,dmstm)

c********** Partition of millable stalks into sugar
			call partsucre1(tempmn, tempmx,gdmstm,dmstm,swdf1,swdf2,
	1                pstrutb,pstrutcroi,pstrufin,pstrudec,gdmsug,dmsug)

c********** Yield and stalk water content calculation
			call rdtsucre2(dmstm,hum_stm,dmsug,yldcan,sug_stm,jat)

c********** Above ground water content calculation
			call humcan1(humaerdeb, humaerfin, humaertb, humaerdec,
	1			hum_aer,tempmx,tempmn,humaersdd,dmbaer)

c********************** INDICATORS  ***************************************
 		call indic2 (jat,bldj,sdjblt,htvddj,sdjhtvdt,swdef,sswdeft,
	1	swdf1,sswdf1t,swdf2,sswdf2t,tmo,stmot,par,spart,pari,sparit,
     2    etm,setmt,etr,setrt)		
! 		call indic2 (jat,bldj,sdjblt,htvddj,sdjhtvdt,swdef,sswdeft,
!	1	swdf1,sswdf1t,swdf2,sswdf2t,tmo,stmot,par,spart,pari,sparit,
 !    2    etm,setmt,etr,setrt)		

		call indic1 (swdef,sswdef,swdf1,sswdf1,swdf2,sswdf2,
	1	pari,spari,astrgo,astrgx,srgo,srgx)
!		call indic1 (swdef,sswdef,swdf1,sswdf1,swdf2,sswdf2,
!	1	tmo,stmo,par,spar,pari,spari,astrgo,astrgx,srgo,srgx)		

	    call indic3 (jat,etm,setm,etr,setr,tmp,stmp,trp,strp,runoff,
	1	srunoff,es,sevap,drain,sdrain,sdraint,srunofft,eeq,setpt,
     2    stmpt,strpt)
!		call indic3 (jat,etm,setm,etr,setr,tmp,stmp,trp,strp,runoff,
!	1	srunoff,es,sevap,drain,sdrain,sdraint,srunofft,eeq,setpt,
!    2    stmpt,strpt)



	END SELECT
	END SUBROUTINE

