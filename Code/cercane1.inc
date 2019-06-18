
c  DECLARATIONS	TYPES          
		  TYPE  paramvar
		    CHARACTER*15  nom
		    REAL          min,max,ecart,val
		    INTEGER       nbre,indice
          END TYPE





c ESSAI ***************************************
		TYPE irrig		!type irrigation  !JFMalloc061116
			integer irryr,irrday,kindirr !année,jour calendaire
			real effirr,irrdose !dose apportée en mm
		END TYPE
		TYPE ferti		!type fertilisation azotée !JFMalloc061116
			integer fertNyr,fertNday
			real fertNdose,dlayN	
		END TYPE
	    TYPE obssol		!type observation du sol !JFMalloc061116
			integer swoyr,swoday !année,jour calendaire
			real praco,dlro(6),swo(11) !prof racin,dens long racin, humidité volumé
		END TYPE
		TYPE trtm   !type traitement jfm020904    !JFMalloc061116
			character name*30,cultivar*15,nomsol*20,iname*25
			character datedeb*10,datefin*10
			character*20 meteo(3),meteosta !noms fichiers meteo jfm20170524
			character*20 pluvio(3) !noms fichiers pluvio
			character*60 meteostr !chaine meteo			
			integer ratoon,nbbgi,metabdeb
			real lati,longi,alti, profrac,ecart,vsolini		
		END TYPE
		TYPE meteoty   !type meteo
			character*20 meteosta
			character*10 meteodat			
			real mprecip,meeq,mtempmx,mtempmn,msolrad,mux,mun,mvt,mtempmm
		END TYPE

		TYPE totot		!type fertilisation azotée
			integer toti
			real totr	
		END TYPE

		integer,parameter :: metabmax=2000, trtmax=10000, irrmax=1000 !JFM140309 10000 <= 100000
		integer,parameter :: fertNmax=1000, obsolmax=1000		
		TYPE (trtm),   DIMENSION (trtmax)     :: trt    !JFM061117
		TYPE(meteoty), DIMENSION (metabmax) :: metab  !tableau donnees meteo
		TYPE(irrig),   DIMENSION (irrmax)     :: irr    !tableau données irrigation                  !JFMalloc061116
		TYPE(ferti),   DIMENSION (fertNmax)   :: fertN  !tableau fertilisation azotée		      !JFMalloc061116
		TYPE(obssol),  DIMENSION (obsolmax)  :: obsol  !tableau observations sol                 !JFMalloc061116
!		TYPE(trtm), ALLOCATABLE::trt(:) !tableau données traitements !jfm020904     !JFMalloc061116  JFM061117
!		TYPE(ferti), ALLOCATABLE::fertN(:) !tableau fertilisation azotée		      !JFMalloc061116
!		TYPE(obssol), ALLOCATABLE::obsol(:) !tableau observations sol                 !JFMalloc061116
!		TYPE(meteoty), ALLOCATABLE::metab(:) !tableau donnees meteo
		        
		integer notrait,nbtrt
		common /alloc/notrait,nbtrt  


c  DECLARATIONS *****************************************************************
        CHARACTER*50 versionmos
		common/version/ versionmos
c parametrage	--------------------------------
		TYPE (paramvar), DIMENSION (12) :: tabpar
		integer, dimension(12) :: par_nbre,par_indice,nvaripar
		real , dimension(12) :: par_min,par_max,par_ecart,par_val ,varipar
		CHARACTER*15, par_nom(12)		
		REAL(8) elapsed_time, elaps1
		common /manip/elapsed_time, elaps1
		
		integer typeparam,nbpar,ipar,nbrecuit,parac,typecalage
		character varparam*15, finparam*3
		integer nbrobs
		REAL rcalc, robs, srcalc, srobs, sdif2, scalcobs 
		REAL sdifabs,srobs2, srcalc2, rmse, mae, alin, blin, r2
c		real , dimension (11) :: dul,ll
		common /manip1/nbrecuit,parac
		common /manip2/typeparam,typecalage,nbpar,ipar,varparam,finparam
		common /param2/par_nom,par_nbre,par_indice,par_min,par_max,par_ecart
		common /param3/par_val,rcalc, robs, nbrobs, srcalc, srobs
		common /param4/sdif2, scalcobs,sdifabs,srobs2, srcalc2
		common /param5/rmse, mae, alin, blin, r2



c  REPERTOIRES **************************************************************
		character repmeteo*100, repdata*100, represult*100
		common /repertoire/ repmeteo, repdata, represult


c  SITUATION **************************************************************

		real lati,longi,alti
		common /situ1/lati,longi,alti




c  **** OBSERVATIONS  **************************************************
	INTEGER obsutil,isobei
	REAL obnbtigv,onnbtigm,obnbtigfl,obhtvd,obhtvdx
	REAL obnbt,obnbv,obltvdsurf,oblai,obei,obrdcan,obmstigusi,obsucre
	REAL obhumtig,obsucretig,obfibretig,obbrixjtig,puretejtig	
	REAL obmslv,obmslt,obmfa,obmsa
	REAL obhumcan,obbrixcan,obcendrecan,obpcihcan,obmsen
	REAL obfibrecan,obsucrecan,obpuretejcan,obbrixjcan
	REAL obstock,obstockfr
	REAL obndftu,obndfae,obholndftu,obholndfae
	REAL oblignndftu,oblignndfae
	
	common /ob0/obsutil,isobei
	common /ob1/obnbtigv,obnbtigm,obnbtigfl,obhtvd,obhtvdx
	common /ob2/obnbt,obnbv,obltvdsurf,oblai,obei,obrdcan,obmstigusi
	common /ob3/obhumtig,obsucretig,obfibretig,obbrixjtig,puretejtig	
	common /ob4/obmslv,obmslt,obmfa,obmsa,obsucre
	common /ob5/obhumcan,obbrixcan,obcendrecan,obpcihcan,obmsen
	common /ob6/obfibrecan,obsucrecan,obpuretejcan,obbrixjcan
	common /ob7/obstock,obstockfr
	common /ob8/obndftu,obndfae,obholndftu,obholndfae
	common /ob9/oblignndftu,oblignndfae


c **** WATBAL ********************************************************************	  
      real flow(20),flux(20),swx(20),stocki,stockf,error !=> cercane1.inc jfm
      REAL LL,drl1,avth,aaa,alf,DRLV(15),RTDEP,cr,runoff !=> cercane1.inc jfm
c N
			real ndef1, ndef2
			real no3,nh4
      COMMON /CLIMT/ TEMPMN,TEMPMX,tempmm,precip,SOLRAD,TMFAC(8)
      COMMON /DATEC/ MO,ND,IYR,IMO,IDD,JDATE,JDATEX,IDIM(12)
			COMMON /STRESS/iwatbal,inbal
      common /etp/ eeq,psi(50)
c      COMMON /GROTH/ GPSM,GPP,GRORT,PTF,LAI,DM,BIOMAS,PLA,SENLA,
c     1  LFWT,SEEDRV,XPLANT,EMAT,PLAY,PLAMX,XNTI,LWMIN,RTWT,STMWT,
c     2  GRNWT,SWMIN,LN,EARWT,TLNO,XN,RWID,SUMP,IDURP,PLAG,GROSTM,
c     3  CARBO,GBLA(35),SLA(35),EARS,REGM,SWMAX
      COMMON /SOILI/ SALB,U,SWCON,DLAYR(15),DUL(15),ll(15),SW(15),aaa,
     1  SAT(15),DEPMAX,TDUL,NLAYR,SMX,WF(10),WR(11),RWU(15),SWEF,CN2
      COMMON /WATER/ SUMES1,SUMES2,T,Tll,PESW,TSW,CUMDEP,ESW(15),
     1  CSD1,CSD2,SI1(6),SI2(6),ICSDUR,ES,EP,ET,EO,CES,CEP,CET,pinf,
     2  RLV(15),CRAIN,DRAIN,IDRSW,RTDEP,SWDF1,SWDF2,STOCKF
     3  TRWU,RWUMX,TPESW,TSAT,CIRR,sperc,spinf,srunof,set,seo,ses,sep
c variables provenant de cerca      
			COMMON /NDFPG/ NDEF1,NDEF2,ndefpod
      COMMON /NSTEMP/ ST(10),ANG,TMN,AMP,TAV,T0(5),ST0,DD,TA,DT,ALX,Z
      COMMON /NGEN2/ STRAW,SDEP,SCN,ROOT,RCN,OC(11)
      COMMON /NSPOOL/ SNH4(10),SNO3(10),NH4(11),NO3(11),FAC(10),BD(15),
     1  PH(11),ABD,ANO3,ANH4,NUMSL
      COMMON /GROTH/ GDMROO,DMROOT,GDMTOP,DMTOPS,GDMVEG,DMVEGE,dmlf,
     1  DMPODS,DMHUSK,DMGRAI,DMPM2,PAI,PAIMAX,endexp,ddmlf
		  COMMON /NNNN/ CNTOT, CEC(11)


      COMMON /WRITS/ AES,AEP,AET,AEO,ASOLR,ATEMX,ATEMN,ARUNOF,
     1  ADRAIN,APRECP,ASWDF1,ASWDF2,IOUTGR,IOUTWA,JHEAD,KHEAD,
     2  RUNOFF
        common /gamma/ gam,alf
        common /crois_rac/ vitenr,pRDPMX,pDSTMX,deproo,pBASET
	common /water1/ep1
	
		
c     ENTREES general, irrigation,  ferti
	character var*15,param*15,cultivar*15,typsol*20,perso*3,testsub*3
	integer ratoon,stopmod,dstopmod,mstopmod
	real rowspac
	real val,irdose
	real Wver, K2ver
	integer iwatbl,ifreq,iobs,nobsol,nfertN,kindirr
	    common/sol/typsol,Wver, K2ver
		common /obssol/nobsol
		COMMON/GENERAL1/var,param,cultivar,perso,testsub	  
		COMMON/GENERAL2/ratoon,rowspac
		COMMON/GENERAL3/val,iwatbl,ifreq,iobs,stopmod,dstopmod,mstopmod
		common/irrig1/nirr,kindirr,effirr,error,irdose
		common /fertN1/nfertN





c	ENTREES PARAMETRES PLANTES *******************************************
	real aa1,aa2,aa3,aa4,aa5
	real p01,p02,p03,p04,p05,p06,p07,p08,p09,p10
	real ind0,ind1,ind2,ind3,ind4,ind5,ind6
	real ind7,ind8,ind9,ind10,ind11,ind12
	real sthydcroi, sthydbio, psilb
	real pblfin,pbldif,pblkdec,slafin,sladif,slakdec
	common /partempo1/p01,p02,p03,p04,p05,p06,p07,p08,p09,p10
	
	
	!** Apparition limbes	 		
	real bltb,blk,blkdec,blvmax
	common /blades01/bltb,blk,blkdec,blvmax
	!** Taille limbes
	real blsurfk,blsurfmax	 	
	common /plbladsurf1/blsurfk,blsurfmax	


c   ******* Apparition senescence tiges
	!** Tallage global
	real taldebtt,talpeaktt,taltb,talpeakval,talfinval !talcas1	 
	!real talto1,talto2,taltm,taltp,taldeb,talkapp,taldebsen !Tallage
	!real talksen,talfin !pltal
	common /pltalglob1/talto1,talto2,taltm,taltp,taldeb,talkapp  !talcas1
	!common /pltalglob2/taldebsen,talksen,talfin
	common /pltalglob0/taldebtt,talpeaktt,taltb,talpeakval,talfinval !talcas1
	!** Germination
	integer nbbgi
	!real gertb,gerto,gertm,gerdeb,gerfin,gerkapp
	!common /talgerm1/gertb,gerto,gertm,gerdeb,gerfin,gerkapp,nbbgi
	!** Tallage
	!real tal2deb,tal2kapp,tal2debsen,tal2ksen,tal2fin
	!common /pltal1/tal2deb,tal2kapp,tal2debsen,tal2ksen,tal2fin

	!** Elongation
	real htvdtb,htvdto,htvdtm,htvdcan,htvdgro1,htvdgro2 !plh1
	common /plh1/htvdtb,htvdto,htvdtm,htvdcan, htvdgro1,htvdgro2


	!** humidite
	real humtigk1,humtigk2
	real humaerdeb, humaerfin, humaertb, humaerdec
	common /plhumtig/humtigk1,humtigk2
	common /plhumaer/ humaerdeb, humaerfin, humaertb, humaerdec
	
	!** feuilles lai ke conv part
	real dsenftig,ksenftig,ksurf,maxsurf,phyll,phyltbase  !plai1
	real laitb,laicroi,laiwksen
	!real laideb,laitb,laidebsen,laitksen,laiwksen,laicroi
	real tbase,debei,keidj,eimax,sdjei,eideb,ke !plai3,plei1
	real ruetk,ruetopt,ebmax,epb,tresent,trescro,ruemax !plconv
	real toptrue,ktrue, prue, conv_a, conv_b, ktemp !plconv
	real debmsa,k_maer,msadeb !plpart
	real pstrutb,pstrutcroi,pstrufin,pstrudec
	real ptigdeb,ptigfin,ptigdec,ptig   !plpart
	real profracvar,ktstru,profrac,vsolini,kcini,kexrn,kcmax
	!* Ligno celluloses
	real tuaera,tuaerb,tuaerk, ndftua,ndftub,ndftuk
	real celtua,celtub,celtuk, hceltua,hceltub,hceltuk
	real ligntua,ligntub,ligntuk
	common /cal/aa1,aa2,aa3,aa4,aa5
	common /inda/ind0,ind1,ind2,ind3,ind4,ind5,ind6
	common /stresshydr/sthydcroi,sthydbio,psilb
	common /indb/ind7,ind8,ind9,ind10,ind11,ind12
	common /plai0/laitb,laicroi,laiwksen
	common /plai1/dsenftig,ksenftig,ksurf,maxsurf,phyll,phyltbase
	!common /plai2/laideb,laitb,laidebsen,laitksen,laiwksen,laicroi
	common /plei1/tbase,debei,keidj,eimax,sdjei,eideb,ke
	common /plconv1/ruetk,ruetopt,ebmax,epb,tresent,trescro,ruemax
	common /plconv1/toptrue,ktrue, prue, conv_a, conv_b, ktemp
	common /plpartrac/debmsa,k_maer,msadeb
	common /plpartblad/pblfin,pbldif,pblkdec,slafin,sladif,slakdec
	common /plparttig/ptigdeb,ptigfin,ptigdec,ptig
	common /plpartsuc/pstrutb,pstrutcroi,pstrufin,pstrudec
	common /rac1/profracvar,ktstru,profrac,vsolini,kcini,kexrn,kcmax
	common /clc1/tuaera,tuaerb,tuaerk, ndftua,ndftub,ndftuk
	common /clc2/celtua,celtub,celtuk, hceltua,hceltub,hceltuk
	common /clc3/ligntua,ligntub,ligntuk


c		 common /fertiN2/fertNyr(20),fertNday(20),qNitr(20),dlayN(20)

c VARIABLES CALCULEES *******************************************************
	    
 	   !** Sorties Astro + ray  JFM120513
		real astlatrad,astdr,astdecl,astomeg,astrgo,astrgx,astdurjour,astdurjourp
	    common /astro1/ astlatrad,astdr,astdecl,astomeg,astrgo,astrgx
		common /astro2/ astdurjour,astdurjourp

		!** Sorties floraison  JFM160204
		integer stflor, niflor
		real iflor, pflor
		common /outflor/iflor,niflor, pflor

		integer nage,jcal,istage, njtot, jat,daytoend 
		
		real depsol,stockfr, llwr,dulwr, satwr,stockfdeb
		real llwt,dulwt, satwt,stock,swdefsat
		real stockm,stockfrm
	    real stmo,srg,spr,setp,pucum,sirr,setm,setmt,setr,setrt,par,etr,etm
		real tmp,stmp,trp,strp,swdef,llw,dulw,satw,p0,setpt,stmpt,strpt
	    real pari,spar,spari,sdj,parim,sparim,sdjlai,sdj1,srgo,srgx
		real blnbtot,blnbv,blsdj,bldj,bltvdsurf,bltvdlong,bltvdlarg,bltigsurf
	    real htvd,lai,glait,glaibiom,ei,ktemprue,htvddj,htvdsdj,dhtvd
	    real gdmtot,dmtot,dmbaer,dmrac,gdmsen,dmsen,dmbaert
		real gdmbaer,sdmbl,gdmrac,gdmsug,gdmstm,gdmndftu	    
		real mst,dmstst,watstm
	    real dmstru,gdmbl,dmtbl,slam
	    real dmstm,dmfibstm,dmredstm,yldcan,dmsug
	    real br_stm,pur_stm,fib_stm,sug_stm,hum_stm,hum_aer, humaersdd
		real ndftu,ndfae,holndftu,holndfae,lignndftu,lignndfae
		real dmndftu,dmndfae,dmceltu,dmcelae,dmhceltu,dmhcelae,dmligntu,dmlignae
	    real cslai,csh,cske,csrue,aNitr,dNitr,ktot
		real sswdf1, sswdf2, sswdef,sdrain,srunoff,sevap
		real sdjblt,sdjhtvdt,sswdeft,sswdf1t,sswdf2t
	    real stmot,stnt,stxt,spart,sparit, sdraint,srunofft
		real stmofm,sparifm,sprfm,swdeffm,swdf1fm,swdf2fm
		real stnfm,stxfm
		real kce,kcp,profsurf,llSurf,seuilruis,pcruis,Ltr          !FAO/Ecotrop
		real StockSurface,ValRSurf,ValRde,ValRfe,cstr,eos          !FAO/Ecotrop
		real CapaRfe,CapaRde,CapaRevap,CapaRuSurf,Caparac,Capatot  !FAO/Ecotrop
		real cum1,cum2,cum3,cum4,cum5                              ! Variables virtuelles
 
 	   !** Sorties Tallage
        real nbtigv, nbtigt,nbtigvl,nbtigtl,nbtigpvl
		real nbtalvl, nbtaltl, talsdj,taldj
		integer pictal,sttal,stlev
        common /outtal1/nbtigv,nbtigt,nbtigvl,nbtigtl,nbtigptl,nbtigpvl
		common /outtal2/nbtalvl,nbtaltl,talsdj,taldj
		common /outtal3/pictal,sttal,stlev

		common /ctmp1/nage,istage,jat,njtot,daytoend, stflor
	    common /cclim1/stmo,srg,spr,setp,pucum,sirr,setm,setmt,setr,setrt,par
		common /cclim2/srgo,srgx,swdef,llw,dulw,satw, p0                     !variables virtuelles
		common /cclim3/cum1,cum2,cum3,cum4,cum5
		common /ccindt1/sdjblt,sdjhtvdt,sswdeft,sswdf1t,sswdf2t
	    common /ccindt2/stmot,stnt,stxt,spart,sparit
		common /ccindfin1/stmofm,sparifm,sprfm,swdeffm,swdf1fm,swdf2fm
		common /ccindfin2/stnfm,stxfm
		
		common /outgro0/pari,spar,spari,sdj,parim,sparim,sdjlai
		common /outgro1/blnbtot,blnbv,blsdj,bldj,bltvdsurf,bltigsurf,bltvdlong,bltvdlarg
	    common /outgro2/htvd,lai,glait,glaibiom,ei,sdj1,htvddj,htvdsdj
        common /outgro3/ktemprue,dhtvd
	    common /outgro4/gdmtot,dmtot,dmbaer,dmrac,gdmsen,dmsen,dmbaert
	    common /outgro5/mst,dmstst,watstm
	    common /outgro6/dmstru,gdmbl,dmtbl,slam
	    common /outgro7/dmstm,dmfibstm,dmredstm,yldcan,dmsug
	    common /outgro8/gdmbaer,sdmbl,gdmrac,gdmsug,gdmstm,gdmndftu
        common /outgro9/br_stm,pur_stm,fib_stm,sug_stm,hum_stm,hum_aer,humaersdd
	    common /outgro11/cslai,csh,cske,csrue
		common /outgro12/sswdf1,sswdf2,sswdef
		common /outgro13/sdrain,srunoff,sevap,sdraint,srunofft
		common /outgro14/ndftu,ndfae,holndftu,holndfae
		common /outgro15/lignndftu,lignndfae
		common /outgro16/dmndftu,dmndfae,dmceltu,dmcelae
		common /outgro17/dmhceltu,dmhcelae,dmligntu,dmlignae
	    common /fertiN2/aNitr,dNitr,ktot
		common /wb11/depsol,stockfr,llwr,dulwr,satwr,stockf,stockfdeb
		common /wb12/llwt,dulwt,satwt,stock,etr,etm,swdefsat
		common /wb13/tmp,stmp,trp,strp,stockm,stockfrm,setpt,stmpt,strpt
		common /wb21/kce,kcp,profsurf,llSurf,seuilruis,pcruis,Ltr
		common /wb22/StockSurface,ValRSurf,ValRde,ValRfe,cstr,eos
		common /wb23/CapaRfe,CapaRde,CapaRevap,CapaRuSurf,Caparac,Capatot		

c ;iwatbl;0-1;0 pas de stress hydrique (swdf1=1et swdf2=1)
c ;ifreq; pas de temps pour les sorties si 0 sortie uniquement aux dates d'observations


c mosicas
	  !real laideb,laimax,klaidj,ksenlai,lailim,ke
	  !real nbtigv,nbtigtot,htvd,ktemprue,gdmtot,gdmbaer
	  !real gdmroot,gdmst
