!    ROUTINES USING PARAMETERS (SOIL and PLANT)
 !	MAJ_INDICES_PARAM : par_indice(ipar) calculation
 !	MAJ_VAL_PARAM : parameters calculation according to par_nom(ipar) and par_val(ipar)
 !	VER_FIN_PARAM : Determines if all parameters combinations have been used => finparam = "OUI"
 

	
	SUBROUTINE MAJ_INDICES_PARAM
	! par_indice(ipar) calculation
        INCLUDE 'cercane1.inc' !JFM000301
		  If (typeparam==1) Then !MAJ des indices des paramètres
	        DO ipar = 1 , nbpar
                If (par_indice(ipar) < par_nbre(ipar)) Then
                    par_indice(ipar) = par_indice(ipar) + 1
                    Exit 
                  Else
                    par_indice(ipar) = 1
                End If
              END DO	        
	       ELSE	        
	        finparam = "OUI"
	    END IF 
	write(*,*)  par_indice(1),par_indice(2)
	END SUBROUTINE MAJ_INDICES_PARAM






      SUBROUTINE MAJ_VAL_PARAM
	! parameters calculation according to par_nom(ipar) and par_val(ipar)
        INCLUDE 'cercane1.inc' !JFM000301
              DO ipar = 1 , nbpar
                SELECT CASE   (adjustl(par_nom(ipar)))
                   CASE ('aa1')
	                aa1=par_val(ipar)
                   CASE ('aa2')
	                aa2=par_val(ipar)
                   CASE ('aa3')
	                aa3=par_val(ipar)
                   CASE ('aa4')
	                aa4=par_val(ipar)
                   CASE ('blk')
	                blk=par_val(ipar)
                   CASE ('blkdec')
	                blkdec=par_val(ipar)
                   CASE ('blsurfk')
	                blsurfk=par_val(ipar)
                   CASE ('blsurfmax')
	                blsurfmax=par_val(ipar)
                   CASE ('bltb')
	                bltb=par_val(ipar)
				 Case ("blvmax")
                      blvmax = par_val(ipar) !Print *,"bhumtig", bhumtig
                   CASE ('debH20')
	                debH20=par_val(ipar)
                   CASE ('deblai')
	                deblai=par_val(ipar)
                   CASE ('finH20')
	                finH20=par_val(ipar)
                   Case ("gertb")
                    gertb = par_val(ipar) !Print *,"gertb", gertb
                   Case ("gertm")
                    gertm = par_val(ipar) !Print *,"gertm", gertm
                   Case ("gerto")
                    gerto = par_val(ipar) !Print *,"gerto", gerto
                   Case ("gerdeb")
                    gerdeb = par_val(ipar) !Print *,"gerdeb", gerdeb  !MajValParam
                   Case ("gerkapp")
                    gerkapp = par_val(ipar) !Print *,"gerkapp", gerkapp
                   Case ("gerfin")
                    gerfin = par_val(ipar) !Print *,"gerfin", gerfin
				 Case ("htvdcan")
                    htvdcan = par_val(ipar)
				 Case ("htvdgro1")
                    htvdgro1 = par_val(ipar)
				 Case ("htvdgro2")
                    htvdgro2 = par_val(ipar)
				 Case ("htvdtb")
                    htvdtb = par_val(ipar)
				 Case ("htvdtm")
                    htvdtm = par_val(ipar)
				 Case ("htvdto")
                    htvdto = par_val(ipar)
                   Case ("humtigk1")
                    humtigk1 = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
                   Case ("humtigk2")
                    humtigk2 = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
				 Case ("humaerdeb")
                    humaerdeb = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
				 Case ("humaerfin")
                    humaerfin = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
				 Case ("humaertb")
                    humaertb = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
				 Case ("humaerdec")
                    humaerdec = par_val(ipar) !;Print *,"htvdrdeb", htvdrdeb
				 CASE ('K2ver')
	                K2ver=par_val(ipar)                   
				 CASE ('nbbgi')
	                nbbgi=par_val(ipar)
                   CASE ('kcmax')
	                kcmax=par_val(ipar)
	             CASE ('ke')
	                ke=par_val(ipar)
	             CASE ('kexrn')
	                kexrn=par_val(ipar)
	             CASE ('klaidj')
	                klaidj=par_val(ipar)
	             CASE ('ksenlai')
	                ksenlai=par_val(ipar)
                   Case ("laitb")
                      laitb = par_val(ipar) !;Print *,"laideb", laideb
                   Case ("laiwksen")
                      laiwksen = par_val(ipar) !;Print *,"laideb", laideb
                   Case ("laicroi")
                      laicroi = par_val(ipar) !;Print *,"laideb", laideb
	             CASE ('p0')
	                p0=par_val(ipar)
	             CASE ('p01')
	                p01=par_val(ipar)
	             CASE ('p02')
	                p02=par_val(ipar)
	             CASE ('p03')
	                p03=par_val(ipar)
	             CASE ('p04')
	                p04=par_val(ipar)
	             CASE ('p05')
	                p05=par_val(ipar)
	             CASE ('p06')
	                p06=par_val(ipar)
	             CASE ('p07')
	                p07=par_val(ipar)
	             CASE ('p08')
	                p08=par_val(ipar)
	             CASE ('p09')
	                p09=par_val(ipar)
	             CASE ('p10')
	                p10=par_val(ipar)
				 Case ("pblfin")
                      pblfin = par_val(ipar) !;Print *,"ph", ph
				 Case ("pbldif")
                      pbldif = par_val(ipar) !;Print *,"ph", ph
                   Case ("pblkdec")
                      pblkdec = par_val(ipar) !;Print *,"ph", ph
                   Case ("pdstmx")
                      pdstmx = par_val(ipar) !;Print *,"ph", ph
                   Case ("pflor")
                      pflor = par_val(ipar) !;Print *,"ph", ph
                   Case ("slafin")
                      slafin = par_val(ipar) !;Print *,"ph", ph
                   Case ("sladif")
                      sladif = par_val(ipar) !;Print *,"ph", ph
                   Case ("slakdec")
                      slakdec = par_val(ipar) !;Print *,"ph", ph
                   Case ("pbaset")
                      pbaset = par_val(ipar) !;Print *,"ph", ph
	             CASE ('plai')
	                plai=par_val(ipar)
                   Case ("psilb")
                     psilb = par_val(ipar) !;Print *,"kracfin", kracfin
				 Case ("pstrutb")
                     pstrutb = par_val(ipar)  !;Print *,"prue ", prue
                   Case ("pstrutcroi")
                     pstrutcroi = par_val(ipar)  !;Print *,"prue ", prue
                   Case ("pstrufin")
                     pstrufin = par_val(ipar)  !;Print *,"prue ", prue
                   Case ("pstrudec")
                     pstrudec = par_val(ipar)  !;Print *,"prue ", prue
				 Case ("ptig")
                     ptig = par_val(ipar) !;Print *,"kracdeb", kracdeb
				 Case ("ptigdeb")
                     ptigdeb = par_val(ipar) !;Print *,"kracdeb", kracdeb
				 Case ("ptigdec")
                     ptigdec = par_val(ipar) !;Print *,"kracdec", kracdec
				 Case ("ptigfin")
                     ptigfin = par_val(ipar) !;Print *,"kracfin", kracfin
	             CASE ('prue')
	                prue=par_val(ipar)
	             CASE ('profrac')
	                profrac=par_val(ipar)
	             CASE ('profracvar')
	                profracvar=par_val(ipar)
	             CASE ('ruemax')
	                ruemax=par_val(ipar)
	             CASE ('ruetk')
	                ruetk=par_val(ipar)
	             CASE ('ruetopt')
	                ruetopt=par_val(ipar)
	             CASE ('rwumx')
	                rwumx=par_val(ipar)
	             CASE ('stock0')
	                stock0=par_val(ipar)
				 Case ("sthydcroi")
                      sthydcroi = par_val(ipar) !;Print *,"ruemax", ruemax
				 Case ("sthydbio")
                      sthydbio = par_val(ipar) !;Print *,"ruemax", ruemax
			  Case ("taldebtt")
                    taldebtt = par_val(ipar) !;Print *,"talfin", talfin
                Case ("talfinval")
                    talfinval = par_val(ipar) !;Print *,"talrdeb", talrdeb
                Case ("talpeaktt")
                    talpeaktt=par_val(ipar) !;Print*,"talrdebsen",talrdebsen
                Case ("talpeakval")
                    talpeakval = par_val(ipar) !;Print *,"talrk", talrk
                Case ("taltb")
                    taltb = par_val(ipar) !;Print *,"talrksen", talrksen
	             CASE ('tbase')
	                tbase=par_val(ipar)
	             CASE ('tbaselai')	                
	                tbaselai=par_val(ipar)
	             CASE ('vitenr')
	                vitenr=par_val(ipar)
	             CASE ('vsolini')
	                vsolini=par_val(ipar)
	             CASE ('Wver')
	                Wver=par_val(ipar)
			  !******** CLC
                Case ("tuaera")
                    tuaera = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("tuaerb")
                    tuaerb = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("tuaerk")
                    tuaerk = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("ndftua")
                    ndftua = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("ndftub")
                    ndftub = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("ndftuk")
                    ndftuk = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("celtua")
                    celtua = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("celtub")
                    celtub = par_val(ipar)  !;Print *,"tbase", tbase                
			  Case ("celtuk")
                    celtuk = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("hceltua")
                    hceltua = par_val(ipar)  !;Print *,"tbase", tbase
                Case ("hceltub")
                    hceltub = par_val(ipar)  !;Print *,"tbase", tbase                
			  Case ("hceltuk")
                    hceltuk = par_val(ipar)  !;Print *,"tbase", tbase
			  Case ("ligntua")
                    ligntua = par_val(ipar)  !;Print *,"tbase", tbase
			  Case ("ligntub")
                    ligntub = par_val(ipar)  !;Print *,"tbase", tbase
			  Case ("ligntuk")
                    ligntuk = par_val(ipar)  !;Print *,"tbase", tbase

                END SELECT                     
              END DO
	    END SUBROUTINE MAJ_VAL_PARAM

      SUBROUTINE VER_FIN_PARAM
	! Determines if all parameters combinations have been used => finparam = "OUI"
        INCLUDE 'cercane1.inc' !JFM000301
            IF (typeparam > 0) THEN
              DO ipar = 1 , nbpar	          
                IF (par_indice(ipar) < par_nbre(ipar)) THEN
                  EXIT
	              END IF 	          
                IF (ipar == nbpar) Then 
                      finparam = "OUI"
	              END IF	          
              END DO
            END IF
	END SUBROUTINE VER_FIN_PARAM
