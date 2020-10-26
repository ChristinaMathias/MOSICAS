C##### ATTENTION : FCN = OBJFUN + autres modif(s) !!

C##### Voir si on devrait pas mettre "PRTVEC" en externe, on l'utilise 
C##### en ligne 180 et + !!  A été mis en externe.

C##### Voir à quel moment on initialise xopt() avec x()

C##### Modifier LB et UB : se sont des éléments qui pour chaque I, ie
C##### pour chaque paramètre, vont normalement dépendre de l'écart en %
C##### etdu nombre d'écarts de part et d'autre de la valeur initiale 
C##### donnée. 1 intervalle pour chaque paramètre donné.

C##### nt = le nombre d'itérations : à donner par nous, à poser(?).
C##### ns = le nombre de cycles : après ns*n, on a ajustement de vm.
C##### n = le nombre de paramètres.   
C##### après nt*ns*n évaluations de fonction, modification de la T.

C##### ATTENTION : il faut définir xopt
C     PROGRAM SIMANN
C	IMPLICIT NONE	
      EXTERNAL OBJFUN
      include 'denopo.inc'               
	include 'cercane1.inc'
	EXTERNAL PRTVEC  ! ajouté en externe car utilisé.
      common /recuit/nob 
C      common /simplex/ nob
	integer nob,npar,i,maxit 
C	integer io1
	integer N
	character  s1*50 
	real, allocatable :: yob(:) 
	real, allocatable :: parametre(:)
	character, allocatable, dimension(:) :: varoptob*15  !JFM130116

	DOUBLE PRECISION, allocatable, dimension(:) :: LB
	DOUBLE PRECISION, allocatable, dimension(:) :: UB
	DOUBLE PRECISION, allocatable, dimension(:) :: PRODUIT
	DOUBLE PRECISION, allocatable, dimension(:) :: X
	DOUBLE PRECISION, allocatable, dimension(:) :: XOPT
	DOUBLE PRECISION, allocatable, dimension(:) :: C
	DOUBLE PRECISION, allocatable, dimension(:) :: VM
	DOUBLE PRECISION, allocatable, dimension(:) :: XP
	INTEGER, allocatable, dimension(:) :: NACP
C##### déclaration pour le Recuit Simulé
!	PARAMETER (N = 5, NEPS = 4)
      INTEGER, PARAMETER :: NEPS=4
      DOUBLE PRECISION     T, EPS, RT, FOPT, FSTAR(NEPS)
      INTEGER  NS, NT, NFCNEV, IER, ISEED1, ISEED2,
     1         MAXEVL, IPRINT, NACC, NOBDS


    ! ou : 
!      DOUBLE PRECISION  LB(N), UB(n),  X(n), XOPT(n), C(n) , VM(n),
!     1                  FSTAR(neps), XP(n) !, T, EPS, RT, FOPT
!      INTEGER  NACP(n) !, NS, NT, NFCNEV, IER, ISEED1, ISEED2,
    ! 1         MAXEVL, IPRINT, NACC, NOBDS

	LOGICAL  MAX
	LOGICAL BEXIST
C##### fin déclaration-recuit


C      LOGICAL FIRST ! voir si "FIRST" aussi ds recuit
c##### déclaration Simplex
c      DOUBLE PRECISION OBJECT, SIMP, STOPCR, P(20),STEP(20),VAR(20)
c	INTEGER IER, IPRINT, IQUAD, MAXF, NLOOP, NOP, J
c      DOUBLE PRECISION ye(500) 
c##### fin déclaration Simplex

	

C***************** Détermination du type de simulation ********************

	!********** 0 normal, 1 sensibilité ou 2 calage (Haus)
	  !write(*,*)'Press any key before simulation'
	  !read(*,*)
	  !BEXIST=.TRUE.
	  !Inquire (FILE='sim.txt',EXIST=BEXIST)
	  !IF (.NOT.BEXIST) then
	  !	write(*,*) 'sim.txt not found!'; read(*,*);stop
	  !END IF


		versionmos='VERSION mosicas-2017.05.29 (J.F. MARTINE)'


 	   OPEN(13,File='sim.txt',err=99)  !16/7/03
		read(13,*,err=99);read(13,*);read(13,*);READ(13,*) !16/7/03
		read(13,*);read(13,*) !jfm150513
	   read(13,*);read(13,*);read(13,*);READ(13,*)i; close(13)	
	goto 100
99	write(*,*) "Lecture de sim.txt impossible  (main_recuit)";	stop	
	   
C        OPEN(13,File='sim.dat'); READ(13,*)s1; READ(13,*)i; close(13)
			!write(*,*) i	
100	IF (i .eq. 2) THEN !calage 
			write(*,*)"calage"
			nbrecuit=0
			elapsed_time= timeF()
			call motsimul !remplit dnobs
	
	OPEN(21,file='recuit.txt')
 !     write (21,*) '########  APPEL DU RECUIT-SIMULE  ########
		else if ( (i .eq. 0) .or. (i .eq. 1)) then !normal ou sensibilité
			nbrecuit=0
		        elapsed_time= timeF()

			call motsimul
		        elapsed_time= timeF()
			write(*,*) 'time s/lect  = ', elaps1
			write(*,*) 'time s/Total = ', elapsed_time, elapsed_time/60
			elaps1= elapsed_time - elaps1
			write(*,*) 'time s/trait = ', elaps1, elaps1/nbtrt
			write(*,*) "END OF SIMULATION"
			write (*,*) versionmos
			!write (*,*) 'press any key to continue' 
!			read(*,*)
			stop 'END OF RUN'
		else
			write(*,*) "type de simulation inconnu! Choisir 0, 1 ou 2"
			stop
	END IF
	CLOSE(13)









C ******************** Si CALAGE / OPTIMISATION *****************************
	!*********   détermine npar et remplit parametre(i)
!!	OPEN(9,file='PARCAL.txt',err=101) 
!!		rewind(9) 
		!old read (9,FMT='(I2,1X,A15)') npar,varopt 
!!		read (9,*,err=101) npar,varopt !new
!!	goto 102
!!101	write(*,*) "Lecture de parcal.txt impossible  (main_recuit)";	stop	
!!102		allocate(parametre(npar)) 
!!		do i=1,npar 
!!		read(9,*)paropt(i),parametre(i),varipar(i),nvaripar(i) !new
!!   		end do 
!!	close(9) 
	OPEN(9,file='PARCAL.txt')  !JFM13017
		rewind(9) 
		read(9,*) nvaropt          		
	    do i=1,nvaropt  
	      read(9,*) varopt(i), varoptpds(i)    
		end do
		read(9,*) npar
		allocate(parametre(npar)) 		        
		do i=1,npar 
		read(9,*)paropt(i),parametre(i),varipar(i),nvaripar(i) 
   		end do 
	close(9)   !JFM13017







	!*********   détermine nob et yob()
!!	nob=0
!!	OPEN(10,file='dnobs.txt')!new !new
!!4			READ (10,FMT='(F13.6)',end=5) valeur !new
!!			nob=nob+1 !new !determine le nb d'observations
!!		 goto 4 !new
!!5			continue	 !new		
!!			allocate(yob(nob)); 
!!              REWIND(10) !new
!!			DO i=1,nob
!!				read(10,FMT='(F13.6)') yob(i)
!!		    END DO !new
	!détermine la valeur des observations et remplit yob()
!!	Close(10) !new
	nob=nobs; allocate(yob(nob)); yob(:)=dnobst(:)   !JFM130117



	!do i=1,nob 	write(*,*)i,dnobst(i),simtempt(i), nvobst(i) end do 

C#####
    !  ds le recuit N représente le nombre de paramètres ie npar.
	N=npar
C  Set underflows to zero on IBM mainframes.
C     CALL XUFLOW(0)
	!ouverture d'un fichier resultat
	OPEN(21,file='recuit.txt')
C  Set input parameters.
      MAX = .FALSE. ! signifie minimisation ; si max, alors .true.
	EPS = 1.0D-6  ! définit un critère d'arrêt sur le nombre de dernières 
				  !	apparitions d'un résultat 
      RT = .85  ! facteur de recuit (la température décroit de façon géométriq)
C##### RT=paramètre que l'on peut ajuster pour améliorer la vitesse et la robustesse.
      ISEED1 = 1 ! 1ère graine pour le générateur de nb aléatoires
      ISEED2 = 2 ! 2ème graine pour le générateur de nb aléatoires
      NS = 20  ! nombre de cycles 20
      NT = 5   ! nombre d'itérations avant la réduction de température 5  // 10 avant JFM130118
C##### NT=paramètre que l'on peut ajuster pour améliorer la vitesse et la robustesse.
      MAXEVL = 1000 !1000 ! nombre max d'évaluations de fonctions
	MAXEVL = (NS*NT*npar)*2
      IPRINT = 1 ! permet d'avoir des détails sur l'affichage, sinon 0,1 ou 2
	allocate(LB(npar)); 	allocate(UB(npar)) 	
	allocate(C(npar));  	allocate(PRODUIT(npar))
	!DO 10, I = 1, npar
       !  LB(I) = 0.1 !-1.0D25  ! pour un I donné, [LB(I), UB(I)] détermine 
	 !  UB(I) = 0.1 ! 1.0D25  ! l'intervalle auquel appartient le paramètre I
       ! C(I) = 2.0       ! le définir !!!!
!10    CONTINUE
C##### Initialisation des domaines d'appartenance des paramètres.
C##### L'intervalle doit être modifier pour chaque paramètre.
C##### On détermine [lb, ub] à partir de sim.dat, de la valeur du
C##### paramètre, du % d'écart et du nombre d'écarts autour de la valeur ini. 
C##### s_err='Ouverture fichier impossible:sim.dat'
!!	OPEN(13,File='sim.txt') !,err=99) 16/7/03
	!READ (13,*) avant 16/7/03 
!!	read(13,*);read(13,*);read(13,*);READ(13,*) ! le 16/7/03
!!	read(13,*);read(13,*);read(13,*) ! le 16/7/03
!!	READ (13,*) typeparam,typecalage !type de simul.
	i=npar; PRODUIT(1:i) = (varipar(1:i) * nvaripar(1:i))/100	!JFM130117
	LB(1:i) = parametre(1:i)-(PRODUIT(1:i) * parametre(1:i))	!JFM130117
	UB(1:i) = parametre(1:i)+(PRODUIT(1:i) * parametre(1:i))	!JFM130117
	C(1:i) = 2.0	!JFM130117
	do i=1,npar !JFM130208
	IF (LB(i)<=0.0) LB(i)=0.000001 !JFM130208
	end do	
	WRITE(*,*) "lb                        ub"
	do i=1,npar !JFM130117
	WRITE(*,*) lb(i), ub(i) 
	end do
!!		DO 15, i=1,npar 
!!			PRODUIT(i) = varipar(i) * nvaripar(i)
!!			PRODUIT(i) = PRODUIT(i) / 100
!!			LB(i) = parametre(i)-(PRODUIT(i) * parametre(i))
!!			UB(i) = parametre(i)+(PRODUIT(i) * parametre(i))
!!			C(i) = 2.0  
!!			WRITE(*,*) "lb                        ub"
!!			WRITE(*,*) lb(i), ub(i)
!!   15		continue

!!      CLOSE(13)
!	write(*,*) npar

  
C####################################
C  Note start at local, but not global, optima of the Judge function.
	allocate(X(npar)); X(:)=parametre(:)  !JFM130117
!!		DO 16, J=1,N
!!			X(J)=parametre(J) ! redéfinit les paramètres pour + d'accords
!!16		continue	
!	end do		! @@@ici problème de tableau memo.

    !  X(1) =  2.354471
    !  X(2) = -0.319186
C  Set input values of the input/output parameters.
      T = 5.0 ! température initiale du recuit
	ALLOCATE(VM(N)); VM(:)=0.1
 !!     DO 20, I = 1, N
 !!        VM(I) = 0.1 ! 1.0 ! le définir !???????????????????!!!
 !!20   CONTINUE



C##### écriture ds le fichier "recuit.txt" des initialisations. 
 !     write(21,1000) N, MAX, T, RT, EPS, NS, NT, NEPS, MAXEVL, IPRINT,
 !    1              ISEED1, ISEED2
      do i=1,npar
	write(21,*) paropt(i) !jfm030929
	end do

!	CALL PRTVEC(X,N,'STARTING VALUES (X)') !old before JFM090930
      write(21,*) "STARTING VALUES (X)"    !new JFM090930
	 DO i = 1, N                         !new JFM090930
		write(21,*) X(i)                 !new JFM090930
      END DO                               !new JFM090930

   !   CALL PRTVEC(VM,N,'INITIAL STEP LENGTH (VM)')
!      CALL PRTVEC(LB,N,'LOWER BOUND (LB)') !old before JFM090930
      write(21,*) "LOWER BOUND (LB)"        !new JFM090930
	 DO i = 1, N                          !new JFM090930
		write(21,*) LB(i)                 !new JFM090930
      END DO                                !new JFM090930

!      CALL PRTVEC(UB,N,'UPPER BOUND (UB)') !old before JFM090930
      write(21,*) "UPPER BOUND (UB)"        !new JFM090930
	 DO i = 1, N                          !new JFM090930
		write(21,*) UB(i)                 !new JFM090930
      END DO                                !new JFM090930



   !   CALL PRTVEC(C,N,'C VECTOR')
 !     write(21,'(/,''  ****   END OF DRIVER ROUTINE OUTPUT   ****''
 !    1          /,''  ****   BEFORE CALL TO SA.             ****'')')
           
!	write(*,*) "avant passage de SA : ok"

      allocate(xopt(n))


C ####### temporaire test  JFM130117
!!	close (21)
!!	write(*,*) ""
!!	write(*,*) "FIN TEST"
!!      STOP
!!      END






C##### appel de la subroutine  de recuit :   !! JFM130117 
      CALL SA(N,X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,C,IPRINT,ISEED1,
     1        ISEED2,T,VM, XOPT,
     2		FOPT,NACC,NFCNEV,NOBDS,IER,FSTAR)


!	write(*,*) " après passage de SA : ok"

C####################################################
 !     write(21,'(/,''  ****   RESULTS AFTER SA   ****   '')'	allocate(Xopt(n)))   
 !     CALL PRTVEC(XOPT,N,'SOLUTION') !old before JFM090930
      write(21,*) "SOLUTION"    !new JFM090930    !! JFM130117
	 DO i = 1, N              !new JFM090930     !! JFM130117
		write(21,*) xopt(i)   !new JFM090930     !! JFM130117
      END DO                    !new JFM090930   !! JFM130117

  !    CALL PRTVEC(VM,N,'FINAL STEP LENGTH')
      write(21,*) 'OPTIMAL FUNCTION VALUE: '   !jfm090930         !! JFM130117
	write(21,*) FOPT !, NFCNEV, NACC, NOBDS, T, IER !jfm090930    !! JFM130117

!1000  FORMAT(/,' NUMBER OF PARAMETERS: ',I3,'   MAXIMAZATION: ',L5,
!     1       /,' INITIAL TEMP: ', G8.2, '   RT: ',G8.2, '   EPS: ',G8.2,
!     2       /,' NS: ',I3, '   NT: ',I2, '   NEPS: ',I2,
!     3       /,' MAXEVL: ',I10, '   IPRINT: ',I1, '   ISEED1: ',I4,
!     4       '   ISEED2: ',I4)
!1001  FORMAT(/,' OPTIMAL FUNCTION VALUE: ',G25.18 ) oldjfm
!1001  FORMAT(/,G25.18 ) !jfm090930
!     1       /,' NUMBER OF FUNCTION EVALUATIONS:     ',I10,
 !    2       /,' NUMBER OF ACCEPTED EVALUATIONS:     ',I10,
  !   3       /,' NUMBER OF OUT OF BOUND EVALUATIONS: ',I10,
   !  4       /,' FINAL TEMP: ', G20.13,'  IER: ', I3)

   !   STOP



C#############################################
C initialisation Simplex
C  NOP is number of variables. NOP <= 20 (see subroutine NELMEAD)
C	  NOP = npar
C  Set up starting values & step sizes.
C      DO J = 1, NOP
C	    P(J) = parametre(j)
C        STEP(J) = 0.25*parametre(j)
C      END DO
C  Set max. no. of function evaluations = 500, print every 20.
C      MAXF = 500
C      IPRINT = 20
C  Set value for stopping criterion.   Stopping occurs when the
C  standard deviation of the values of the objective function at
C  the points of the current simplex < stopcr.
C      STOPCR = 1.D-5
C  Fit a quadratic surface to be sure a minimum has been found.
C      IQUAD = 1
C  As function value is being evaluated in DOUBLE PRECISION, it
C  should be accurate to about 15 decimals.   If we set simp = 1.d-6,
C  we should get about 9 dec. digits accuracy in fitting the surface.
C      SIMP = 1.D-6
C appel du Simplex
C         CALL NELMEAD (P, STEP, NOP, OBJECT, MAXF, IPRINT, STOPCR, 
C     1      NLOOP, IQUAD, SIMP, VAR, OBJFUN, IER)
C##############################################################
   !   CLOSE (3)
C     close(18)
	CLOSE (21)														!! JFM130117
	write(*,*) "F à la fin du main="								!! JFM130117
!	write(*,*) f													!! JFM130117
	elapsed_time= timeF()
	write(*,*) 'time (s/mn) = ', elapsed_time, elapsed_time/60		!! JFM130117
	write(*,*) 'Nb recuits = ', nbrecuit							!! JFM130117
      STOP															!! JFM130117
      END															!! JFM130117
C	END 








C#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C##### Fonction objectif pour le Recuit-Simulé
      !subroutine OBJFUN(NOP,P,FUNC)
	subroutine OBJFUN(N,P,F)
    !  IMPLICIT NONE
	include 'denopo.inc'
      common /recuit/nob 
!      common /simplex/ nob
!	 cf. N=nop
      DOUBLE PRECISION P(N), F, ye(nob), yob(nob)
	!DOUBLE PRECISION FAC
      INTEGER N,I,J,NOB
!	write(*,*) "nob vaut="
!	write(*,*) nob
	!*** Ecrire les valeurs parametre() et noms dans un fichier parcal.txt()
!!	OPEN(9,file='PARCAL.txt') !JFM090930 PARCAL.txt PARCALtemp.txt Intérêt puisque ce n'est plus lu autre part ????
!!        rewind(9)
!!		write(9,*)N,achar(9),adjustl(varopt)
!!		do i=1,N
!!	 	write(9,*) trim(adjustl(paropt(i))),P(i)
!!		end do
!!	close(9)

	OPEN(9,file='PARCAL.txt') !JFM130118
		rewind(9)
		write(9,FMT='(i1)')nvaropt
		do i=1,nvaropt
			write(9,*) adjustl(varopt(i)), varoptpds(i)
		end do
		write(9,FMT='(i1)')N
		do i=1,N
	 		write(9,*) trim(adjustl(paropt(i))),P(i)
		end do
	close(9)


	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	!NOB EST INCONNU!!!!!!!!!! : essai après
 
!!	OPEN(10,file='dnobs.txt')!new !new
!!			DO i=1,nob
!!				read(10,FMT='(F13.6)') yob(i)
!!			!	write(*,*) "test sur yob"
!!			!	write(*,*) yob(i)
!!		    END DO !new
	!détermine la valeur des observations et remplit yob()
!!	Close(10) !new
	yob(:)=dnobst(:)


	!NOB EST INCONNU!!!!!!!!!! : essai avant
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		
C*** Appeler motsimul et inclure dans cercane 
	  
	  call motsimul

	!*** Lecture de simtemp() à mettre dans ye() et à écrire dans 3,dsim
!!	OPEN(12,file='simtemp.txt')!new !new
!!			DO i=1,nob
!!			   read(12,FMT='(F13.6)')ye(i); 
			  ! write (3,*) ye(i) 
!!			END DO !new
			!write(*,*) ye(1) 
!!	Close(12) !new
	ye(:)=simtempt(:)



c calcule la valeur de la fonction objectif
c ye = simulées; yob = observées
      F=0.
!	write(*,*) ye(1), yob(1)
      do 100 j=1,nob
         F = F + (ye(j)-yob(j))**2
 100     continue
         F = F/nob
C        write (18,*) 'Objective func.',func,ye(5),yob(5)

	return
      end

C	end program
	

C********************************************************


C##### contient la subroutine de recuit et d'autres pour ecriture
C##### ATTENTION : 	EXTERNAL OBJFUN ou pas ????????

C##### Faut-il redefinir le neps ici ? ie sa valeur ? pas fait mais ds
C##### fait ds le main_recuit.for.

C##### Paramètres à modifier pour un algo + rapide et + robustes : NT et RT.


C ABSTRACT:
C   Simulated annealing is a global optimization method that distinguishes
C   between different local optima. Starting from an initial point, the
C   algorithm takes a step and the function is evaluated. When minimizing a
C   function, any downhill step is accepted and the process repeats from this
C   new point. An uphill step may be accepted. Thus, it can escape from local
C   optima. This uphill decision is made by the Metropolis criteria. As the
C   optimization process proceeds, the length of the steps decline and the
C   algorithm closes in on the global optimum. Since the algorithm makes very
C   few assumptions regarding the function to be optimized, it is quite
C   robust with respect to non-quadratic surfaces. The degree of robustness
C   can be adjusted by the user. In fact, simulated annealing can be used as
C   a local optimizer for difficult functions.
C
C   This implementation of simulated annealing was used in "Global Optimization
C   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
C   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
C   65-100. Briefly, we found it competitive, if not superior, to multiple
C   restarts of conventional optimization routines for difficult optimization
C   problems.
C
C   For more information on this routine, contact its author:
C   Bill Goffe, bgoffe@whale.st.usm.edu
C
C  This file is an example of the Corana et al. simulated annealing
C  algorithm for multimodal and robust optimization as implemented
C  and modified by Goffe, Ferrier and Rogers. Counting the above line
C  ABSTRACT as 1, the routine itself (SA), with its supplementary
C  routines, is on lines 232-990. A multimodal example from Judge et al.
C  (FCN) is on lines 150-231. The rest of this file (lines 1-149) is a
C  driver routine with values appropriate for the Judge example. Thus, this
C  example is ready to run.
C
C  To understand the algorithm, the documentation for SA on lines 236-
C  484 should be read along with the parts of the paper that describe
C  simulated annealing. Then the following lines will then aid the user
C  in becomming proficient with this implementation of simulated
C  annealing.
C
C  Learning to use SA:
C      Use the sample function from Judge with the following suggestions
C  to get a feel for how SA works. When you've done this, you should be
C  ready to use it on most any function with a fair amount of expertise.
C    1. Run the program as is to make sure it runs okay. Take a look at
C       the intermediate output and see how it optimizes as temperature
C       (T) falls. Notice how the optimal point is reached and how
C       falling T reduces VM.
C    2. Look through the documentation to SA so the following makes a
C       bit of sense. In line with the paper, it shouldn't be that hard
C       to figure out. The core of the algorithm is described on pp. 68-70
C       and on pp. 94-95. Also see Corana et al. pp. 264-9.
C    3. To see how it selects points and makes decisions about uphill
C       and downhill moves, set IPRINT = 3 (very detailed intermediate
C       output) and MAXEVL = 100 (only 100 function evaluations to limit
C       output).
C    4. To see the importance of different temperatures, try starting
C       with a very low one (say T = 10E-5). You'll see (i) it never
C       escapes from the local optima (in annealing terminology, it
C       quenches) & (ii) the step length (VM) will be quite small. This
C       is a key part of the algorithm: as temperature (T) falls, step
C       length does too. In a minor point here, note how VM is quickly
C       reset from its initial value. Thus, the input VM is not very
C       important. This is all the more reason to examine VM once the
C       algorithm is underway.
C    5. To see the effect of different parameters and their effect on
C       the speed of the algorithm, try RT = .95 & RT = .1. Notice the
C       vastly different speed for optimization. Also try NT = 20. Note
C       that this sample function is quite easy to optimize, so it will
C       tolerate big changes in these parameters.
C
C	  !!!!!!!!!!!!!  Important et très utile : !!!!!!!!!!!!!!!!
C	  RT and NT are the parameters one should adjust to modify the 
C       runtime of the algorithm and its robustness.
C
C    6. Try constraining the algorithm with either LB or UB.


      SUBROUTINE SA(N,X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,C,IPRINT,
     1              ISEED1,ISEED2,T,VM,XOPT, FOPT,NACC,NFCNEV,NOBDS,IER,
     2              FSTAR)


C  Version: 3.2
C  Date: 1/22/94.
C  Differences compared to Version 2.0:
C     1. If a trial is out of bounds, a point is randomly selected
C        from LB(i) to UB(i). Unlike in version 2.0, this trial is
C        evaluated and is counted in acceptances and rejections.
C        All corresponding documentation was changed as well.
C  Differences compared to Version 3.0:
C     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
C        The idea is that if T is high relative to LB & UB, most
C        points will be accepted, causing VM to rise. But, in this
C        situation, VM has little meaning; particularly if VM is
C        larger than the acceptable region. Setting VM to this size
C        still allows all parts of the allowable region to be selected.
C  Differences compared to Version 3.1:
C     1. Test made to see if the initial temperature is positive.
C     2. WRITE statements prettied up.
C     3. References to paper updated.
C
C  Synopsis:
C  This routine implements the continuous simulated annealing global
C  optimization algorithm described in Corana et al.'s article
C  "Minimizing Multimodal Functions of Continuous Variables with the
C  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
C  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
C  Software.
C
C  A very quick (perhaps too quick) overview of SA:
C  But de l'algo : SA tries to find the global optimum of an N 
C  dimensional function.
C  It moves both up and downhill and as the optimization process
C  proceeds, it focuses on the most promising area.
C     To start, it randomly chooses a trial point within the step length
C  VM (a vector of length N) of the user selected starting point. The
C  function is evaluated at this trial point and its value is compared
C  to its value at the initial point.
C     In a maximization problem, all uphill moves are accepted and the
C  algorithm continues from that trial point. Downhill moves may be
C  accepted; the decision is made by the Metropolis criteria. It uses T
C  (temperature) and the size of the downhill move in a probabilistic
C  manner. The smaller T and the size of the downhill move are, the more
C  likely that move will be accepted. If the trial is accepted, the
C  algorithm moves on from that point. If it is rejected, another point
C  is chosen instead for a trial evaluation.
C     Each element of VM periodically adjusted so that half of all
C  function evaluations in that direction are accepted.
C     A fall in T is imposed upon the system with the RT variable by
C  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
C  downhill moves are less likely to be accepted and the percentage of
C  rejections rise. Given the scheme for the selection for VM, VM falls.
C  Thus, as T declines, VM falls and SA focuses upon the most promising
C  area for optimization.
C
C  The importance of the parametre T:
C     The parametre T is crucial in using SA successfully. It influences
C  VM, the step length over which the algorithm searches for optima. For
C  a small intial T, the step length may be too small; thus not enough
C  of the function might be evaluated to find the global optima. The user
C  should carefully examine VM in the intermediate output (set IPRINT =
C  1) to make sure that VM is appropriate. The relationship between the
C  initial temperature and the resulting step length is function
C  dependent.
C     To determine the starting temperature that is consistent with
C  optimizing a function, it is worthwhile to run a trial run first. Set
C  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
C  rises as well. Then select the T that produces a large enough VM.
C
C  For modifications to the algorithm and many details on its use,
C  (particularly for econometric applications) see Goffe, Ferrier
C  and Rogers, "Global Optimization of Statistical Functions with
C  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2, 
C  Jan./Feb. 1994, pp. 65-100.
C  For more information, contact 
C              Bill Goffe
C              Department of Economics and International Business
C              University of Southern Mississippi 
C              Hattiesburg, MS  39506-5072 
C              (601) 266-4484 (office)
C              (601) 266-4920 (fax)
C              bgoffe@whale.st.usm.edu (Internet)
C
C  As far as possible, the parameters here have the same name as in
C  the description of the algorithm on pp. 266-8 of Corana et al.
C
C  In this description, SP is single precision, DP is double precision,
C  INT is integer, L is logical and (N) denotes an array of length n.
C  Thus, DP(N) denotes a double precision array of length n.
C
C  Input Parameters:
C    Note: The suggested values generally come from Corana et al. To
C          drastically reduce runtime, see Goffe et al., pp. 90-1 for
C          suggestions on choosing the appropriate RT and NT.
C    N - Number of variables in the function to be optimized. (INT)
C    X - The starting values for the variables of the function to be
C        optimized. (DP(N))
C    MAX - Denotes whether the function should be maximized or
C          minimized. A true value denotes maximization while a false
C          value denotes minimization. Intermediate output (see IPRINT)
C          takes this into account. (L)
C    RT - The temperature reduction factor. The value suggested by
C         Corana et al. is .85. See Goffe et al. for more advice. (DP)
C    EPS - Error tolerance for termination. If the final function
C          values from the last neps temperatures differ from the
C          corresponding value at the current temperature by less than
C          EPS and the final function value at the current temperature
C          differs from the current optimal function value by less than
C          EPS, execution terminates and IER = 0 is returned. (EP)
C    NS - Number of cycles. After NS*N function evaluations, each
C         element of VM is adjusted so that approximately half of
C         all function evaluations are accepted. The suggested value
C         is 20. (INT)
C    NT - Number of iterations before temperature reduction. After
C         NT*NS*N function evaluations, temperature (T) is changed
C         by the factor RT. Value suggested by Corana et al. is
C         MAX(100, 5*N). See Goffe et al. for further advice. (INT)
C    NEPS - Number of final function values used to decide upon termi-
C           nation. See EPS. Suggested value is 4. (INT)
C    MAXEVL - The maximum number of function evaluations. If it is
C             exceeded, IER = 1. (INT)
C    LB - The lower bound for the allowable solution variables. (DP(N))
C    UB - The upper bound for the allowable solution variables. (DP(N))
C         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
C         I = 1, N, a point is from inside is randomly selected. This
C         This focuses the algorithm on the region inside UB and LB.
C         Unless the user wishes to concentrate the search to a par-
C         ticular region, UB and LB should be set to very large positive
C         and negative values, respectively. Note that the starting
C         vector X should be inside this region. Also note that LB and
C         UB are fixed in position, while VM is centered on the last
C         accepted trial set of variables that optimizes the function.
C    C - Vector that controls the step length adjustment. The suggested
C        value for all elements is 2.0. (DP(N))
C    IPRINT - controls printing inside SA. (INT)
C             Values: 0 - Nothing printed.
C                     1 - Function value for the starting value and
C                         summary results before each temperature
C                         reduction. This includes the optimal
C                         function value found so far, the total
C                         number of moves (broken up into uphill,
C                         downhill, accepted and rejected), the
C                         number of out of bounds trials, the
C                         number of new optima found at this
C                         temperature, the current optimal X and
C                         the step length VM. Note that there are
C                         N*NS*NT function evalutations before each
C                         temperature reduction. Finally, notice is
C                         is also given upon achieveing the termination
C                         criteria.
C                     2 - Each new step length (VM), the current optimal
C                         X (XOPT) and the current trial X (X). This
C                         gives the user some idea about how far X
C                         strays from XOPT as well as how VM is adapting
C                         to the function.
C                     3 - Each function evaluation, its acceptance or
C                         rejection and new optima. For many problems,
C                         this option will likely require a small tree
C                         if hard copy is used. This option is best
C                         used to learn about the algorithm. A small
C                         value for MAXEVL is thus recommended when
C                         using IPRINT = 3.
C             Suggested value: 1
C             Note: For a given value of IPRINT, the lower valued
C                   options (other than 0) are utilized.
C    ISEED1 - The first seed for the random number generator RANMAR.
C             0 .LE. ISEED1 .LE. 31328. (INT)
C    ISEED2 - The second seed for the random number generator RANMAR.
C             0 .LE. ISEED2 .LE. 30081. Different values for ISEED1
C             and ISEED2 will lead to an entirely different sequence
C             of trial points and decisions on downhill moves (when
C             maximizing). See Goffe et al. on how this can be used
C             to test the results of SA. (INT)
C
C  Input/Output Parameters:
C    T - On input, the initial temperature. See Goffe et al. for advice.
C        On output, the final temperature. (DP)
C    VM - The step length vector. On input it should encompass the
C         region of interest given the starting value X. For point
C         X(I), the next trial point is selected is from X(I) - VM(I)
C         to  X(I) + VM(I). Since VM is adjusted so that about half
C         of all points are accepted, the input value is not very
C         important (i.e. is the value is off, SA adjusts VM to the
C         correct value). (DP(N))
C
C  Output Parameters:
C    XOPT - The variables that optimize the function. (DP(N))
C    FOPT - The optimal value of the function. (DP)
C    NACC - The number of accepted function evaluations. (INT)
C    NFCNEV - The total number of function evaluations. In a minor
C             point, note that the first evaluation is not used in the
C             core of the algorithm; it simply initializes the
C             algorithm. (INT).
C    NOBDS - The total number of trial function evaluations that
C            would have been out of bounds of LB and UB. Note that
C            a trial point is randomly selected between LB and UB.
C            (INT)
C    IER - The error return number. (INT)
C          Values: 0 - Normal return; termination criteria achieved.
C                  1 - Number of function evaluations (NFCNEV) is
C                      greater than the maximum number (MAXEVL).
C                  2 - The starting value (X) is not inside the
C                      bounds (LB and UB).
C                  3 - The initial temperature is not positive.
C                  99 - Should not be seen; only used internally.
C
C  Work arrays that must be dimensioned in the calling routine:
C       RWK1 (DP(NEPS))  (FSTAR in SA)
C       RWK2 (DP(N))     (XP    "  " )
C       IWK  (INT(N))    (NACP  "  " )
C
C  Required Functions (included):
C    EXPREP - Replaces the function EXP to avoid under- and overflows.
C             It may have to be modified for non IBM-type main-
C             frames. (DP)
C    RMARIN - Initializes the random number generator RANMAR.
C    RANMAR - The actual random number generator. Note that
C             RMARIN must run first (SA does this). It produces uniform
C             random numbers on [0,1]. These routines are from
C             Usenet's comp.lang.fortran. For a reference, see
C             "Toward a Universal Random Number Generator"
C             by George Marsaglia and Arif Zaman, Florida State
C             University Report: FSU-SCRI-87-50 (1987).
C             It was later modified by F. James and published in
C             "A Review of Pseudo-random Number Generators." For
C             further information, contact stuart@ads.com. These
C             routines are designed to be portable on any machine
C             with a 24-bit or more mantissa. I have found it produces
C             identical results on a IBM 3081 and a Cray Y-MP.
C
C  Required Subroutines (included):
C    PRTVEC - Prints vectors.
C    PRT1 ... PRT10 - Prints intermediate output.
C    FCN - Function to be optimized. The form is
C            SUBROUTINE FCN(N,X,F)
C            INTEGER N
C            DOUBLE PRECISION  X(N), F
C            ...
C            function code with F = F(X)
C            ...
C            RETURN
C            END
C          Note: This is the same form used in the multivariable
C          minimization algorithms in the IMSL edition 10 library.
C
C  Machine Specific Features:
C    1. EXPREP may have to be modified if used on non-IBM type main-
C       frames. Watch for under- and overflows in EXPREP.
C    2. Some FORMAT statements use G25.18; this may be excessive for
C       some machines.
C    3. RMARIN and RANMAR are designed to be portable; they should not
C       cause any problems.
    !  IMPLICIT NONE


C################################################################
C  Type all external variables.



    !  DOUBLE PRECISION  X(*), LB(*), UB(*), C(*), VM(*), FSTAR(*), 
    !					  !fstar: pas externe? même chose pour xopt, xp.
    ! 1                  XOPT(*), XP(*), T, EPS, RT, FOPT
    !
    !  INTEGER   N, NS, NT, NEPS, NACC, MAXEVL, IPRINT,
    ! 1         NOBDS, IER, NFCNEV, ISEED1, ISEED2, NACP(*)
  
	

C ou bien :


    !  DOUBLE PRECISION, allocatable, dimension(:) :: X
    !  DOUBLE PRECISION, allocatable, dimension(:) :: LB
    !  DOUBLE PRECISION, allocatable, dimension(:) :: UB
    !  DOUBLE PRECISION, allocatable, dimension(:) :: C
    !  DOUBLE PRECISION, allocatable, dimension(:) :: VM
    !  DOUBLE PRECISION, allocatable, dimension(:) :: FSTAR
    !  DOUBLE PRECISION, allocatable, dimension(:) :: XOPT
    !  DOUBLE PRECISION, allocatable, dimension(:) :: XP
    !  INTEGER, ALLOCATABLE, dimension(:) :: nacp
    !  INTEGER   N, NS, NT, NEPS, NACC, MAXEVL, IPRINT,
    ! 1         NOBDS, IER, NFCNEV, ISEED1, ISEED2
!	double precision :: T, EPS, RT, FOPT


C ou bien :

      DOUBLE PRECISION, intent (inout) :: X(n)
      DOUBLE PRECISION, intent (in) :: LB(n)
      DOUBLE PRECISION, intent (in) :: UB(n)
      DOUBLE PRECISION, intent (in) :: C(n)
      DOUBLE PRECISION, intent (inout) :: VM(n)
      DOUBLE PRECISION, intent (inout) :: FSTAR(NEPS)
      DOUBLE PRECISION:: XOPT(n)
      DOUBLE PRECISION :: XP(n) !
      INTEGER, allocatable :: nacp(:)
      INTEGER :: N, NS, NT, NEPS, NACC, MAXEVL, IPRINT,
     1         NOBDS, IER, NFCNEV, ISEED1, ISEED2
	DOUBLE PRECISION :: T, EPS, RT, FOPT



C#############################################################

C##### ajouté car utilisé ici, mais pas nécessaire si déjà fait
C##### ds le main?
!	PARAMETER (NEPS = 4)
!     INTEGER, PARAMETER :: neps=4
!	EXTERNAL OBJFUN
      LOGICAL  MAX

C  Type all internal variables.
      DOUBLE PRECISION  F, FP, P, PP, RATIO
      INTEGER  NUP, NDOWN, NREJ, NNEW, LNOBDS, H, I, J, M
      LOGICAL  QUIT

C  Type all functions.
      DOUBLE PRECISION  EXPREP
      REAL  RANMAR

C  Initialize the random number generator RANMAR.
      CALL RMARIN(ISEED1,ISEED2)
	OPEN(21,file='recuit.txt')
C  Set initial values.
      NACC = 0
      NOBDS = 0
      NFCNEV = 0
      IER = 99
!	write(*,*) "n="
!	write(*,*) n
C      MAX = .FALSE. ! signifie minimisation
C##### Initialisation du vecteur x optimal aux valeurs du début de recherche.
    !  allocate(XOPT(N))
	allocate(nacp(N))
	DO 10, i = 1, N
	 	XOPT(i) = X(i)      
	    NACP(i) = 0  !définir le sens.
10    CONTINUE

      DO 20, I = 1, NEPS
         FSTAR(I) = 1.0D+20
20    CONTINUE 

C  If the initial temperature is not positive, notify the user and 
C  return to the calling routine.
C##### Test sur la température de départ, si <0, alors erreur.  
    !  IF (T .LE. 0.0) THEN
    !     write(21,'(/,''  THE INITIAL TEMPERATURE IS NOT POSITIVE. ''
    ! 1             /,''  RESET THE VARIABLE T. ''/)')
    !     IER = 3
    !     RETURN
    !  END IF


C  If the initial value is out of bounds, notify the user and return
C  to the calling routine.
C##### Test sur l'appartenance ou non des valeurs de départ à l'intervalle
C##### [LB(I), UB(I)].
      DO 30, I = 1, N
         IF ((X(I) .GT. UB(I)) .OR. (X(I) .LT. LB(I))) THEN
!            CALL PRT1
            IER = 2
            RETURN
         END IF
30    CONTINUE

C  Evaluate the function with input X and return value as F.
C##### appel de la subroutine qui va calculer un F à partir des fichiers
C##### Dnobs (yob) et Simptemp (ye).
	
	CALL OBJFUN(N,X,F)
	!CALL FCN(N,X,F)
!		write(*,*) "F après 1er appel ="
!		write(*,*) F
C  If the function is to be minimized, switch the sign of the function.
C  Note that all intermediate and final output switches the sign back
C  to eliminate any possible confusion for the user.
C##### Test si on a une maximisation, si oui, on change le signe de F.      
	IF(.NOT. MAX) F = -F
      NFCNEV = NFCNEV + 1 !ds ce cas, on compte une évaluation.
      FOPT = F !initialisation de fopt à f, point de départ.
      FSTAR(1) = F	!	@@@@@@@@@@@ ici, problème @@@@@@@@@@@@@@@@
      IF(IPRINT .GE. 1) CALL PRT2(MAX,N,X,F)

C  Start the main loop. Note that it terminates if (i) the algorithm
C  succesfully optimizes the function or (ii) there are too many
C  function evaluations (more than MAXEVL).
100   NUP = 0
      NREJ = 0
      NNEW = 0
      NDOWN = 0
      LNOBDS = 0
C--------------------------------------------------------------
      DO 400, M = 1, NT
         DO 300, J = 1, NS
            DO 200, H = 1, N

C  Generate XP, the trial value of X. Note use of VM to choose XP.
C##### Boucle qui va expliciter les XP à l'aide de la fonction aléatoire
C##### "RANMAR" ; on utilise aussi pour cela le vecteur VM.
      
    !  allocate(XP(N))

C---------------------------------------------------------------
	         DO 110, I = 1, N
                  IF (I .EQ. H) THEN
                     XP(I) = X(I) + (RANMAR()*2.- 1.) * VM(I)
                  ELSE
                     XP(I) = X(I)
                  END IF

C  If XP is out of bounds, select a point in bounds for the trial.
C##### Test sur l'appartenance du XP(I) à l'intervalle [LB(I), UB(I)].
C##### Si on sort de cet intervalle alors, on détermine un nouveau XP() de façon
C##### aléatoire et qui appartient à celui-ci.
                  IF((XP(I) .LT. LB(I)) .OR. (XP(I) .GT. UB(I))) THEN
                    XP(I) = LB(I) + (UB(I) - LB(I))*RANMAR()
				  ! cf. signification de çà @@@@@@@@@@@@@
				  LNOBDS = LNOBDS + 1
C##### On incrémente le nombre d'évaluations                    
				  NOBDS = NOBDS + 1
C##### Test sur le type d'affichage (3 = détails max)
C##### Voir l'initialisation et la définition du FP.!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                  IF(IPRINT .GE. 3) CALL PRT3(MAX,N,XP,X,FP,F)
                  END IF
110            CONTINUE


C  Evaluate the function with the trial point XP and return as FP.
C##### Evaluation de la fonction objectif, prend en entrée XP et ressort FP,
C##### ie une valeur intermédiaire de la fonction optimisée. 
               CALL OBJFUN(N,XP,FP)
		!	 CALL FCN(N,XP,FP)
C##### Test si maximisation, alors changement de signes
               IF(.NOT. MAX) FP = -FP
C##### Dans ce cas, on compte une évaluation de fonction supplémentaire.
               NFCNEV = NFCNEV + 1
C##### Test sur le type d'affichage et appel PRT4. 
   !            IF(IPRINT .GE. 3) CALL PRT4(MAX,N,XP,X,FP,F)
C  If too many function evaluations occur, terminate the algorithm.
C##### Test si le nombre total d'évaluations de fonction est > au nombre autorisé               
			 IF(NFCNEV .GE. MAXEVL) THEN
C##### Si tel est le cas, appel de la subroutine qui envoie un message d'erreur.
    !              CALL PRT5
C##### Test sur le type d'optimisation.
                  IF (.NOT. MAX) FOPT = -FOPT
                  IER = 1
                  RETURN
               END IF


C---------------------------------------------------------------------
C  Accept the new point if the function value increases.
               IF(FP .GE. F) THEN
   !               IF(IPRINT .GE. 3) THEN
   !                  write(21,'(''  POINT ACCEPTED'')')
   !               END IF
                  DO 120, I = 1, N
                     X(I) = XP(I)
120               CONTINUE
                  F = FP
                  NACC = NACC + 1
                  NACP(H) = NACP(H) + 1
                  NUP = NUP + 1

C  If greater than any other point, record as new optimum.
                  IF (FP .GT. FOPT) THEN
  !                IF(IPRINT .GE. 3) THEN
  !                   write(21,'(''  NEW OPTIMUM'')')
  !                END IF
                     DO 130, I = 1, N
                        XOPT(I) = XP(I)
130                  CONTINUE
                     FOPT = FP
                     NNEW = NNEW + 1
                  END IF

C---------------------------------------------------------------------
C  If the point is lower, use the Metropolis criteria to decide on
C  acceptance or rejection.
			 ELSE
C##### on calcule sa proba d'acceptation et on la compare à un aléatoire : 
                  P = EXPREP((FP - F)/T) ! proba d'acceptation 
									   ! entre 0 et 1 car (FP-F)<0	 
                  PP = RANMAR() ! aléatoire
                  IF (PP .LT. P) THEN
   !                  IF(IPRINT .GE. 3) CALL PRT6(MAX)
                     DO 140, I = 1, N
                        X(I) = XP(I)
140                  CONTINUE
                     F = FP
                     NACC = NACC + 1
                     NACP(H) = NACP(H) + 1
                     NDOWN = NDOWN + 1
                  ELSE
                     NREJ = NREJ + 1
   !                  IF(IPRINT .GE. 3) CALL PRT7(MAX)
                  END IF
               END IF

200         CONTINUE
300      CONTINUE
C------------------------------------------------------------------

C##### voir précisement le sens de çà et de VM.
C  Adjust VM so that approximately half of all evaluations are accepted.
         DO 310, I = 1, N
            RATIO = DFLOAT(NACP(I)) /DFLOAT(NS)
            IF (RATIO .GT. .6) THEN
               VM(I) = VM(I)*(1. + C(I)*(RATIO - .6)/.4)
            ELSE IF (RATIO .LT. .4) THEN
               VM(I) = VM(I)/(1. + C(I)*((.4 - RATIO)/.4))
            END IF
            IF (VM(I) .GT. (UB(I)-LB(I))) THEN
               VM(I) = UB(I) - LB(I)
            END IF
310      CONTINUE

    !     IF(IPRINT .GE. 2) THEN
    !        CALL PRT8(N,VM,XOPT,X)
    !     END IF

         DO 320, I = 1, N
            NACP(I) = 0
320      CONTINUE

400   CONTINUE
!	DEALLOCATE(XP)
   !   IF(IPRINT .GE. 1) THEN
   !      CALL PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
   !   END IF

C  Check termination criteria.
C##### Critères d'arrêt :
      QUIT = .FALSE.
      FSTAR(1) = F
C##### Critère d'arrêt vérifié, on est plus petit que eps, limite inferieure.
      IF ((FOPT - FSTAR(1)) .LE. EPS) QUIT = .TRUE.
      DO 410, I = 1, NEPS
C##### Critère d'arrêt non vérifié, continue.
         IF (ABS(F - FSTAR(I)) .GT. EPS) QUIT = .FALSE.
410   CONTINUE
      
	!	@@@@@@@@@@@@@@@@@@@ ici @@@@@@@@@@@@@@@@@@

C  Terminate SA if appropriate.
      IF (QUIT) THEN
         DO 420, I = 1, N
            X(I) = XOPT(I)
420      CONTINUE
         IER = 0
         IF (.NOT. MAX) FOPT = -FOPT
    !     IF(IPRINT .GE. 1) CALL PRT10
         RETURN
      END IF

C  If termination criteria is not met, prepare for another loop.
      T = RT*T
      DO 430, I = NEPS, 2, -1
         FSTAR(I) = FSTAR(I-1)
430   CONTINUE
      F = FOPT
      DO 440, I = 1, N
         X(I) = XOPT(I)
440   CONTINUE

!	write(*,*) "F à la fin de SA ="
!	write(*,*) F

C  Loop again.
      GO TO 100
 
      write(21,*) "Solutions :"
      DO i = 1, N
		write(21,*) xopt(i)
      END DO


      END
C########## fin de la subroutine SA



C#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c#########
      FUNCTION  EXPREP(RDUM)
C  This function replaces exp to avoid under- and overflows and is
C  designed for IBM 370 type machines. It may be necessary to modify
C  it for other machines. Note that the maximum and minimum values of
C  EXPREP are such that they has no effect on the algorithm.

      DOUBLE PRECISION  RDUM, EXPREP

      IF (RDUM .GT. 174.) THEN
         EXPREP = 3.69D+75
      ELSE IF (RDUM .LT. -180.) THEN
         EXPREP = 0.0
      ELSE
         EXPREP = EXP(RDUM)
      END IF

      RETURN
      END





c###########
      subroutine RMARIN(IJ,KL)
C  This subroutine and the next function generate random numbers. See
C  the comments for SA for more information. The only changes from the
C  orginal code is that (1) the test to make sure that RMARIN runs first
C  was taken out since SA assures that this is done (this test didn't
C  compile under IBM's VS Fortran) and (2) typing ivec as integer was
C  taken out since ivec isn't used. With these exceptions, all following
C  lines are original.

C This is the initialization routine for the random number generator
C     RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
      real U(97), C, CD, CM
      integer I97, J97
      common /raset1/ U, C, CD, CM, I97, J97
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.
     *    KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value
     *between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and
     *30081'
            stop
      endif
      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169)
      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         U(ii) = s
2     continue
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
      return
      end !de la subroutine rmarin






c#######################
      function ranmar()
      real U(97), C, CD, CM
      integer I97, J97
      common /raset1/ U, C, CD, CM, I97, J97
         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         uni = uni - C
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RANMAR = uni
      return
      END !de la fonction ranmar





C#####
      SUBROUTINE PRT1
C  This subroutine prints intermediate output, as does PRT2 through
C  PRT10. Note that if SA is minimizing the function, the sign of the
C  function value and the directions (up/down) are reversed in all
C  output to correspond with the actual function optimization. This
C  correction is because SA was written to maximize functions and
C  it minimizes by maximizing the negative a function.

      write(21,'(/,''  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS ''
     1          /,''  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY''
     2          /,''  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  ''
     3          /,''  LB(I) .LT. X(I) .LT. UB(I), I = 1, N. ''/)')

      RETURN
      END





C#####
      SUBROUTINE PRT2(MAX,N,X,F)

      DOUBLE PRECISION  X(*), F
      INTEGER  N
      LOGICAL  MAX

      write(21,'(''  '')')
  !    CALL PRTVEC(X,N,'INITIAL X')
      IF (MAX) THEN
         write(21,'(''  INITIAL F: '',/, G25.18)') F
      ELSE
         write(21,'(''  INITIAL F: '',/, G25.18)') -F
      END IF

      RETURN
      END




C##### Pourquoi il y a introduction du FP, on ne s'en sert pas !!

      SUBROUTINE PRT3(MAX,N,XP,X,FP,F)


      DOUBLE PRECISION  XP(*), X(*), FP, F
      INTEGER  N
      LOGICAL  MAX

      write(21,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
         write(21,'(''  CURRENT F: '',G25.18)') F
      ELSE
         write(21,'(''  CURRENT F: '',G25.18)') -F
      END IF
      CALL PRTVEC(XP,N,'TRIAL X')
      write(21,'(''  POINT REJECTED SINCE OUT OF BOUNDS'')')

      RETURN
      END





C#####
      SUBROUTINE PRT4(MAX,N,XP,X,FP,F)

      DOUBLE PRECISION  XP(*), X(*), FP, F
      INTEGER  N
      LOGICAL  MAX

      write(21,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
         write(21,'(''  CURRENT F: '',G25.18)') F
         CALL PRTVEC(XP,N,'TRIAL X')
         write(21,'(''  RESULTING F: '',G25.18)') FP
      ELSE
         write(21,'(''  CURRENT F: '',G25.18)') -F
         CALL PRTVEC(XP,N,'TRIAL X')
         write(21,'(''  RESULTING F: '',G25.18)') -FP
      END IF

      RETURN
      END





C#####
      SUBROUTINE PRT5

      write(21,'(/,''  TOO MANY FUNCTION EVALUATIONS; CONSIDER ''
     1          /,''  INCREASING MAXEVL OR EPS, OR DECREASING ''
     2          /,''  NT OR RT. THESE RESULTS ARE LIKELY TO BE ''
     3          /,''  POOR.'',/)')

      RETURN
      END



C#####
      SUBROUTINE PRT6(MAX)

      LOGICAL  MAX

      IF (MAX) THEN
         write(21,'(''  THOUGH LOWER, POINT ACCEPTED'')')
      ELSE
         write(21,'(''  THOUGH HIGHER, POINT ACCEPTED'')')
      END IF

      RETURN
      END




C#####
      SUBROUTINE PRT7(MAX)

      LOGICAL  MAX

      IF (MAX) THEN
         write(21,'(''  LOWER POINT REJECTED'')')
      ELSE
         write(21,'(''  HIGHER POINT REJECTED'')')
      END IF

      RETURN
      END



C#####
      SUBROUTINE PRT8(N,VM,XOPT,X)

      DOUBLE PRECISION  VM(*), XOPT(*), X(*)
      INTEGER  N

      write(21,'(/,
     1  '' INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT'',/)')
      CALL PRTVEC(VM,N,'NEW STEP LENGTH (VM)')
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(X,N,'CURRENT X')
      write(21,'('' '')')

      RETURN
      END



C#####
      SUBROUTINE PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)

      DOUBLE PRECISION  XOPT(*), VM(*), T, FOPT
      INTEGER  N, NUP, NDOWN, NREJ, LNOBDS, NNEW, TOTMOV
      LOGICAL  MAX

      TOTMOV = NUP + NDOWN + NREJ

      write(21,'(/,
     1  '' INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION'',/)')
      write(21,'(''  CURRENT TEMPERATURE:            '',G12.5)') T
      IF (MAX) THEN
         write(21,'(''  MAX FUNCTION VALUE SO FAR:  '',G25.18)') FOPT
         write(21,'(''  TOTAL MOVES:                '',I8)') TOTMOV
         write(21,'(''     UPHILL:                  '',I8)') NUP
         write(21,'(''     ACCEPTED DOWNHILL:       '',I8)') NDOWN
         write(21,'(''     REJECTED DOWNHILL:       '',I8)') NREJ
         write(21,'(''  OUT OF BOUNDS TRIALS:       '',I8)') LNOBDS
         write(21,'(''  NEW MAXIMA THIS TEMPERATURE:'',I8)') NNEW
      ELSE
         write(21,'(''  MIN FUNCTION VALUE SO FAR:  '',G25.18)') -FOPT
         write(21,'(''  TOTAL MOVES:                '',I8)') TOTMOV
         write(21,'(''     DOWNHILL:                '',I8)')  NUP
         write(21,'(''     ACCEPTED UPHILL:         '',I8)')  NDOWN
         write(21,'(''     REJECTED UPHILL:         '',I8)')  NREJ
         write(21,'(''  TRIALS OUT OF BOUNDS:       '',I8)')  LNOBDS
         write(21,'(''  NEW MINIMA THIS TEMPERATURE:'',I8)')  NNEW
      END IF
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(VM,N,'STEP LENGTH (VM)')
      write(21,'('' '')')

      RETURN
      END



C#####
      SUBROUTINE PRT10

      write(21,'(/,'' SA ACHIEVED TERMINATION CRITERIA. IER = 0. '',/)')

      RETURN
      END



C#####
      SUBROUTINE PRTVEC(VECTOR,NCOLS,NAME)
C  This subroutine prints the double precision vector named VECTOR.
C  Elements 1 thru NCOLS will be printed. NAME is a character variable
C  that describes VECTOR. Note that if NAME is given in the call to
C  PRTVEC, it must be enclosed in quotes. If there are more than 10
C  elements in VECTOR, 10 elements will be printed on each line.

      INTEGER NCOLS
      DOUBLE PRECISION VECTOR(NCOLS)
      CHARACTER *(*) NAME

      write(21,1001) NAME !old

      IF (NCOLS .GT. 10) THEN
         LINES = INT(NCOLS/10.)

         DO 100, I = 1, LINES
            LL = 10*(I - 1)
            write(21,1000) (VECTOR(J),J = 1+LL, 10+LL)
  100    CONTINUE

         write(21,1000) (VECTOR(J),J = 11+LL, NCOLS)
      ELSE
         write(21,1000) (VECTOR(J),J = 1, NCOLS)
      END IF

 1000 FORMAT( 10(G12.5,1X))

 1001 FORMAT(/,25X,A)

      RETURN
      END




