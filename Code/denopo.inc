c      double precision theta(100),tempe(100),nconc(100) 
      character varopt(12)*15,paropt(12)*15, nvobst(500)*15  !JFM130116
	  integer nvaropt, nobs   !JFM130116
	  real varoptpds(12), dnobst(500), simtempt(500)
	  common /varcalee/varopt,paropt,nvobst
	  common /varval2/nvaropt,varoptpds,dnobst,simtempt,nobs
c	  common /forcage/ theta,tempe,nconc 
c	  common /sol/ por
	