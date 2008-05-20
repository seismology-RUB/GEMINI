c--------------------------------------------------------------
c  $Id: isoptbelpar.f,v 1.5 2003/05/30 06:37:14 wolle Exp $
c
c  Routine which for given parameter index
c  computes change of complex elastic parameters
c  at fixed r.
c  Assumes that layer_nlactive is correctly set.
c  Assumes that model is isotropic !
c  Assumes that perturbations are in solid medium !
c
c  dispersion correction was not OK 
c    - zkap,zmu are already dispersion corrected
c    - real(zkap-4/3 zmue)/ro ergab nicht das alfa aus der Tabelle !
c    - zkap=kappa*cmplx(1,qk)*dexp( 2*atan(qk)/pi*dlog(fratio) )
c    - um alfa zu bekommen muss der dexp(..)-Term rausdividiert werden
c----------------------------------------------------------------------
	subroutine isoptbelpar(ipar,ro,zkap,zmu,ron,zpan,zpcn,zpfn,zpln,zpnn,ptb)
	include 'pi.h'
	include 'flnm.h'
	integer ipar
	double complex zkap,zmu,zom
	double complex zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun
	double precision ro,ron,ptb,del,om,fratio
	double precision alf,bet,qk,qm,eta,kappa,mue
	common/omega/zom,om
c
	if(flnm_seldamp.eq.1) then
		fratio = 1.d0
	else if(flnm_seldamp.eq.3) then
		fratio=om/(2.*pi)/flnm_fref
	else
		print *,'illegal value of flnm_seldamp'
		stop
	endif
c
c  get tabulated velocities alf and bet, qk and qm
c
	qk=dimag(zkap)/real(zkap)
	qm=dimag(zmu)/real(zmu)
	kappa=real(zkap)/dexp( 2.*datan(qk)/pi*dlog(fratio) )
	mue=real(zmu)/dexp( 2.*datan(qm)/pi*dlog(fratio) )
	alf=dsqrt( (kappa+4.d0*mue/3.d0)/ro )
	bet=dsqrt( mue/ro )
	eta=1.d0
c
c  apply perturbation
c
	ron=ro
	if(ipar.eq.1) then
		del=ro*ptb
		ron=ro+del
		call zelcon(fratio,ro+del,alf,alf,bet,bet,eta,qk,qm,zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun)
	else if(ipar.eq.2) then
		del=alf*ptb
		call zelcon(fratio,ro,alf+del,alf+del,bet,bet,eta,qk,qm,zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun)
	else if(ipar.eq.3) then
		del=bet*ptb
		call zelcon(fratio,ro,alf,alf,bet+del,bet+del,eta,qk,qm,zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun)
	else if(ipar.eq.4) then
		del=qk*ptb
		call zelcon(fratio,ro,alf,alf,bet,bet,eta,qk+del,qm,zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun)
	else if(ipar.eq.5) then
		del=qm*ptb
		call zelcon(fratio,ro,alf,alf,bet,bet,eta,qk,qm+del,zpan,zpcn,zpfn,zpln,zpnn,zkapn,zmun)
	endif
c
	return
	end
