c-----------------------------------------------------------------------------------
c add a mode to the total displacement for given frequency and receiver
c-----------------------------------------------------------------------------------
	subroutine add_tor_displ(i,f,p,cg,qatt,psi,rearth,wr,zdis)
	include 'sardim.h'
	include 'sar.h'
	real f,p,cg,qatt,rearth,wr,azi,saz,caz
	complex df_pot,dt_pot,zspl,zspt,zspn,zspe
	integer i
	complex psi(0:2),zdis

	if(sar_comp(i).eq.'Z' .or. sar_comp(i).eq.'R' .or. sar_comp(i).eq.'H') then
		return
	else if(sar_comp(i).eq.'L') then
		call potential(sar_dis(i),sar_phi(i),f,p,cg,qatt,psi,rearth,df_pot,3)
		zdis=zdis-wr*df_pot
	else if(sar_comp(i).eq.'T') then
		call potential(sar_dis(i),sar_phi(i),f,p,cg,qatt,psi,rearth,dt_pot,2)
		zdis=zdis+wr*dt_pot
	else if(sar_comp(i).eq.'N' .or.sar_comp(i).eq.'E') then
		call potential(sar_dis(i),sar_phi(i),f,p,cg,qatt,psi,rearth,dt_pot,2)
		call potential(sar_dis(i),sar_phi(i),f,p,cg,qatt,psi,rearth,df_pot,3)
		zspl = -wr*df_pot
		zspt = +wr*dt_pot
		call sar_propdir(sar_xsta(i),sar_ysta(i),azi)
		caz=cos(azi)
		saz=sin(azi)
		if(sar_csys.eq.'S') then
			zspn=-zspl*caz+zspt*saz
			zspe=+zspl*saz+zspt*caz
		else
			zspn=+zspl*caz-zspt*saz
			zspe=+zspl*saz+zspt*caz
		endif
		if(sar_comp(i).eq.'N') zdis=zdis+zspn
		if(sar_comp(i).eq.'E') zdis=zdis+zspe
	endif
	return
	end


