c------------------------------------------------
c  computes for real number r integer number ceil
c  such that ceil >= r.
c------------------------------------------------
	integer function ceil(r)
	real r
	integer ir
c
	ir=int(r)
	if(real(ir).eq.r) then
		ceil=ir
	else
		ceil=ir+1
	endif
c
	return
	end
