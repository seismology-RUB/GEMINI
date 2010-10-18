c---------------------------------------------------------
c  trim a string from left, return index k of first
c  non--blank character, returns k=len(string)+1 if 
c  string is all blanks
c---------------------------------------------------------
	subroutine trimleft(string,k)
	character string*(*)
	integer k,nc
c
	k=1
	nc=len(string)
 1	if(string(k:k).ne.' ') return
	if(k.lt.nc) then
		k=k+1
		goto 1
	else
		k=k+1
	endif
	return
	end
c----------------------------------------------------------
c  trim a string from right and return index k of last
c  non--blank chracter, returns k=0 if string is all blanks
c----------------------------------------------------------
	subroutine trimright(string,k)
	character string*(*)
	integer k
c
	k=len(string)
 1	if(string(k:k).ne.' ') return
	if(k.gt.1) then
		k=k-1
		goto 1
	else
		k=k-1
	endif
	return
	end
