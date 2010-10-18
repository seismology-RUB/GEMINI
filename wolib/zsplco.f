c----------------------------------------------------------------
c  $Id: zsplco.f,v 1.2 2003/04/03 12:09:38 wolle Exp $
c
c  Spline complex function values by splining real and imaginary
c  part separately. Calls dsplco.f: see comments there.
c  work: workspace of dimension n
c
	subroutine zsplco(x,zy,n,zyp1,zypn,zy2,work,work2)
	integer n,i
	double precision x(n),work(n),work2(n)
	double complex zy(n),zy2(n),zyp1,zypn,zi
c 
	zi=dcmplx(0.d0,1.d0)
	do i=1,n
		work(i)=real(zy(i))
c		print *, work(i)
	enddo
	call dsplco(x,work,n,real(zyp1),real(zypn),work2)
	do i=1,n
		zy2(i)=dcmplx(work2(i),0.d0)
c		print *,work2(i)
	enddo
	do i=1,n
		work(i)=dimag(zy(i))
c		print *,work(i)
	enddo
	call dsplco(x,work,n,dimag(zyp1),dimag(zypn),work2)
	do i=1,n
		zy2(i)=zy2(i)+zi*work2(i)
	enddo
	return
	end

