c-------------------------------------------------------
c  $Id: splitstring.f,v 1.2 2003/10/25 12:16:30 wolle Exp $
c-------------------------------------------------------
c  Take a character string and split it into
c  an array of floats using <splitc> as splitting
c  delimiting character.
c-------------------------------------------------------
	subroutine fsplitstring(string,splitc,n,f)
	integer n,ich1,ich2,i,j,ia
	real f(n)
	character string*(*),splitc*1
c
	call trimleft(string,ich1)
	call trimright(string,ich2)
	if(ich2.eq.0) then
		write(0,*) '<fsplitstring>: empty string !!'
		stop
	endif
	n=0
	i=ich1
 2	continue
		j=0
		ia=i
 1		continue
			if(i.gt.ich2) then
				n=n+1
				read(string(ia:(ia+j-1)),*) f(n)
				goto 3
			endif
			if(string(i:i).ne.splitc) then
				j=j+1
				i=i+1
				goto 1
			else if(string(i:i).eq.splitc.and.j.gt.0) then
				n=n+1
				read(string(ia:(ia+j-1)),*) f(n)
				i=i+1
				goto 2
			else
				i=i+1
				goto 2
			endif
 3	continue
	return
	end
c-------------------------------------------------------
c  Take a character string and split it into
c  an array of integers using <splitc> as splitting
c  delimiting character.
c-------------------------------------------------------
	subroutine isplitstring(string,splitc,n,k)
	integer n,ich1,ich2,i,j,ia
	integer k(n)
	character string*(*),splitc*1
c
	call trimleft(string,ich1)
	call trimright(string,ich2)
	if(ich2.eq.0) then
		print *,'<fsplitstring>: empty string !!'
		stop
	endif
	n=0
	i=ich1
 2	continue
		j=0
		ia=i
 1		continue
			if(i.gt.ich2) then
				n=n+1
				read(string(ia:(ia+j-1)),*) k(n)
				goto 3
			endif
			if(string(i:i).ne.splitc) then
				j=j+1
				i=i+1
				goto 1
			else if(string(i:i).eq.splitc.and.j.gt.0) then
				n=n+1
				read(string(ia:(ia+j-1)),*) k(n)
				i=i+1
				goto 2
			else
				i=i+1
				goto 2
			endif
 3	continue
	return
	end
c-------------------------------------------------------
c  Take a character string and split it into
c  an array of words using <splitc> as splitting
c  delimiting character.
c-------------------------------------------------------
	subroutine csplitstring(string,splitc,n,word)
	integer n,ich1,ich2,i,j,ia
	character*(*) word(n)
	character string*(*),splitc*1
c
	call trimleft(string,ich1)
	call trimright(string,ich2)
	print *,string(ich1:ich2)
	if(ich2.eq.0) then
		print *,'<csplitstring>: empty string !!'
		stop
	endif
	n=0
	i=ich1
 2	continue
		j=0
		ia=i
 1		continue
			if(i.gt.ich2) then
				n=n+1
				word(n)=string(ia:(ia+j-1))
				goto 3
			endif
			if(string(i:i).ne.splitc) then
				j=j+1
				i=i+1
				goto 1
			else if(string(i:i).eq.splitc.and.j.gt.0) then
				n=n+1
				word(n)=string(ia:(ia+j-1))
				i=i+1
				goto 2
			else
				i=i+1
				goto 2
			endif
 3	continue
	return
	end
