      subroutine sft(x,y,n,is)
c     schnelle fouriertransformation der 2**n komplexen werte (x,y)
c     is negativ - in den frequenzbereich / positiv - in den zeitbereich
c     normierung - abs(is) =1 - ohne / 2 - mit 1/ng / 3 - mit 1/sqrt(ng)
c                            4 - ohne normierung,ohne kontrollausdruck
      implicit double precision(a-h,o-z)
	implicit integer(i-n)
      dimension x(7),y(7)
c     tabelle der zweierpotenzen
      integer   zh(21)
      zh(1)=1
      do 1 l=1,n
    1 zh(l+1)=2*zh(l)
      ng=zh(n+1)
      gn=1./float(ng)
      piz=8.*datan(1.d0)
c     kernprogramm.dreifache schleife ueber schritt/index/teilserie
      do 2 m=1,n
      nar=zh(m)
      lar=ng/nar
      larh=lar/2
      alpha =  piz/dble(isign(lar,is))
      do 3 jr=1,larh
      beta=alpha*dble(jr-1)
      excos = dcos(beta)
      exsin = dsin(beta)
      ja=jr-lar
      do 4 nr=1,nar
      ja=ja+lar
      jb=ja+larh
      zx = x(ja)-x(jb)
      zy = y(ja)-y(jb)
      x(ja) = x(ja)+x(jb)
      y(ja) = y(ja)+y(jb)
      x(jb) = zx*excos-zy*exsin
      y(jb) = zx*exsin+zy*excos
    4 continue
    3 continue
c     normierung
    2 continue
      if(iabs(is).eq.1.or.iabs(is).eq.4) go to 10
      if(iabs(is).eq.3) gn=dsqrt(gn)
    5 do 6  l=1,ng
      y(l) = y(l)*gn
    6 x(l)=x(l)*gn
c     umordnung nach @bitreversed@ indizes
   10 do 7 j=1,ng
      js=j-1
      k=1
      nny=n+1
      do 8 ny=1,n
      nny=nny-1
      if(js.lt.zh(nny)) go to 8
      js=js-zh(nny)
      k=k+zh(ny)
    8 continue
      if(j-k) 9,7,7
    9 zx = x(j)
      zy = y(j)
      x(j) = x(k)
      y(j) = y(k)
      x(k) = zx
      y(k) = zy
    7 continue
      if(iabs(is).eq.4) return
      return
      end

