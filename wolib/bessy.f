      FUNCTION bessy(n,x)
      INTEGER n
      REAL bessy,x
CU    USES bessy0,bessy1
      INTEGER j
      REAL by,bym,byp,tox,bessy0,bessy1
      if(n.lt.2)pause 'bad argument n in bessy'
      tox=2./x
      by=bessy1(x)
      bym=bessy0(x)
      do 11 j=1,n-1
        byp=j*tox*by-bym
        bym=by
        by=byp
11    continue
      bessy=by
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.
