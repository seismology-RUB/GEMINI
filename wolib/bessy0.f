      FUNCTION bessy0(x)
      REAL bessy0,x
CU    USES bessj0
      REAL xx,z,bessj0
      REAL p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.e0,-.1098628627e-2,.2734510407e-4,
     *-.2073370639e-5,.2093887211e-6/, q1,q2,q3,q4,q5/-.1562499995e-1,
     *.1430488765e-3,-.6911147651e-5,.7621095161e-6,-.934945152e-7/
      DATA r1,r2,r3,r4,r5,r6/-2957821389.e0,7062834065.e0,
     *-512359803.6e0,10879881.29e0,-86327.92757e0,228.4622733e0/,s1,s2,
     *s3,s4,s5,s6/40076544269.e0,745249964.8e0,7189466.438e0,
     *47447.26470e0,226.1030244e0,1.e0/
      if(x.lt.8.)then
        y=x**2
        bessy0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))+.636619772*bessj0(x)*log(x)
      else
        z=8./x
        y=z**2
        xx=x-.785398164
        bessy0=sqrt(.636619772/x)*(sin(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))+z*cos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.

