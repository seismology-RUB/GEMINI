      FUNCTION bessj0(x)
      REAL bessj0,x
      REAL ax,xx,z
      REAL p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.e0,-.1098628627e-2,.2734510407e-4,
     *-.2073370639e-5,.2093887211e-6/, q1,q2,q3,q4,q5/-.1562499995e-1,
     *.1430488765e-3,-.6911147651e-5,.7621095161e-6,-.934945152e-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.e0,-13362590354.e0,
     *651619640.7e0,-11214424.18e0,77392.33017e0,-184.9052456e0/,s1,s2,
     *s3,s4,s5,s6/57568490411.e0,1029532985.e0,9494680.718e0,
     *59272.64853e0,267.8532712e0,1.e0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.

