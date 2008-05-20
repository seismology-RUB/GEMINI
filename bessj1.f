      FUNCTION bessj1(x)
      REAL bessj1,x
      REAL ax,xx,z
      REAL p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.e0,-7895059235.e0,
     *242396853.1e0,-2972611.439e0,15704.48260e0,-30.16036606e0/,s1,s2,
     *s3,s4,s5,s6/144725228442.e0,2300535178.e0,18583304.74e0,
     *99447.43394e0,376.9991397e0,1.e0/
      DATA p1,p2,p3,p4,p5/1.e0,.183105e-2,-.3516396496e-4,
     *.2457520174e-5,-.240337019e-6/, q1,q2,q3,q4,q5/.04687499995e0,
     *-.2002690873e-3,.8449199096e-5,-.88228987e-6,.105787412e-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.

