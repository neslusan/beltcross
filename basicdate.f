c Calculation of the date when the asteroid crosses the meteoroid stream
c Version: September 6, 2022

      implicit integer (i,j)
      implicit real*8 (a-h,k-z)
      character*15 name,namx
      character*80 txt
      character*110 ret110
      dimension isol(1500),iau(1500),jnn(1500)
      dimension qs(1500),es(1500),sos(1500),gos(1500),sks(1500)

      pi=4.0d0*datan(1.0d0)
      pi2=2.0d0*pi
      pi180=pi/1.80d2
      gauss=1.7202098950d-2

      open(unit=10,file='allshowers.d',access='sequential')
      open(unit=11,file='object.dat',access='sequential')
      open(unit=12,file='relativeV.dw',access='sequential')
      open(unit=15,file='basicdate.dw',access='sequential')
      read(10,*) ret110

      j=0
  10  continue
      j=j+1
      read(10,*,end=30) ipor,iau(j),isol(j),ls,al,dl,vg,a,qs(j),es(j),so
     *s(j),gos(j),sks(j),jnn(j)
      if(sos(j).lt.-0.5d0.or.gos(j).lt.-0.5d0.or.sks(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(dabs(a+1.0d0).lt.1.0d-5.and.qs(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(dabs(a+1.0d0).lt.1.0d-5.and.es(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(qs(j).lt.-0.5d0.and.es(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(qs(j).lt.-0.5d0) then
        qs(j)=a*(1.0d0-es(j))
        goto 10
      endif
      if(es(j).lt.-0.5d0) es(j)=1.0d0-qs(j)/a
      goto 10
  30  continue
      jsh=j-1
      close(unit=10)
      
  50  continue
      read(11,*) txt
      read(11,*) name
      read(11,*) txt
      read(11,*) a
      read(11,*) txt
      read(11,*) ecc
      read(11,*) txt
      read(11,*) som
      read(11,*) txt
      read(11,*) gom
      read(11,*) txt
      read(11,*) sk
      read(11,*) txt
      read(11,*) irk
      read(11,*) txt
      read(11,*) imes
      read(11,*) txt
      read(11,*) day
      read(11,*) txt
      read(11,*) mmepo
      read(11,*) txt
      read(11,*) jyear
      read(11,*) txt
      read(11,*) rcrit
      sdp=gauss/a/dsqrt(a)

      call julian(irk,imes,day,t19,t20)
      tjd=3.6525d4*t20+2.451545d6
  70  continue
      jast=j-1
      close(unit=11)

  90  continue
      read(12,*,end=120) namx,fa,dr,rast,vrel,angl
      read(12,*) iaux,isx,jnnw,iarc
      iok=0
      do ii=1,jsh
        if(isx.eq.isol(ii).and.iaux.eq.iau(ii)) then
          iok=1
          qx=qs(ii)
          ex=es(ii)
          sox=sos(ii)
          gox=gos(ii)
          skx=sks(ii)
        endif
      enddo
      if(iok.eq.0) then
        write(*,*) 'Unknown ERROR (shower was not identified).'
        write(*,*) 'Program terminated.'
        stop
      endif
      if(dr.le.rcrit) then
        tt=tjd-mmepo*pi180/sdp
        arg=dsqrt((1.0d0-ecc)/(1.0d0+ecc))*dtan(fa*pi180/2.0d0)
        ee=2.0d0*datan(arg)
        mm=ee-ecc*dsin(ee)
        tcros=mm/sdp+tt
        pp=pi2/sdp
      write(15,300) namx,iaux,isx,jnnw,fa,tcros,pp,iarc,dr,qx,ex,sox,gox
     *,skx,rast,vrel,angl
      endif
      goto 90
 120  continue
      close(unit=12)
      close(unit=15)
      stop

 300  format(a15,i5,i4,i6,f9.3,2e14.7,i4,3f8.4,3f9.3,2f7.2,f8.2)
      end
C --------------------------------------------------------------------
      SUBROUTINE JULIAN(IROK,IMES,DEN,T19,T20)
C --------------------------------------------------------------------
C     - CALCULATION OF JULIAN DATE FROM CIVIL DATE
C
C INPUT:
C    IROK, IMES, DEN - YEAR, MONTH, AND DAY OF THE DATE
C
C OUTPUT:
C    T19 - GIVEN TIME ACCOUNTED FROM THE BEGINNING OF EPOCH 1900.0
C          (JD = 2415020.0) AND EXPRESSED IN JULIAN CENTURIES
C    T20 - GIVEN TIME ACCOUNTED FROM THE BEGINNING OF EPOCH 2000.0
C          (JD = 2451545.0) AND EXPRESSED IN JULIAN CENTURIES
C
C    JULIAN DATE CAN BE CALCULATED FROM "T20" AS:
C    JD=36525.0*T20+2451545.0
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (A-H,K-Z)
      IMPLICIT INTEGER (I,J)
      DIMENSION JMESI(13)
      DATA JMESI/0,31,28,31,30,31,30,31,31,30,31,30,31/
C
      IF(IROK.GT.1582) GOTO 101
      IF(IROK.EQ.1582) GOTO 102
      IF(IROK.GT.0.AND.IROK.LT.1582) GOTO 103
      IF(IROK.GE.-4713.AND.IROK.LT.0) GOTO 104
C
      T20=(IROK+4713)*365
      T20=T20-INT((-IROK-4712)/4)
  100 T0JD=-0.5D0
      GOTO 210
  101 T20=(IROK-1583)*365
      T20=T20+INT((IROK-1581)/4)
      T20=T20-INT((IROK-1501)/100)+INT((IROK-1201)/400)
      T0JD=2.2992385D6
      GOTO 210
  102 T20=0.0D0
      IF(IMES.EQ.10.AND.DEN.GT.1.49999999D1) T20=-1.0D1
      IF(IMES.GT.10) T20=-1.0D1
      T0JD=2.2988835D6
      GOTO 210
  103 T20=(IROK-1)*365
      T20=T20+INT((IROK-1)/4)
      T0JD=1.7214235D6
      GOTO 210
  104 T20=(IROK+4713)*365
      T20=T20+INT((IROK+4715)/4)
      T0JD=-0.5D0
  210 YY=IROK/1.00D2-INT(IROK/100)
      ZZ=IROK/4.00D2-INT(IROK/400)
      IF(IROK.GT.1582.AND.YY.EQ.0.0D0.AND.ZZ.GT.0.0D0) GOTO 220
      YY=IROK/4.0D0-INT(IROK/4)
      IF(YY.EQ.0.0D0) JMESI(3)=29
      IF(IROK.EQ.-1) JMESI(3)=29
  220 CONTINUE
      DO 225 JJM=1,IMES-1
      T20=T20+JMESI(JJM+1)
  225 CONTINUE
      T20=T20+DEN-1.0D0+T0JD
      JMESI(3)=28
C
      T19=(T20-2.415020D6)/3.6525D4
      T20=(T20-2.451545D6)/3.6525D4
      RETURN
C
      END
