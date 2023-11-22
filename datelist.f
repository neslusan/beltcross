c Creation of the list of passages of asteroids through the meteoroid
c   streams
c Version: September 6, 2022

      implicit integer (i,j)
      implicit real*8 (a-h,k-z)
      character*1 cf(10)
      character*4 month(12)
      character*11 arc
      character*15 name
      character*16 filnam
      character*80 txt
      
      cf(1)='0'
      cf(2)='1'
      cf(3)='2'
      cf(4)='3'
      cf(5)='4'
      cf(6)='5'
      cf(7)='6'
      cf(8)='7'
      cf(9)='8'
      cf(10)='9'

      month(1)='Jan.'
      month(2)='Feb.'
      month(3)='Mar.'
      month(4)='Apr.'
      month(5)='May '
      month(6)='June'
      month(7)='July'
      month(8)='Aug.'
      month(9)='Sep.'
      month(10)='Oct.'
      month(11)='Nov.'
      month(12)='Dec.'
      
      open(unit=10,file='basicdate.dw',access='sequential')
      open(unit=11,file='object.dat',access='sequential')

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
      read(11,*) mm
      read(11,*) txt
      read(11,*) jyear
      read(11,*) txt
      read(11,*) rcrit
      open(unit=16,file='year.d',access='sequential')
      write(16,320) jyear
      close(unit=16)

      it1=1
      it2=1
      it3=1
      it4=1
      if(jyear.lt.1.or.jyear.gt.9999) then
      write(*,*) 'The calculation can be done only for interval of years
     *'
        write(*,*) 'from 1 to 9999.'
        write(*,*) 'Program terminated.'
        stop
      endif
      do jj=1,jyear
        it1=it1+1
        if(it1.lt.11) goto 20
        it1=1
        it2=it2+1
        if(it2.lt.11) goto 20
        it2=1
        it3=it3+1
        if(it3.lt.11) goto 20
        it3=1
        it4=it4+1
        if(it4.lt.11) goto 20
          write(*,*) 'Unknown error; program terminated.'
          stop
  20  continue
      enddo
      filnam='unarranged'//cf(it4)//cf(it3)//cf(it2)//cf(it1)//'.d'
      open(unit=15,file=filnam,access='sequential')
      write(15,330) name
      write(15,340)

  30  continue
      read(10,300,end=50) name,iaux,isx,jnn,f,tcros,pp,iarc,dr,qx,ex,sox
     *,gox,skx,rast,vrel,angl
      if(iarc.eq.1) then
        arc=' postperih.'
      else
        arc=' preperih.'
      endif
      tjd0=(tcros-2.451545d6)/3.6525d4
      call INVJUL(tjd0,irk0,ims,den)
      if(jyear.eq.irk0) then
        dni=tcros
        goto 40
      endif
      if(jyear.lt.irk0) then
        jjgo=-1
        j=-1
        dni=tcros-pp
      else
        jjgo=1
        j=0
        dni=tcros
      endif
  40  continue
        tjd=(dni-2.451545d6)/3.6525d4
        call INVJUL(tjd,irk,ims,den)
  45  continue
        if(irk.eq.jyear) then
      if(iaux.le.9) write(15,350) iaux,isx,jnn,arc,dr,rast,vrel,angl,irk
     *,ims,den
      if(iaux.ge.10.and.iaux.le.99) write(15,360) iaux,isx,jnn,arc,dr,ra
     *st,vrel,angl,irk,ims,den
      if(iaux.ge.100.and.iaux.le.999) write(15,370) iaux,isx,jnn,arc,dr,
     *rast,vrel,angl,irk,ims,den
      if(iaux.ge.1000) write(15,380) iaux,isx,jnn,arc,dr,rast,vrel,angl,
     *irk,ims,den
        else
          j=j+jjgo
          dni=tcros+j*pp
          if(dni.le.1.7215235d6.or.dni.gt.3.1823945d6) goto 30
          goto 40
        endif
      goto 30
  50  continue
      close(unit=10)
      close(unit=15)
      stop
  
 300  format(a15,i5,i4,i6,f9.3,2e14.7,i4,3f8.4,3f9.3,2f7.2,f8.2)
 320  format(i5)
 330  format('OBJECT NAME: ',a15,/) 
 340  format('IAUNo. Sol.  n      arc        MOID   r_obj  v_rel   angle
     *        date')
 350  format('#000',i1,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 360  format('#00',i2,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 370  format('#0',i3,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 380  format('#',i4,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)

      end
C ---------------------------------------------------------------------
      SUBROUTINE INVJUL(TJD,IRK,IMS,DNI)
C ---------------------------------------------------------------------
C     - CALCULATION OF CIVIL DATE FROM JULIAN DATE (INVERSION OPERATION
C       TO THAT EXECUTED BY SUBROUTINE "JULIAN")
C
C INPUT:
C    TJD - TIME FROM THE BEGINNING OF EPOCH 2000.0 (JD = 2451545.0)
C          IN JULIAN CENTURIES (IT CAN BE CALCULATED FROM JULIAN DATE
C          "JD" IN DAYS AS: TJD=(JD-2451545.0)/36525.0)
C
C OUTPUT:
C    IRK, IMS, DNI - YEAR, MONTH, AND DAY OF THE CORRESPONDING CIVIL
C                    DATE
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IMPLICIT REAL*8 (D,T)
      IMPLICIT INTEGER (I)
      DIMENSION IMESI(13)
      DATA IMESI/0,31,28,31,30,31,30,31,31,30,31,30,31/
C
      DNI=TJD*3.6525D4+2.451545D6+1.0D0
      IF(DNI.LT.0.5D0) GOTO 16
      IF(DNI.GE.0.5D0.AND.DNI.LT.1.7214245D6) GOTO 15
      IF(DNI.GE.1.7214245D6.AND.DNI.LT.2.2991615D6) GOTO 14
      IF(DNI.GE.2.2991615D6.AND.DNI.LT.2.2992395D6) GOTO 13
      IF(DNI.GE.2.2992395D6.AND.DNI.LT.2.4515455D6) GOTO 12
C
C  FROM 2000-JAN-1.0
  11  CONTINUE
      DNI=DNI-2.4515445D6
      IRK=INT(DNI/3.65242198D2)
      DN0=365*IRK+INT((IRK+3)/4.0D0)-INT((IRK-1)/1.0D2)+INT((IRK-1)/4.0D
     *2)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 30
      IRK=IRK+2000
      GOTO 20
  30  IRK=IRK-1
      DNI=DNI+3.65D2
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DH100=IRK/1.0D2-INT(IRK/1.0D2)
      IF(DABS(DH4).LT.1.0D-6.AND.DH100.GT.0.0D0) DNI=DNI+1.0D0
      DH400=IRK/4.0D2-INT(IRK/4.0D2)
      IF(DABS(DH400).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK+2000
      GOTO 20
C
C   FROM 1583-JAN-1.0 TO 2000-JAN-1.0
  12  CONTINUE 
      DNI=DNI-2.2992385D6
      IRK=INT(DNI/3.65242198D2)
      DN0=365*IRK+INT((IRK+2)/4.0D0)-INT((IRK+82)/1.0D2)+INT((IRK+382)/4
     *.0D2)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 35
      IRK=IRK+1583
      GOTO 20
  35  IRK=IRK-1
      DH4=(IRK+3)/4.0D0-INT((IRK+3)/4.0D0)
      DH100=(IRK+83)/1.0D2-INT((IRK+83)/1.0D2)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6.AND.DH100.GT.0.0D0) DNI=DNI+1.0D0
      DH400=(IRK+383)/4.0D2-INT((IRK+383)/4.0D2)
      IF(DABS(DH400).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK+1583
      GOTO 20
C
C   FROM 1582-OCT-15.0 TO 1583-JAN-1.0
  13  CONTINUE 
      IRK=1582
      DNI=DNI-2.2991605D6+1.4D1
      IF(DNI.GE.3.2D1) GOTO 40
      IMS=10
      GOTO 26
  40  CONTINUE
      IF(DNI.GE.6.2D1) GOTO 45
      DNI=DNI-3.1D1
      IMS=11
      GOTO 26
  45  DNI=DNI-6.1D1
      IMS=12
      GOTO 26
C
C   FROM 1-JAN-1.0 TO 1582-OCT-15.0
  14  CONTINUE
      DNI=DNI-1.7214235D6
      IRK=INT(DNI/3.6525D2)+1
      DN0=3.65D2*(IRK-1)+INT((IRK-1)/4.0D0)
      DNI=DNI-DN0
      IF(DNI.GE.1.0D0) GOTO 20
      IRK=IRK-1
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      GOTO 20
C
C   FROM -4713-JAN-1.0 TO 1-JAN-1.0
  15  CONTINUE 
      DNI=DNI+0.5D0
      IRK=INT(DNI/3.6525D2)
      DN0=365*IRK+INT((IRK+2)/4.0D0)
      IF(DNI.GT.1.721118D6) DNI=DNI-1.0D0
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 50
      IRK=IRK-4713
      GOTO 20
  50  IRK=IRK-1
      DH4=(IRK+3)/4.0D0-INT((IRK+3)/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK-4713
      GOTO 20
C
C   TO -4713-JAN-1.0
  16  CONTINUE
      DNI=DNI+1.9307115D6
      IRK=INT(DNI/3.6525D2)
      DN0=365*IRK+INT(IRK/4.0D0)
      DNI=DNI-DN0
      IF(DNI.LT.1.0D0) GOTO 55
      IRK=IRK-9999
      GOTO 20
  55  IRK=IRK-1
      DH4=(IRK+1)/4.0D0-INT((IRK+1)/4.0D0)
      DNI=DNI+3.65D2
      IF(DABS(DH4).LT.1.0D-6) DNI=DNI+1.0D0
      IRK=IRK-9999
C
  20  CONTINUE
      DH4=IRK/4.0D0-INT(IRK/4.0D0)
      DH100=IRK/1.0D2-INT(IRK/1.0D2)
      DH400=IRK/4.0D2-INT(IRK/4.0D2)
      IF(DABS(DH4).LT.1.0D-6) IMESI(3)=29
      IF(IRK.GT.1582.AND.DABS(DH100).LT.1.0D-6) IMESI(3)=28
      IF(IRK.GT.1582.AND.DABS(DH400).LT.1.0D-6) IMESI(3)=29
      IF(IRK.EQ.0) IRK=-1
      IMS=1
  23  CONTINUE
      IF(IMS.NE.13) GOTO 60
      IMS=12
      DNI=DNI-1.0D0
  60  CONTINUE
      IF(IMESI(IMS+1)+1.GT.DNI) GOTO 26
      DNI=DNI-IMESI(IMS+1)
      IMS=IMS+1
      GOTO 23
  26  CONTINUE
      IMESI(3)=28
      RETURN
C
      END
