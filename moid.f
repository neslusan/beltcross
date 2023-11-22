c Calculation of the MOID between the mean orbit of meteor shower and
c  main-belt asteroid
c Version: August 17, 2022

      implicit integer (i,j)
      implicit real*8 (a-h,k-z)
      character*15 name
      character*80 txt
      character*110 ret110
      dimension ea(5),eb(5),q(1500),ecc(1500),jnn(1500)
      dimension som(1500),gom(1500),sk(1500),iau(1500),isol(1500)
      
      open(unit=10,file='allshowers.d',access='sequential')
      open(unit=11,file='object.dat',access='sequential')
      open(unit=15,file='moid.dw',access='sequential')
      read(10,*) ret110

c      rcrit=1.0d-2
      
      j=0
  10  continue
      j=j+1
      read(10,*,end=30) ilp,iau(j),isol(j),ls,al,dl,vg,ash,q(j),ecc(j),s
     *om(j),gom(j),sk(j),jnn(j)
      if(som(j).lt.-0.5d0.or.gom(j).lt.-0.5d0.or.sk(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(dabs(ash+1.0d0).lt.1.0d-5.and.q(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(dabs(ash+1.0d0).lt.1.0d-5.and.ecc(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(q(j).lt.-0.5d0.and.ecc(j).lt.-0.5d0) then
        j=j-1
        goto 10
      endif
      if(q(j).lt.-0.5d0) then
        q(j)=ash*(1.0d0-ecc(j))
        goto 10
      endif
      if(ecc(j).lt.-0.5d0) ecc(j)=1.0d0-q(j)/ash
      goto 10
  30  continue
      jsh=j-1
      close(unit=10)
c       write(*,*) jsh
      
  50  continue
      read(11,*) txt
      read(11,*) name
      read(11,*) txt
      read(11,*) a
      read(11,*) txt
      read(11,*) eb(2)
      read(11,*) txt
      read(11,*) eb(3)
      read(11,*) txt
      read(11,*) eb(4)
      read(11,*) txt
      read(11,*) eb(5)
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
c       write(*,*) rcrit
c       read(*,*) jjyy
      eb(1)=a*(1.0d0-eb(2))
      do i=1,jsh
        ea(1)=q(i)
        ea(2)=ecc(i)
        ea(3)=som(i)
        ea(4)=gom(i)
        ea(5)=sk(i)
        call approach(ea,eb,fkpost,fapost,rpost,fkpre,fapre,rpre)
        if(rpost.le.rcrit.or.rpre.le.rcrit) then
          write(15,300) name,fapost,fapre,rpost,rpre
          write(15,320) iau(i),isol(i),jnn(i),fkpost,fkpre
          write(15,340) q(i),ecc(i),som(i),gom(i),sk(i)
        endif
c        write(*,*) i
      enddo
c      goto 50
  70  continue
      close(unit=11)
      close(unit=15)
      stop

 300  format(a15,2f9.3,2f8.4)
 320  format(i5,i4,i6,2f9.3)
 340  format(f7.4,f8.4,3f9.3)
      end
C ---------------------------------------------------------------------
      subroutine approach(ea,eb,fkpost,fapost,rpost,fkpre,fapre,rpre)
C ---------------------------------------------------------------------
C      - CALCULATION OF minimum orbit intersection distance (MOID)
C        BETWEEN TWO KEPLERIAN ORBITS
C        
C INPUT:
C    ORBITAL ELEMENTS OF THE FIRST ORBIT
C       ea(1) - PERIHELION DISTANCE [AU]
C       ea(2) - ECCENTRICITY
C       ea(3) - ARGUMENT OF PERIHELION [DEGREES]
C       ea(4) - LONGITUDE OF NODE [DEGREES]
C       ea(5) - INCLINATION TO ECLIPTIC [DEGREES]
C    ORBITAL ELEMENTS OF THE SECOND ORBIT
C       eb(1) - PERIHELION DISTANCE [AU]
C       eb(2) - ECCENTRICITY
C       eb(3) - ARGUMENT OF PERIHELION [DEGREES]
C       eb(4) - LONGITUDE OF NODE [DEGREES]
C       eb(5) - INCLINATION TO ECLIPTIC [DEGREES]
C
C OUTPUT:
C   THE OUTPUT QUANTITIES AT THE POST-PERIHELION ARC ARE:
C      fkpost - TRUE ANOMALY OF THE FIRST BODY IN THE MOMENT OF ITS
C               NEAREST POST-PERIHELION APPROACH TO THE ORBIT OF THE
C               SECOND BODY [DEGREES]
C      fapost - TRUE ANOMALY OF THE SECOND BODY IN THE MOMENT OF ITS
C               NEAREST APPROACH TO THE POST-PERIHELION ARC OF THE
C               ORBIT OF THE FIRST BODY [DEGREES]
C      rpost - MOID AT THE POST-PERIHELION ARC [AU]
C   THE OUTPUT QUANTITIES AT THE PRE-PERIHELION ARC ARE ANALOGOUSLY
C   STORED IN "fkpre", "fapre", AND "rpost"
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      implicit real*8 (a-h,k-z)
      implicit integer (i,j)
      dimension ea(5),eb(5)
C
      pi=4.0d0*datan(1.0d0)
      pi180=pi/0.18d3
      gauss=0.17202098950d-1
C
      cgok=dcos(pi180*ea(4))
      sgok=dsin(pi180*ea(4))
      cskk=dcos(pi180*ea(5))
      sskk=dsin(pi180*ea(5))
      prk=ea(1)*(1.0d0+ea(2))
        cgoa=dcos(pi180*eb(4))
        sgoa=dsin(pi180*eb(4))
        cska=dcos(pi180*eb(5))
        sska=dsin(pi180*eb(5))
        pra=eb(1)*(1.0d0+eb(2))
C
C  === 1-ST ARC ===
C
      rpost=1.0d6
      if(ea(2).gt.0.999d0) then
        ivkup=175
      else
        ivkup=179
      endif
      do ivk=0,ivkup
        men=1.0d0+ea(2)*dcos(pi180*ivk)
        rk=1.0d6
        if(dabs(men).lt.1.0d-9) goto 100
        rk=prk/men
 100  continue
        c2q=dcos(pi180*(ivk+ea(3)))
        s2q=dsin(pi180*(ivk+ea(3)))
        xk=rk*(c2q*cgok-s2q*cskk*sgok)
        yk=rk*(c2q*sgok+s2q*cskk*cgok)
        zk=rk*s2q*sskk
C
        do iva=0,359
          men=1.0d0+eb(2)*dcos(pi180*iva)
          ra=pra/men
          c2q=dcos(pi180*(iva+eb(3)))
          s2q=dsin(pi180*(iva+eb(3)))
          xa=ra*(c2q*cgoa-s2q*cska*sgoa)
          ya=ra*(c2q*sgoa+s2q*cska*cgoa)
          za=ra*s2q*sska
          dr=dsqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya)+(zk-za)*(zk-za))
          if(dr.gt.1.0d2) goto 240
          if(dr.ge.rpost) goto 40
          rpost=dr
          ifkpost=ivk
          ifapost=iva
  40  continue
        enddo
      enddo
c more detail:
      do ivk=1,200
        vk=(ifkpost-1)*1.0d0+ivk*1.0d-2
        if(vk.lt.0.0d0) vk=0.0d0
        men=1.0d0+ea(2)*dcos(pi180*vk)
        rk=prk/men
        c2q=dcos(pi180*(vk+ea(3)))
        s2q=dsin(pi180*(vk+ea(3)))
        xk=rk*(c2q*cgok-s2q*cskk*sgok)
        yk=rk*(c2q*sgok+s2q*cskk*cgok)
        zk=rk*s2q*sskk
C
        do iva=1,200
          va=(ifapost-1)*1.0d0+iva*1.0d-2
          men=1.0d0+eb(2)*dcos(pi180*va)
          ra=pra/men
          c2q=dcos(pi180*(va+eb(3)))
          s2q=dsin(pi180*(va+eb(3)))
          xa=ra*(c2q*cgoa-s2q*cska*sgoa)
          ya=ra*(c2q*sgoa+s2q*cska*cgoa)
          za=ra*s2q*sska
          dr=dsqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya)+(zk-za)*(zk-za))
          if(dr.ge.rpost) goto 140
          rpost=dr
          fkpost=vk
          fapost=va
c            xaw1=xa
c            yaw1=ya
c            zaw1=za
c            xkw1=xk
c            ykw1=yk
c            zkw1=zk
 140  continue
        enddo
      enddo

 240  continue
C
C  === 2-ND ARC ===
C
      rpre=1.0d6
      do ivk=359,ivkup,-1
        men=1.0d0+ea(2)*dcos(pi180*ivk)
        rk=1.0d6
        if(dabs(men).lt.1.0d-9) goto 300
        rk=prk/men
 300  continue
        c2q=dcos(pi180*(ivk+ea(3)))
        s2q=dsin(pi180*(ivk+ea(3)))
        xk=rk*(c2q*cgok-s2q*cskk*sgok)
        yk=rk*(c2q*sgok+s2q*cskk*cgok)
        zk=rk*s2q*sskk
C
        do iva=0,359
          men=1.0d0+eb(2)*dcos(pi180*iva)
          ra=pra/men
          c2q=dcos(pi180*(iva+eb(3)))
          s2q=dsin(pi180*(iva+eb(3)))
          xa=ra*(c2q*cgoa-s2q*cska*sgoa)
          ya=ra*(c2q*sgoa+s2q*cska*cgoa)
          za=ra*s2q*sska
          dr=dsqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya)+(zk-za)*(zk-za))
          if(dr.gt.1.0d2) goto 440
          if(dr.ge.rpre) goto 340
          rpre=dr
          ifkpre=ivk
          ifapre=iva
 340  continue
        enddo
      enddo
c more detail:
      do ivk=1,200
        vk=(ifkpre-1)*1.0d0+ivk*1.0d-2
        if(vk.lt.0.0d0) vk=0.0d0
        men=1.0d0+ea(2)*dcos(pi180*vk)
        rk=prk/men
        c2q=dcos(pi180*(vk+ea(3)))
        s2q=dsin(pi180*(vk+ea(3)))
        xk=rk*(c2q*cgok-s2q*cskk*sgok)
        yk=rk*(c2q*sgok+s2q*cskk*cgok)
        zk=rk*s2q*sskk
C
        do iva=1,200
          va=(ifapre-1)*1.0d0+iva*1.0d-2
          men=1.0d0+eb(2)*dcos(pi180*va)
          ra=pra/men
          c2q=dcos(pi180*(va+eb(3)))
          s2q=dsin(pi180*(va+eb(3)))
          xa=ra*(c2q*cgoa-s2q*cska*sgoa)
          ya=ra*(c2q*sgoa+s2q*cska*cgoa)
          za=ra*s2q*sska
          dr=dsqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya)+(zk-za)*(zk-za))
          if(dr.ge.rpre) goto 350
          rpre=dr
          fkpre=vk
          fapre=va
c            xkw2=xk
c            ykw2=yk
c            zkw2=zk
c            xaw2=xa
c            yaw2=ya
c            zaw2=za
 350  continue
        enddo
      enddo

 440  continue
c        write(*,*) xkw1,ykw1,zkw1
c        write(*,*) xaw1,yaw1,zaw1
c        fz1=datan(yaw1/xaw1)/pi180+1.80d2
c        write(*,*) fz1
c        write(*,*) xkw2,ykw2,zkw2
c        write(*,*) xaw2,yaw2,zaw2
c        fz2=datan(yaw2/xaw2)/pi180+3.60d2
c        write(*,*) fz2
      return

      end
