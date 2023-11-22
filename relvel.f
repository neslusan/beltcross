c Calculation of the relative velocity of asteroid to a mean meteoroid
c   in the positions of the closest approach; calculation of the angle
c   between the heliocentric velocity vectors
c Version: September 7, 2023

      implicit integer (i,j)
      implicit real*8 (a-h,k-z)
      character*15 name
      character*80 txt

      pi=4.0d0*datan(1.0d0)
      pi180=pi/1.80d2
      gauss=1.7202098950d-2
      aukm=1.49597870d8
      day=2.4d1*3.6d3
      prem=aukm/day
            
      open(unit=10,file='object.dat',access='sequential')
      open(unit=11,file='moid.dw',access='sequential')
      open(unit=15,file='relativeV.dw',access='sequential')
      
        read(10,*) txt
        read(10,*) name
        read(10,*) txt
        read(10,*) a
        read(10,*) txt
        read(10,*) ecc
        read(10,*) txt
        read(10,*) som
        read(10,*) txt
        read(10,*) gom
        read(10,*) txt
        read(10,*) sk
        read(10,*) txt
        read(10,*) irk
        read(10,*) txt
        read(10,*) imes
        read(10,*) txt
        read(10,*) day
        read(10,*) txt
        read(10,*) mm
        read(10,*) txt
        read(10,*) jyear
        read(10,*) txt
        read(10,*) rcrit
      close(unit=10)
      
  10  continue
      read(11,300,end=50) name,fapost,fapre,rpost,rpre
      read(11,*) iau,isol,jnn,fkpost,fkpre
      read(11,*) qsh,esh,sosh,gosh,sksh
      jarc1=0
      if(rpost.le.rcrit) then
        jarc1=1
        men2=1.0d0+ecc*ecc+2.0d0*ecc*dcos(fapost*pi180)
        cvfi=-ecc*dsin(fapost*pi180)/dsqrt(men2)
        svfi=(1.0d0+ecc*dcos(fapost*pi180))/dsqrt(men2)
        if(cvfi.gt.1.0d-9) then
          vfi=datan(svfi/cvfi)
          vfi=dabs(vfi)
          if(cvfi.lt.0.0d0.and.svfi.ge.0.0d0) vfi=pi-vfi
          if(cvfi.lt.0.0d0.and.svfi.lt.0.0d0) vfi=pi+vfi
          if(cvfi.ge.0.0d0.and.svfi.lt.0.0d0) vfi=2.0d0*pi-vfi
        else
          vfi=pi/2.0d0
        endif
        rast1=a*(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(fapost*pi180))
        vast1=gauss*dsqrt(2.0d0/rast1-1.0d0/a)
        c3q=dcos((som+gom)*pi180-vfi)
        s3q=dsin((som+gom)*pi180-vfi)
        cgo=dcos(gom*pi180)
        sgo=dsin(gom*po180)
        csk=dcos(sk*pi180)
        ssk=dsin(sk*pi180)
        vax1=-vast1*(c3q*cgo-s3q*csk*sgo)
        vay1=-vast1*(c3q*sgo+s3q*csk*cgo)
        vaz1=-vast1*s3q*ssk
        
        men2sh=1.0d0+esh*esh+2.0d0*esh*dcos(fkpost*pi180)
        cvfish=-esh*dsin(fkpost*pi180)/dsqrt(men2sh)
        svfish=(1.0d0+esh*dcos(fkpost*pi180))/dsqrt(men2sh)
        if(cvfish.gt.1.0d-9) then
          vfish=datan(svfish/cvfish)
          vfish=dabs(vfish)
          if(cvfish.lt.0.0d0.and.svfish.ge.0.0d0) vfish=pi-vfish
          if(cvfish.lt.0.0d0.and.svfish.lt.0.0d0) vfish=pi+vfish
          if(cvfish.ge.0.0d0.and.svfish.lt.0.0d0) vfish=2.0d0*pi-vfish
        else
          vfish=pi/2.0d0
        endif
        rsh1=qsh*(1.0d0+esh)/(1.0d0+esh*dcos(fkpost*pi180))
        vsh1=gauss*dsqrt(2.0d0/rsh1-(1.0d0-esh)/qsh)
        c3qsh=dcos((sosh+gosh)*pi180-vfish)
        s3qsh=dsin((sosh+gosh)*pi180-vfish)
        cgosh=dcos(gosh*pi180)
        sgosh=dsin(gosh*po180)
        csksh=dcos(sksh*pi180)
        ssksh=dsin(sksh*pi180)
        vshx1=-vsh1*(c3qsh*cgosh-s3qsh*csksh*sgosh)
        vshy1=-vsh1*(c3qsh*sgosh+s3qsh*csksh*cgosh)
        vshz1=-vsh1*s3qsh*ssksh
        vrelx1=vshx1-vax1
        vrely1=vshy1-vay1
        vrelz1=vshz1-vaz1
        vrel1=dsqrt(vrelx1*vrelx1+vrely1*vrely1+vrelz1*vrelz1)*prem
        cangl1=(vax1*vshx1+vay1*vshy1+vaz1*vshz1)/vast1/vsh1
        if(cangl1.gt.1.0d0) then
          write(*,*) 'Unknown error; |cos(ANGLE)| > 1.'
          write(*,*) 'Program terminated.'
          stop
        endif
        sangl1=dsqrt(1.0d0-cangl1*cangl1)
        if(dabs(cangl1).gt.1.0d-9) then
          angl1=datan(sangl1/cangl1)
          if(angl1.lt.0.0d0) angl1=angl1+pi
          angl1=angl1/pi180
        else
          angl1=9.0d1
        endif
        write(15,350) name,fapost,rpost,rast1,vrel1,angl1
        write(15,360) iau,isol,jnn,1
      endif

      jarc2=0
      if(rpre.le.rcrit) then
        jarc2=1
        men2=1.0d0+ecc*ecc+2.0d0*ecc*dcos(fapre*pi180)
        cvfi=-ecc*dsin(fapre*pi180)/dsqrt(men2)
        svfi=(1.0d0+ecc*dcos(fapre*pi180))/dsqrt(men2)
        if(cvfi.gt.1.0d-9) then
          vfi=datan(svfi/cvfi)
          vfi=dabs(vfi)
          if(cvfi.lt.0.0d0.and.svfi.ge.0.0d0) vfi=pi-vfi
          if(cvfi.lt.0.0d0.and.svfi.lt.0.0d0) vfi=pi+vfi
          if(cvfi.ge.0.0d0.and.svfi.lt.0.0d0) vfi=2.0d0*pi-vfi
        else
          vfi=pi/2.0d0
        endif
        rast2=a*(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(fapre*pi180))
        vast2=gauss*dsqrt(2.0d0/rast2-1.0d0/a)
        c3q=dcos((som+gom)*pi180-vfi)
        s3q=dsin((som+gom)*pi180-vfi)
        cgo=dcos(gom*pi180)
        sgo=dsin(gom*po180)
        csk=dcos(sk*pi180)
        ssk=dsin(sk*pi180)
        vax2=-vast2*(c3q*cgo-s3q*csk*sgo)
        vay2=-vast2*(c3q*sgo+s3q*csk*cgo)
        vaz2=-vast2*s3q*ssk
        men2sh=1.0d0+esh*esh+2.0d0*esh*dcos(fkpre*pi180)
        cvfish=-esh*dsin(fkpre*pi180)/dsqrt(men2sh)
        svfish=(1.0d0+esh*dcos(fkpre*pi180))/dsqrt(men2sh)
        if(cvfish.gt.1.0d-9) then
          vfish=datan(svfish/cvfish)
          vfish=dabs(vfish)
          if(cvfish.lt.0.0d0.and.svfish.ge.0.0d0) vfish=pi-vfish
          if(cvfish.lt.0.0d0.and.svfish.lt.0.0d0) vfish=pi+vfish
          if(cvfish.ge.0.0d0.and.svfish.lt.0.0d0) vfish=2.0d0*pi-vfish
        else
          vfish=pi/2.0d0
        endif
        rsh2=qsh*(1.0d0+esh)/(1.0d0+esh*dcos(fkpre*pi180))
        vsh2=gauss*dsqrt(2.0d0/rsh2-(1.0d0-esh)/qsh)
        c3qsh=dcos((sosh+gosh)*pi180-vfish)
        s3qsh=dsin((sosh+gosh)*pi180-vfish)
        cgosh=dcos(gosh*pi180)
        sgosh=dsin(gosh*po180)
        csksh=dcos(sksh*pi180)
        ssksh=dsin(sksh*pi180)
        vshx2=-vsh2*(c3qsh*cgosh-s3qsh*csksh*sgosh)
        vshy2=-vsh2*(c3qsh*sgosh+s3qsh*csksh*cgosh)
        vshz2=-vsh2*s3qsh*ssksh
        vrelx2=vshx2-vax2
        vrely2=vshy2-vay2
        vrelz2=vshz2-vaz2
        vrel2=dsqrt(vrelx2*vrelx2+vrely2*vrely2+vrelz2*vrelz2)*prem
        cangl2=(vax2*vshx2+vay2*vshy2+vaz2*vshz2)/vast2/vsh2
        if(cangl2.gt.1.0d0) then
          write(*,*) 'Unknown error; |cos(ANGLE)| > 1.'
          write(*,*) 'Program terminated.'
          stop
        endif
        sangl2=dsqrt(1.0d0-cangl2*cangl2)
        if(dabs(cangl2).gt.1.0d-9) then
          angl2=datan(sangl2/cangl2)
          if(angl2.lt.0.0d0) angl2=angl2+pi
          angl2=angl2/pi180
        else
          angl2=9.0d1
        endif
        write(15,350) name,fapre,rpre,rast2,vrel2,angl2
        write(15,360) iau,isol,jnn,-1
      endif
      goto 10
  50  continue
      close(unit=11)
      close(unit=15)
      stop
      
 300  format(a15,2f9.3,2f8.4)
 350  format(a15,f9.3,f8.4,2f7.2,f8.2)
 360  format(i5,i4,i6,i4)
      end
