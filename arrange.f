c Arranging the result in the order of increasind date
c Version: September 8, 2023

      implicit integer (i,j)
      implicit real*8 (a-h,k-z)
      character*1 cf(10),cc
      character*4 cislo
      character*11 arc,arcx
      character*13 ret13
      character*15 nam0,name
      character*16 filin,filout
      character*80 txt
      dimension iau(1000),isol(1000),jnn(1000),arc(1000),moid(1000)
      dimension ra(1000),vr(1000),angl(1000)
      dimension irk(1000),im(1000),den(1000)
      
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

      open(unit=10,file='object.dat',access='sequential')
      
        read(10,*) txt
        read(10,*) nam0
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
        read(10,*) iyr
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
      cislo=cf(it4)//cf(it3)//cf(it2)//cf(it1)
      filin='unarranged'//cislo//'.d'
      open(unit=11,file=filin,access='sequential')
        read(11,300) ret13,name
        if(nam0.ne.name) then
      write(*,*) 'Unknown reason of error; different names of the object
     *'
          write(*,*) 'in files <object.dat> and <unarranged',cislo,'.d.'
          write(*,*) 'Program terminated.'
          stop
        endif
        read(11,*) txt
        j=0
  30  continue
        j=j+1
        if(j.gt.1000) then
          write(*,*) 'Capacity of arrays exceeded.'
          write(*,*) 'The number of approaches is extremely large.'
          write(*,*) 'Program terminated.'
          write(*,*) ' '
      write(*,*) 'Try to run the program for a smaller value of critical
     *'
      write(*,*) 'MOID. Or, do not run code <arrange.exes>. Then, see'
          write(*,*) 'the result in file <unarranged????.d>, where'
          write(*,*) 'the approaches are not arranged by date, however.'
          stop
        endif
      read(11,320,end=50) cc,iau(j),isol(j),jnn(j),arc(j),moid(j),ra(j),
     *vr(j),angl(j),irk(j),im(j),den(j)
        goto 30
  50  continue
        jall=j-1
      close(unit=11)

      filout='datelist'//cislo//'.dat'
      open(unit=15,file=filout,access='sequential')
      
      jrep=0
  60  continue
      jrep=jrep+1
      if(jrep.gt.10000) then
        write(*,*) 'Unknown error; the arranging does not converge.'
        write(*,*) 'Program terminated.'
        stop
      endif
      jyes=0
      do i=1,jall-1
        do j=i+1,jall
          if(irk(i).gt.irk(j)) goto 70
          if(irk(i).eq.irk(j).and.im(i).gt.im(j)) goto 70
      if(irk(i).eq.irk(j).and.im(i).eq.im(j).and.den(i).gt.den(j)) goto 
     *70
          goto 90
  70  continue
          iaux=iau(i)
          isolx=isol(i)
          jnnx=jnn(i)
          arcx=arc(i)
          moidx=moid(i)
          rax=ra(i)
          vrx=vr(i)
          anglx=angl(i)
          irkx=irk(i)
          imx=im(i)
          denx=den(i)
        iau(i)=iau(j)
        isol(i)=isol(j)
        jnn(i)=jnn(j)
        arc(i)=arc(j)
        moid(i)=moid(j)
        ra(i)=ra(j)
        vr(i)=vr(j)
        angl(i)=angl(j)
        irk(i)=irk(j)
        im(i)=im(j)
        den(i)=den(j)
          iau(j)=iaux
          isol(j)=isolx
          jnn(j)=jnnx
          arc(j)=arcx
          moid(j)=moidx
          ra(j)=rax
          vr(j)=vrx
          angl(j)=anglx
          irk(j)=irkx
          im(j)=imx
          den(j)=denx
          jyes=1
  90  continue
        enddo
      enddo
      if(jyes.eq.1) goto 60

      write(15,330) name
      write(15,340)
      do ii=1,jall
      if(iau(ii).le.9) write(15,350) iau(ii),isol(ii),jnn(ii),arc(ii),mo
     *id(ii),ra(ii),vr(ii),angl(ii),irk(ii),im(ii),den(ii)
      if(iau(ii).ge.10.and.iau(ii).le.99) write(15,360) iau(ii),isol(ii)
     *,jnn(ii),arc(ii),moid(ii),ra(ii),vr(ii),angl(ii),irk(ii),im(ii),de
     *n(ii)
      if(iau(ii).ge.100.and.iau(ii).le.999) write(15,370) iau(ii),isol(i
     *i),jnn(ii),arc(ii),moid(ii),ra(ii),vr(ii),angl(ii),irk(ii),im(ii),
     *den(ii)
      if(iau(ii).ge.1000) write(15,380) iau(ii),isol(ii),jnn(ii),arc(ii)
     *,moid(ii),ra(ii),vr(ii),angl(ii),irk(ii),im(ii),den(ii)
      enddo
      write(15,*) ' '
      write(15,*) 'Explanation:'
      write(15,*) 'IAUNo. - IAU number of meteor shower'
      write(15,*) 'Sol. - its solution'
      write(15,*) 'n - number of meteors in the solution'
      write(15,*) '    (if n=-1, then the number is unknown)'
      write(15,*) 'arc - arc (postperihelion or preperihelion) of'
      write(15,*) '      meteor-stream orbit crossed by the object'
      write(15,*) 'MOID - MOID between the mean orbit of stream'
      write(15,*) '       and nominal orbit of object [au]'
      write(15,*) 'r_obj - heliocentric distance of object in the'
      write(15,*) '        position of the closest approach [au]'
      write(15,*) 'v_rel - relative velocity of object and meteoroid,'
      write(15,*) '        which moves with the average velocity, in'
      write(15,*) '        the positions of the closest approach [km/s]'
      write(15,*) 'angle - angle between the heliocentric velocity'
      write(15,*) '        vectors of asteroid and meteoroid moving with
     *'
      write(15,*) '        the average velocity, at the closest approach
     *'
      write(15,*) '        [deg]'
      write(15,*) 'date - date of the closest approach of object to'
      write(15,*) '       the mean orbit of stream [year month day]'
      close(unit=15)
      stop

 300  format(a13,a15)
 320  format(a1,i4,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 330  format('OBJECT NAME: ',a15,/) 
 340  format('IAUNo. Sol.  n      arc        MOID   r_obj  v_rel   angle
     *        date')
 350  format('#000',i1,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 360  format('#00',i2,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 370  format('#0',i3,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
 380  format('#',i4,i4,i6,1x,a11,1x,f8.4,2f7.2,f8.2,2x,i6,i3,f6.2)
      end
