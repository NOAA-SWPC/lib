      subroutine weidelt_dst(model,ndst,dstvec,est,ist)
      implicit real*8 (a-h,o-z)
      parameter (mt=46,nd=416352,nc=8766) !413424, 8766 h/a * 5a
c mt are stuetz points for spline, nd is only for test
      complex*16 qq,qr,expo
      dimension tt(mt),qt(mt),c(3,mt)
      dimension fc(0:nc),tau(5),ztau(5)

      dimension dst(0:nd),fext(0:nd),finn(0:nd)
      dimension thext(0:nd),thinn(0:nd)

!      dimension dst(0:ndst-1),fext(0:ndst-1),finn(0:ndst-1)
!      dimension thext(0:ndst-1),thinn(0:ndst-1)

      REAL*8 ist(0:ndst)
      common /conduct/r(21),sig(20),nl
      common /param/n
      external q

!      open (6,file='weidelt_dst.res')
      pi=4*datan(1.d0)
     
      if (model.eq.1) then      
c     4-layer model of N. Olsen, p. 307 of GJI, v. 133 (1998)
         nl=4
      
         r(1)=6.371d6
         r(2)=r(1)-5.20d5
         r(3)=r(1)-6.46d5 
         r(4)=r(1)-7.91d5  
         r(5)=r(1)-2.9d6       

         sig(1)=1./76.
         sig(2)=1./4.5
         sig(3)=1./5.9
         sig(4)=1./0.662

      elseif (model.eq.4) then
c        Shankland lab2
         nl = 8

         r(1)=6.371d6
         r(2)=r(1)-1.0d5
         r(3)=r(1)-2.0d5
         r(4)=r(1)-3.0d5
         r(5)=r(1)-4.1d5
         r(6)=r(1)-5.2d5
         r(7)=r(1)-6.6d5
         r(8)=r(1)-8.0d5
         r(9)=r(1)-2.9d6

         sig(1) = 10.0**(-2.15)
         sig(2) = 10.0**(-1.8)
         sig(3) = 10.0**(-2.05)
         sig(4) = 10.0**(-2.2)
         sig(5) = 10.0**(-1.05)
         sig(6) = 10.0**(-0.9)
         sig(7) = 10.0**(-0.6)
         sig(8) = 10.0**(0.3)

      elseif (model.eq.5) then
c        Shankland lab3
         nl = 9

         r(1)=6.371d6
         r(2)=r(1)-1.0d5
         r(3)=r(1)-2.0d5
         r(4)=r(1)-3.0d5
         r(5)=r(1)-4.1d5
         r(6)=r(1)-5.2d5
         r(7)=r(1)-6.6d5
         r(8)=r(1)-8.0d5
         r(9)=r(1)-1.5d6
         r(10)=r(1)-2.9d6

         sig(1) = 10.0**(-2.15)
         sig(2) = 10.0**(-1.8)
         sig(3) = 10.0**(-2.05)
         sig(4) = 10.0**(-2.2)
         sig(5) = 10.0**(-0.79)
         sig(6) = 10.0**(-0.56)
         sig(7) = 10.0**(-0.25)
         sig(8) = 10.0**(0.3)
         sig(9) = 10.0**(0.8)

      elseif (model.eq.7) then
c        Utada model b
         nl = 15

         r(1)=6.371d6
         r(2)=r(1)-1.0d5
         r(3)=r(1)-1.5d5
         r(4)=r(1)-2.0d5
         r(5)=r(1)-2.5d5
         r(6)=r(1)-3.0d5
         r(7)=r(1)-3.5d5
         r(8)=r(1)-4.0d5
         r(9)=r(1)-4.5d5
         r(10)=r(1)-5.0d5
         r(11)=r(1)-5.5d5
         r(12)=r(1)-6.0d5
         r(13)=r(1)-6.5d5
         r(14)=r(1)-8.5d5
         r(15)=r(1)-1.5d6
         r(16)=r(1)-2.9d6

         sig(1) = 10.0**(-3.15)
         sig(2) = 10.0**(-2.8)
         sig(3) = 10.0**(-2.7)
         sig(4) = 10.0**(-2.6)
         sig(5) = 10.0**(-2.5)
         sig(6) = 10.0**(-2.4)
         sig(7) = 10.0**(-2.3)
         sig(8) = 10.0**(-0.95)
         sig(9) = 10.0**(-0.75)
         sig(10) = 10.0**(-0.5)
         sig(11) = 10.0**(-0.4)
         sig(12) = 10.0**(-0.35)
         sig(13) = 10.0**(0)
         sig(14) = 10.0**(0.2)
         sig(15) = 10.0**(0.2)

      elseif (model.eq.8) then
c        Utada model c
         nl = 15

         r(1)=6.371d6
         r(2)=r(1)-1.0d5
         r(3)=r(1)-1.5d5
         r(4)=r(1)-2.0d5
         r(5)=r(1)-2.5d5
         r(6)=r(1)-3.0d5
         r(7)=r(1)-3.5d5
         r(8)=r(1)-4.1d5
         r(9)=r(1)-4.5d5
         r(10)=r(1)-5.2d5
         r(11)=r(1)-5.5d5
         r(12)=r(1)-6.0d5
         r(13)=r(1)-6.5d5
         r(14)=r(1)-8.5d5
         r(15)=r(1)-1.5d6
         r(16)=r(1)-2.9d6

         sig(1) = 10.0**(-3.1)
         sig(2) = 10.0**(-2.7)
         sig(3) = 10.0**(-2.7)
         sig(4) = 10.0**(-2.55)
         sig(5) = 10.0**(-2.5)
         sig(6) = 10.0**(-2.4)
         sig(7) = 10.0**(-2.3)
         sig(8) = 10.0**(-1.35)
         sig(9) = 10.0**(-1.3)
         sig(10) = 10.0**(-0.05)
         sig(11) = 10.0**(-0.1)
         sig(12) = 10.0**(-0.1)
         sig(13) = 10.0**(-0.1)
         sig(14) = 10.0**(0.2)
         sig(15) = 10.0**(0.2)

      elseif (model.eq.2) then
c     D+model of N. Olsen, p. 307 of GJI, v. 133 (1998)
         nl=10
         ntau=5
         tau(1)=  1.65d3
         tau(2)=  7.42d3
         tau(3)=  39.8d3
         tau(4)= 318.4d3
         tau(5)=1040.0d3
         ztau(1)=0.d0
         ztau(2)=3.51d5
         ztau(3)=6.09d5
         ztau(4)=8.53d5
         ztau(5)=1.322d6

c     correct tau and ztau for spherical geometry:
c     -------------------------------------------
         a=6.371d6
         n=1
         expo=1.d0/(2*n+1)
         do it=1,ntau
            z=ztau(it)         
            rr=a*((1-(n+1)*(z/a))/(1+n*(z/a)))**expo
            ztau(it)=a-rr
            f=((n+1)*(a/rr)**n+n*(rr/a)**(n+1))/(2*n+1)
            tau(it)=tau(it)/f**2         
         enddo

        
         dd=2.d4
         do il=1,nl
            sig(il)=1.d-3
         enddo
           
         sig(1)=tau(1)/dd
         sig(3)=tau(2)/dd
         sig(5)=tau(3)/dd
         sig(7)=tau(4)/dd
         sig(9)=tau(5)/dd
       
         r(1)=6.371d6
         r(2)=r(1)-dd
         r(3)=r(1)-ztau(2)+dd/2
         r(4)=r(3)-dd
         r(5)=r(1)-ztau(3)+dd/2
         r(6)=r(5)-dd
         r(7)=r(1)-ztau(4)+dd/2
         r(8)=r(7)-dd
         r(9)=r(1)-ztau(5)+dd/2
         r(10)=r(9)-dd
         r(11)=r(1)-2.9d6
      elseif (model.eq.3) then      
c     7-layer model of Constable and Constable (2000 in G^3)

         nl=7
      
         r(1)=6.371d6
         r(2)=r(1) - 4000.0
         r(3)=r(1) - 220000.0
         r(4)=r(1) - 440000.0
         r(5)=r(1) - 670000.0 
         r(6)=r(1) - 960000.0 
         r(7)=r(1) -1250000.0
         r(8)=r(1) -2800000.0

         sig(1)=1./1.2997
         sig(2)=1./44.4947
         sig(3)=1./148.4817
         sig(4)=1./75.7809
         sig(5)=1./4.9898
         sig(6)=1./0.1883
         sig(7)=1./0.01

      endif    
 
!      write (6,*)
!      write (6,*) 'conductivity structure'
!      write (6,*) '----------------------'
!      do il=1,nl
!         write (6,100) r(il),r(il+1),sig(il)
!      enddo

!      write (6,110) r(nl+1)

      n=1
      dt=3600.d0
!      write (6,120) n
!      write (6,130) dt,nc
      ta=dt/2

c     calculation of step-on response
c     ------------------------------- 
c     at logarithmically distributed time lags tt
c     -------------------------------------------
c     with 10 or 15 data points per decade
c     ------------------------------------  
c     first time lag is ta=dt/2
c     -------------------------          
      call frt10(ta,mt,tt,qt,q)           
      
c     filter coefficients:
c     -------------------
!      write (6,*)
!      write (6,*) 'filter coefficients:'
!      write (6,*) '------------------- ' 
!      write (6,135)

      call spline1(mt,qt,c)
      
      tl1=dlog(tt(1))
      tl2=dlog(tt(mt))
      do ic=0,nc
         t=ic*dt
         if (ic.eq.0) then
            fc(0)=qt(1)
            fp=fc(0)
         else           
            tlint=dlog(t+dt/2)
            fm=fp
            call spline2(mt,tlint,tl1,tl2,c,fp)
            fc(ic)=fp-fm           
         endif
!         write (6,140) ic,t,fc(ic)
      enddo 

!c     testing filter procedure with perioc signals:
!c     ----------------------------------------------
!      write (6,*)
!
!      write (6,*) 'test case with periodic signal '
!      write (6,*) '------------------------------ '
!c     nd = number of data (time step dt)       
!      period=24*dt
!      ampl=1.d0
!
!c     theoretical q-ratio:
!      qr=qq(n,period)
!      write (6,150) period,qr
!
!c     generate test data dst(t)=ampl*sin(2*pi*t/period)
!      do j=0,nd
!         arg=2*pi*j*dt/period
!         expo=cdexp((0.d0,1.d0)*arg)
!         dst(j)=ampl*dimag(expo)
!         thext(j)=ampl*dimag(1/(1+qr)*expo)
!         thinn(j)=dst(j)-thext(j)          
!      enddo   
!      
!      call filt(nc,fc,nd,dst,fext,finn)
!
!      write (6,*)   
!      cext=1/(1+dble(qr))
!      cinn=dble(qr)*cext
!      write (6,*) 'cext=1/[1+Re(q)]=', cext,', ext_app=cext*dst'
!      write (6,*) 'cinn=Re(q)/[1+Re(q)]=',cinn,
!     &', int_app=cinn*dst'
!      write (6,*)
! 
!      write (6,160)  
!      do j=0,nd
!         write (6,170) j,dst(j),thext(j),fext(j),cext*dst(j),
!     &                          thinn(j),finn(j),cinn*dst(j) 
!      enddo 
               
c     *********************************************************** 

 100  format (' from r = ',1pe9.3,' m to r = ', 1pe9.3,
     & ' m: sigma = ', 1pe10.3,' S/m')
 110  format (' from r = ',1pe9.3,
     &' m to r = 0 m: sigma = infinity',/)
 120  format (' spherical degree n = ', i2)
 130  format (' time step dt = ', f8.1,' s, nc = ',i3,
     & ' coefficients')
 135  format ('  #     lag [s]   coefficient')
 140  format (i5,1pe12.4,2x,1pe12.4)  
 150  format (' period =',1pe11.4,0p,'s, q = (',f6.4,',',f6.4,')')
 160  format ('  #    d_st   ext_th ext_flt ext_app',
     &        '   int_th int_flt int_app') 
 170  format (i3,f8.4,2(1x,3f8.4))


c     nd = number of data (time step dt = 1h)       

      call filt(nc,fc,ndst,dstvec,est,ist)
      return
      end 
      

c     ***********************************************************

      subroutine weidelt_qq(period,a,b)
      complex*16 qq,cval
      real*8 period,a,b
      n = 1                     !SH degree
      cval = qq(n,period)
      a = dble(cval)
      b = dimag(cval)
      return
      end


      complex*16 function qq(n,period)

c     function program for the computation of the ratio of
c     internal to external part for a spherically layered 
c     conductor
c     input parameter:
c     n = spherical degree
c     period = period in seconds
c     the conductivity structure is transferred via the 
c     common block 'conduct' defined in the calling program:

c     nl = number of uniform layers (excluding a perfectly 
c     conducting  core)
c     r(1)    > r > r(2) :   sigma = sig(1)
c     r(2)    > r > r(3) :   sigma = sig(2)
           
c                    until

c     r(nl)   > r > r(nl+1): sigma = sig(nl)
c     r(nl+1) > r > 0 :      sigma = infinity
c     r in meter, sig in Siemens/meter

c     the program should work also for very small periods, 
c     where it models the response of a layered plane conductor  
        

      implicit real*8 (a-h,o-z)
 
      common /conduct/r(21),sig(20),nl

      complex*16 p(2),pd(2),q(2),qd(2),dp,dq
      complex*16 b,d,e,i,k,rm,rmd,rp,rpd,z(2),zz
      complex*16 v1,v2,v3,v4,v5,v6      
      data pi,eps,i,zlimit/3.141592654d0,1.d-7,(0.d0,1.d0),3.0d0/

      fac1=1.d0
      do j=1,n
         fac1=fac1*(2*j+1)
      enddo
      fac2=(-1)**n*fac1/(2*n+1)

      do il=nl,1,-1
         k=cdsqrt(8.d-7*i*pi*pi*sig(il)/period)
         z(1)=k*r(il)
         z(2)=k*r(il+1)
        
c        calculate spherical bessel functions with small argument
c        by power series (abramowitz & Stegun 10.2.5, 10.2.6 
c        and 10.2.4):

         if (cdabs(z(1)).lt.zlimit) then                  
            do m=1,2
               p(m)=1.d0 
               q(m)=1.d0
               pd(m)=n
               qd(m)=-(n+1)
               zz=z(m)**2/2
               j=0
               dp=1.d0
               dq=1.d0 
          
 10            j=j+1             
               dp=dp*zz/dfloat(j)/dfloat(2*j+1+2*n)
               dq=dq*zz/dfloat(j)/dfloat(2*j-1-2*n)
               p(m)=p(m)+dp
               q(m)=q(m)+dq
               pd(m)=pd(m)+dp*(2*j+n)
               qd(m)=qd(m)+dq*(2*j-n-1) 
               if (cdabs(dp).gt.eps.or.cdabs(dq).gt.eps) goto 10
            
               p(m)=p(m)*z(m)**n/fac1
               q(m)=q(m)*z(m)**(-n-1)*fac2
               q(m)=(-1)**(n+1)*pi/2*(p(m)-q(m))
               pd(m)=pd(m)*z(m)**(n-1)/fac1
               qd(m)=qd(m)*z(m)**(-n-2)*fac2
               qd(m)=(-1)**(n+1)*pi/2*(pd(m)-qd(m))
            enddo

            v1= p(2)/p(1)
            v2=pd(1)/p(1)
            v3=pd(2)/p(1)    
            v4= q(1)/q(2)
            v5=qd(1)/q(2)
            v6=qd(2)/q(2)
         endif

c        calculate spherical bessel functions with large argument 
c        the exponential behaviour is split off and is treated 
c        separately (abramowitz & stegun 10.2.9 and 10.2.15)
    
         if (cdabs(z(1)).ge.zlimit) then
            do m=1,2
               zz=2*z(m)
               rm=1.d0
               rp=1.d0
               rmd=1.d0
               rpd=1.d0 
               d=1.d0
               sg=1.d0 
               do j=1,n
                  d=d*dfloat((n+1-j)*(n+j))/dfloat(j)/zz
                  sg=-sg
                  rp=rp+d              
                  rm=rm+sg*d
                  rmd=rmd+sg*d*(j+1)
                  rpd=rpd+d*(j+1)
               enddo                
               e=cdexp(-2*z(m))
               p(m)=(rm-sg*rp*e)/zz                
               q(m)=(pi/zz)*rp
               pd(m)=(rm+sg*rp*e)/zz-2*(rmd-sg*rpd*e)/zz**2
               qd(m)=-q(m)-2*pi*rpd/zz**2
            enddo
           
            e=cdexp(-(z(1)-z(2)))
            v1= p(2)/p(1)*e
            v2=pd(1)/p(1)
            v3=pd(2)/p(1)*e    
            v4= q(1)/q(2)*e
            v5=qd(1)/q(2)*e
            v6=qd(2)/q(2)
         endif 

         if (il.eq.nl) then
            b=k*(v2-v5*v1)/(1-v4*v1)
         else
            b=k*((v2-v5*v1)*b+k*(v5*v3-v2*v6))/
     &          ((1- v4*v1)*b+k*(v4*v3-   v6))
         endif 
      enddo
      qq=n/dfloat(n+1)
      qq=qq*(1-(2*n+1)/(r(1)*b+n+1))
      return
      end        

c     ***************************************************************
   
      subroutine frt15(ta,mt,t,ft,fun)
c     subroutine zur berechnung der fouriertransformierten der 
c     funktion fun(omega), definiert im function-unterprogramm 
c     fun(omega) die transformierte funktion wird an logarith-
c     misch aequidistanten stuetzstellen berechnet mit 10 
c     stuetzstellen pro dekade.

c     bedeutung der parameter:
c     ta = erste stuetzstelle im t-bereich (input)
c     mt = anzahl der stuetzstellen (input)
c     t = bereich der logarithmisch aequidistanten stuetzst.
c     (output)
c     ft = bereich der transformierten funktion (output)
c     fun = name des function-unterprogramms  zur berechnung der 
c     frequenzabhaengigen funktion (=external im rufenden 
c     programm)

      implicit real*8 (a-h,o-z)
      dimension h(261),t(*),ft(*)
      complex*16 fun,fomega
      data pi,q/3.141592653589793d0,1.16591440118d0/
      data nc0,nc,(h(i),i=1,261)/131,261,     
     * 1.22479786541788E-14, 1.54192967639227E-14, 1.94117374688725E-14,
     * 2.44379391826438E-14, 3.07655296045224E-14, 3.87315247584115E-14,
     * 4.87600766518288E-14, 6.13853323440408E-14, 7.72795102586419E-14,
     * 9.72891998222996E-14, 1.22479763644888E-13, 1.54192998764276E-13,
     * 1.94117332378691E-13, 2.44379449340851E-13, 3.07655217862637E-13,
     * 3.87315353862113E-13, 4.87600622048614E-13, 6.13853519826186E-13,
     * 7.72794835628161E-13, 9.72892361114413E-13, 1.22479714315010E-12,
     * 1.54193065821159E-12, 1.94117241224488E-12, 2.44379573251901E-12,
     * 3.07655049423362E-12, 3.87315582831118E-12, 4.87600310798131E-12,
     * 6.13853942926516E-12, 7.72794260484020E-12, 9.72893142940276E-12,
     * 1.22479608037012E-11, 1.54193210290835E-11, 1.94117044838711E-11,
     * 2.44379840210162E-11, 3.07654686531942E-11, 3.87316076129889E-11,
     * 4.87599640229296E-11, 6.13854854468548E-11, 7.72793021373535E-11,
     * 9.72894827333032E-11, 1.22479379068006E-10, 1.54193521541317E-10,
     * 1.94116621738380E-10, 2.44380415354301E-10, 3.07653904706080E-10,
     * 3.87317138909871E-10, 4.87598195532535E-10, 6.13856818326323E-10,
     * 7.72790351790917E-10, 9.72898456247223E-10, 1.22478885769237E-09,
     * 1.54194192110152E-09, 1.94115710196348E-09, 2.44381654464787E-09,
     * 3.07652220313322E-09, 3.87319428599931E-09, 4.87595083027718E-09,
     * 6.13861049329638E-09, 7.72784600349519E-09, 9.72906274505842E-09,
     * 1.22477822989255E-08, 1.54195636806913E-08, 1.94113746338574E-08,
     * 2.44384324047403E-08, 3.07648591399132E-08, 3.87760267728197E-08,
     * 4.88128263741473E-08, 6.13614586713336E-08, 7.72129587894784E-08,
     * 9.72957975465685E-08, 1.22542558523773E-07, 1.54218844496570E-07,
     * 1.94047949414646E-07, 2.44347729158895E-07, 3.07688910323352E-07,
     * 3.87395426942958E-07, 4.87545930506262E-07, 6.13817876642722E-07,
     * 7.72748769037344E-07, 9.73034376362735E-07, 1.22472917373381E-06,
     * 1.54198565562436E-06, 1.94095560343693E-06, 2.44407807349090E-06,
     * 3.07630759894132E-06, 3.87354743595132E-06, 4.87534530819265E-06,
     * 6.13932190719988E-06, 7.72696181044333E-06, 9.73039790393848E-06,
     * 1.22459098842999E-05, 1.54219178017462E-05, 1.94081087330868E-05,
     * 2.44429483141394E-05, 3.07586488196733E-05, 3.87404881736564E-05,
     * 4.87473381508339E-05, 6.14020175804902E-05, 7.72557456005821E-05,
     * 9.73193205963702E-05, 1.22434967068984E-04, 1.54247602338967E-04,
     * 1.94032567531162E-04, 2.44476408791169E-04, 3.07491837697500E-04,
     * 3.87483530664877E-04, 4.87280120178207E-04, 6.14131950758203E-04,
     * 7.72146864755635E-04, 9.73311869459581E-04, 1.22343779255314E-03,
     * 1.54242843814657E-03, 1.93818519482363E-03, 2.44390216121152E-03,
     * 3.06963949011437E-03, 3.87097274627399E-03, 4.85916079220826E-03,
     * 6.12753973006850E-03, 7.68490745325353E-03, 9.68823557105223E-03,
     * 1.21336164257623E-02, 1.52846174793362E-02, 1.90987890260545E-02,
     * 2.40156395862743E-02, 2.98915604118627E-02, 3.74477389774394E-02,
     * 4.62895386766490E-02, 5.75627504992849E-02, 7.02642576389800E-02,
     * 8.61063199481585E-02, 1.02635894673126E-01, 1.22138909449955E-01,
     * 1.38831147894879E-01, 1.55250983600674E-01, 1.58685462350600E-01,
     * 1.51999745985231E-01, 1.11620805416968E-01, 4.79636092085399E-02,
     *-6.59958907891261E-02,-1.75114294240188E-01,-2.72429084088151E-01,
     *-2.05512509493180E-01,-8.95229587687801E-03, 3.15188860127114E-01,
     * 2.37178375962863E-01,-1.67129454539411E-01,-3.89486718623514E-01,
     * 4.69395996090898E-01,-1.81704320197439E-01,-4.79333224758986E-02,
     * 1.22756324360416E-01,-1.16936501005281E-01, 9.10625148146839E-02,
     *-6.70167357658165E-02, 4.87676580356703E-02,-3.55323199801816E-02,
     * 2.59727315593579E-02,-1.90338807769771E-02, 1.39717419500800E-02,
     *-1.02658003799853E-02, 7.54695165040799E-03,-5.54984985828324E-03,
     * 4.08190272675497E-03,-3.00250008140396E-03, 2.20863797464272E-03,
     *-1.62471694082797E-03, 1.19519049147996E-03,-8.79224178055185E-04,
     * 6.46791147025135E-04,-4.75806021520080E-04, 3.50022446653354E-04,
     *-2.57490755236966E-04, 1.89421182988102E-04,-1.39346442298395E-04,
     * 1.02508856775501E-04,-7.54097078313750E-05, 5.54748565814861E-05,
     *-4.08096661134060E-05, 3.00210367411730E-05,-2.20848258928823E-05,
     * 1.62467609869699E-05,-1.19516440312774E-05, 8.79197554324395E-06,
     *-6.46795540161166E-06, 4.75817754903626E-06,-3.50009454613693E-06,
     * 2.57482374988796E-06,-1.89435683550768E-06, 1.39348990804983E-06,
     *-1.02494082372593E-06, 7.54124172617477E-07,-5.54879450596859E-07,
     * 4.08025014807835E-07,-3.00110564680988E-07, 2.20951460250530E-07,
     *-1.62526982476834E-07, 1.19516631393177E-07,-8.79214841634993E-08,
     * 6.46787588254737E-08,-4.75804279580232E-08, 3.50021732911949E-08,
     *-2.57490776709214E-08, 1.89421095480929E-08,-1.39346161721808E-08,
     * 1.02508924559328E-08,-7.54098963650564E-09, 5.54747061705591E-09,
     *-4.08095379127959E-09, 3.00212204736251E-09,-2.20848783108477E-09,
     * 1.62465696700589E-09,-1.19516631393177E-09, 8.79214841634989E-10,
     *-6.46787588254742E-10, 4.75804279580229E-10,-3.50021732911947E-10,
     * 2.57490776709216E-10,-1.89421095480930E-10, 1.39346161721808E-10,
     *-1.02508924559327E-10, 7.54098963650572E-11,-5.54747061705581E-11,
     * 4.08095379127958E-11,-3.00212204736254E-11, 2.20848783108476E-11,
     *-1.62465696700591E-11, 1.19516631393176E-11,-8.79214841634984E-12,
     * 6.46787588254741E-12,-4.75804279580233E-12, 3.50021732911941E-12,
     *-2.57490776709219E-12, 1.89421095480930E-12,-1.39346161721807E-12,
     * 1.02508924559326E-12,-7.54098963650558E-13, 5.54747061705595E-13,
     *-4.08095379127967E-13, 3.00212204736254E-13,-2.20848783108478E-13,
     * 1.62465696700586E-13,-1.19516631393176E-13, 8.79214841635006E-14,
     *-6.46787588254739E-14, 4.75804279580232E-14,-3.50021732911949E-14,
     * 2.57490776709211E-14,-1.89421095480929E-14, 1.39346161721806E-14,
     *-1.02508924559328E-14, 7.54098963650577E-15,-5.54747061705592E-15,
     * 4.08095379127955E-15,-3.00212204736252E-15, 2.20848783108471E-15,
     *-1.62465696700590E-15, 1.19516631393179E-15,-8.79214841635001E-16,
     * 6.46787588254735E-16,-4.75804279580229E-16, 3.50021732911938E-16/

      ncmt=nc+mt-1
      do it=1,mt
         t(it)=ta*q**(it-1)
         ft(it)=0.d0
      enddo
      do nn=1,ncmt
         n=-nc+nc0+nn
         omega=q**(-(n-1))/ta
         fomega=dsqrt(omega)*fun(omega)
         ita=max0(1,nn-nc+1)
         ite=min0(mt,nn)
         do it=ita,ite
            itn=nc-nn+it
            ft(it)=ft(it)-dimag(fomega)*h(itn)
         enddo
      enddo
      do it=1,mt
         ft(it)=dsqrt(8.d0*pi/t(it))*ft(it)
      enddo
      return
      end
            
c     ************************************************************

      subroutine frt10(ta,mt,t,ft,fun)
c     subroutine zur berechnung der fouriertransformierten der 
c     funktion fun(omega), definiert im function-unterprogramm 
c     fun(omega) die transformierte funktion wird an logarith-
c     misch aequidistanten stuetzstellen berechnet mit 10 
c     stuetzstellen pro dekade.

c     bedeutung der parameter:
c     ta = erste stuetzstelle im t-bereich (input)
c     mt = anzahl der stuetzstellen (input)
c     t = bereich der logarithmisch aequidistanten stuetzst.
c     (output)
c     ft = bereich der transformierten funktion (output)
c     fun = name des function-unterprogramms  zur berechnung der 
c     frequenzabhaengigen funktion (=external im rufenden 
c     programm)

      implicit real*8 (a-h,o-z)
      dimension h(80),t(*),ft(*)
      complex*16 fun,fomega
      data pi,q/3.141592653589793d0,1.258925411794167d0/
      data nc0,nc,(h(i),i=1,80)/40,80,
     * 2.59526236e-07, 3.66544843e-07, 5.17830795e-07, 7.31340622e-07,
     * 1.03322805e-06, 1.45918500e-06, 2.06161065e-06, 2.91137793e-06,
     * 4.11357863e-06, 5.80876420e-06, 8.20798075e-06, 1.15895083e-05,
     * 1.63778560e-05, 2.31228459e-05, 3.26800649e-05, 4.61329334e-05,
     * 6.52101085e-05, 9.20390575e-05, 1.30122935e-04, 1.83620431e-04,
     * 2.59656626e-04, 3.66311982e-04, 5.18141184e-04, 7.30717340e-04,
     * 1.03392184e-03, 1.45742714e-03, 2.06292302e-03, 2.90599911e-03,
     * 4.11471902e-03, 5.79042763e-03, 8.20004722e-03, 1.15192930e-02,
     * 1.63039133e-02, 2.28257757e-02, 3.22249222e-02, 4.47864328e-02,
     * 6.27329625e-02, 8.57059100e-02, 1.17418314e-01, 1.53632655e-01,
     * 1.97717964e-01, 2.28849849e-01, 2.40311038e-01, 1.65409220e-01,
     * 2.84701476e-03,-2.88016057e-01,-3.69097406e-01,-2.50107514e-02,
     * 5.71811256e-01,-3.92261572e-01, 7.63280044e-02, 5.16233994e-02,
     *-6.48012082e-02, 4.89047141e-02,-3.26936331e-02, 2.10539842e-02,
     *-1.33862549e-02, 8.47124695e-03,-5.35123972e-03, 3.37796651e-03,
     *-2.13174466e-03, 1.34513833e-03,-8.48749612e-04, 5.35531006e-04,
     *-3.37898780e-04, 2.13200109e-04,-1.34520273e-04, 8.48765787e-05,
     *-5.35535069e-05, 3.37899801e-05,-2.13200365e-05, 1.34520337e-05,
     *-8.48765949e-06, 5.35535110e-06,-3.37899811e-06, 2.13200368e-06,
     *-1.34520338e-06, 8.48765951e-07,-5.35535110e-07, 3.37899811e-07/

      ncmt=nc+mt-1
      do it=1,mt
         t(it)=ta*q**(it-1)
         ft(it)=0.d0
      enddo
      do nn=1,ncmt
         n=-nc+nc0+nn
         omega=q**(-(n-1))/ta
         fomega=dsqrt(omega)*fun(omega)
         ita=max0(1,nn-nc+1)
         ite=min0(mt,nn)
         do it=ita,ite
            itn=nc-nn+it
            ft(it)=ft(it)-dimag(fomega)*h(itn)
         enddo
      enddo
      do it=1,mt
         ft(it)=dsqrt(8.d0*pi/t(it))*ft(it)
      enddo
      return
      end

c     ************************************************************
  
      complex*16 function q(omega)
      implicit real*8 (a-h,o-z)
      common /param/n
      complex*16 qq,ii

      pi=4*datan(1.d0)
      ii=(0.d0,1.d0)
     
      period=2*pi/omega

c     unit step function excitation
     
      q=1/(2*pi*ii*omega)*qq(n,period)
      return
      end 

 
c     ************************************************************

      subroutine spline1(n,y,c)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     cubic spline interpolation for n equidistant ordinates y(1),  c
c     ...,y(n) with boundary conditions y"(1)=y"(n)=0 at the end-   c
c     points. the output array c(3,n) contains the interpolation    c
c     coefficients to be used in the subroutine spline2.            c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)
      dimension y(n),c(3,n)
      n1=n-1
      do i=2,n1
         c(1,i)=y(i+1)-2.d0*y(i)+y(i-1)
      enddo

      c(2,1)=0.d0
      c(3,1)=0.d0
      do i=2,n1
         p=4.d0+c(2,i-1)
         c(2,i)=-1.d0/p
         c(3,i)=(c(1,i)-c(3,i-1))/p
      enddo

      c(1,n)=0.d0
      do ii=2,n1
         i=n+1-ii
         c(1,i)=c(2,i)*c(1,i+1)+c(3,i)
      enddo
      c(1,1)=0.d0

      do i=1,n1
         c(2,i)=y(i+1)-y(i)-c(1,i+1)+c(1,i)
         c(3,i)=y(i)-c(1,i)
      enddo
      c(3,n)=y(n)
      return
      end

c     ***************************************************************

      subroutine spline2(n,xint,x1,x2,c,yint)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     this subroutine calculates for the abscissa xint the func-    c
c     tion yint from the coefficients c(3,n) determined in spline1. c
c     x1 and x2 are the abscissae of y(1) and y(n), respectively.   c
c     outside this range yint is linearly extrapolated.             c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)
      dimension c(3,n)
      h=(x2-x1)/dfloat(n-1)
      if (xint.lt.x1) goto 10
      if (xint.ge.x2) goto 20

      u=(xint-x1)/h
      i=1+int(u)
      p=u-i+1
      q=1.d0-p
      yint=c(1,i)*q**3+c(1,i+1)*p**3+c(2,i)*p+c(3,i)
      return

   10 p=(xint-x1)/h
      yint=c(2,1)*p+c(3,1)
      return

   20 p=(xint-x2)/h
      yint=c(2,n-1)*p+c(3,n)
      return
      end
  
c     ************************************************************

      subroutine filt(nc,fc,nd,dst,ext,inn)
      implicit real*8 (a-h,o-z)
      real*8 inn
!      dimension fc(0:nc),dst(0:nd),ext(0:nd),inn(0:nd)               ! BIG BUG!!! Either C or Fortran order, not both!
      dimension fc(0:nc-1),dst(0:nd-1),ext(0:nd-1),inn(0:nd-1)

      c0=1.d0/(1.d0+fc(0))

      do j=0,nd-1
c         if (j/10.eq.j/10.0) 
!         write(6,*)"filtering ih =",j
         ext(j)=dst(j)
         if (j.gt.0) then 
            kmax=min0(j,nc)        
            do k=1,kmax
               ext(j)=ext(j)-fc(k)*ext(j-k)
            enddo
         endif 
         ext(j)=c0*ext(j)
         inn(j)=dst(j)-ext(j)
!         write(6,*)"filtering ih =",j, dst(j), ext(j), inn(j)
      enddo 
      return
      end
          
c     ************************************************************
!      call filt(nc,fc,ndst,dstvec,est,ist)
