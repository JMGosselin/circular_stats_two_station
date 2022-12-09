c Copyright 2022 by Jeremy M. Gosselin
c Department of Geoscience
c University of Calgary

c This program is free software: 
c you can redistribute it and/or modify it under the terms of the 
c GNU General Public License as published by the Free Software 
c Foundation, either version 3 of the License, or (at your option) 
c any later version.

c This program is distributed in the hope that it will be useful, but 
c WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
c or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
c for more details.

c You should have received a copy of the GNU General Public License along 
c with this program. 
c If not, see <https://www.gnu.org/licenses/>.



c========================================================================

      SUBROUTINE SWD_calc(nlay,nfreq,H,Vp,Vs,rho,
     &                    freqs,vels,nmodes,group,ierr)

c========================================================================
      implicit none
      integer nlay, nfreq, nmodes, group, ierr, ii
      real H(nlay),Vp(nlay),Vs(nlay),rho(nlay),freqs(nfreq),vels(nfreq)
 
      real omega, c, u, ek, y0(3), yij(15), cmn, cmx, dc, tol 
      integer ier(0:3), itr

      cmn = 0.1
      cmx = 5.0
      dc  = 0.01
      tol = 0.001
      itr = 100

      do ii=1,nfreq 
            omega = freqs(ii)*2.*3.14159 ! Convert to angular frequency
                       
            call raydsp(H,rho,Vp,Vs,nlay,omega,
     &                  cmn,cmx,dc,tol,itr,c,u,ek,y0,yij,ier)
           
            if (ier(0) .ne. 0) then
                  ierr = 1
                  return
            endif

            if (group .eq. 1) then
                  vels(ii) = u
            elseif (group .eq. 0) then
                  vels(ii) = c
            endif

            if ((vels(ii) .lt. 1.5) .or. (vels(ii) .gt. 10.0)) then
               ierr = 1
               return
            endif  
      enddo

      ierr = 0

      return
      end



c=============================================================================      
      Subroutine raydsp(h,rho,vp,vs,l,w,cmn,cmx,
     &   dc,tol,itr,c,u,ek,y0,yij,ier)
c=============================================================================
c
c     rayleigh wave dispersion interpolation
c
c     input
c       diffeq : subroutine to integrate eq. of motion
c              rayrkg or raymrx
c       h    : layer thickness.  h(l) is arbitrary
c       rho  : density
c       vp   : compressional wave velocity
c       vs   : shear wave velocity
c       l    : number of layers including the bottom half space
c       w    : angular frequency
c       cmn  : lower limit of phase velocity
c       cmx  : upper limit of phase velocity
c       dc   : increment of phase velocity
c       tol  : relative accuracy of phase velocity
c       itr  : maximum number of iterations
c     output
c       c    : phase velocity
c       u    : group velocity by differentiation
c       ek   : kinetic energy by differentiation, 2*w**2*i1
c       y0   : surface values of eigenfunction
c              y0(1) ; y1 (scale factor)
c              y0(2) ; y2/abs(y1), dispersion function
c              y0(3) ; y3/y1
c       yij  : surface values of yij (compounds),
c              c*dyij/dc, and w*dyij/dw
c       ier(0):   return code
c                 =-2 ; iteration error       (rua 17/03/95)
c                       Either l<0 or w<0 or c<0
c                        ier(1): i, layer where error occured   
c                        ier(2): rho(i), densiy in bad layer
c                        ier(3): vp(i), p-wave velocity in bad layer
c                 =-1 ; raydsp input error    (rua 17/03/95)
c                       Either l<0 or w<0 or c<0
c                        ier(1): l, # of layers   
c                        ier(2): w, angular frequency
c                        ier(3): c, phase velocity
c                 = 0 ; no error
c                 = 1 ; slow convergence
c                 = 2 ; root not found
c        When ier(0)>=0, ier(1)=ier(2)=ier(3)=0 (rua 17/03/95)  
c              
c
c     subroutine : diffeq (raymrx or rayrkg)
c
c     disper-80
c
c
c******************************************************************
c
c
c
      dimension  h(l),rho(l),vp(l),vs(l),y0(3),yij(15),ier(0:3)
c
c     initialization
c
      if( l.le.2 .or. w.le.0. )  go to  90
      tol1 = (1. + tol) - 1.
      ier1 = 0
c
c =================================================================
c
c     Inserted by rua 17/03/95
c
      ier(1)=0
      ier(2)=0
      ier(3)=0
c
c =================================================================
c
c     Removed by rua 17/03/95
c
csed      write(6,1)  w
csed    1   format(/23x,'rayleigh wave'//7x,'w'/
csed     *         1pe18.6//7x,'c',17x,'y2',16x,'y3')
c
c =================================================================
c
      c3   = cmn
      if( c3.le.0. )  go to  92
	call  raymrx(h,rho,vp,vs,l,w,c3,1,u,ek,y0,yij,ier)

      if( ier(0).ne.0 )  then
	return
	endif

      f3   = y0(2)
c
c =================================================================
c
c     Removed by rua 17/03/95
c
csed      write(6,2)  c3,f3,y0(3)
csed    2   format(1p3e18.6)
c
c =================================================================
c
      if( f3.eq.0. .and. tol1.gt.0. )  go to  12
c
c     find a zero-cross
c
      kx   = 0
      if( dc.ne.0. )  kx = (cmx - cmn)/dc + 0.1
      if( kx.le.0 )  go to  4
c
      do  3  k=1,kx
        c1   = c3
        f1   = f3
        c3   = cmn + float(k)*dc
        if( c3.le.0. )  go to  92
        call  raymrx(h,rho,vp,vs,l,w,c3,1,u,ek,y0,yij,ier)

       if( ier(0).ne.0 )  then
	return
	endif

        f3   = y0(2)
c
c =================================================================
c
c     Removed by rua 17/03/95
c
csed        write(6,2)  c3,f3,y0(3)
c
c =================================================================
c
        if( tol1.le.0. )  go to  3
        if( f3*f1.le.0. )  go to  5
    3 continue
c
    4 ier(0)  = 2
      if( tol1.le.0. )  return
      go to  94
c
c     interpolation
c
    5 if( f3.eq.0. )  go to  12
      c2   = c3
      f2   = f3
      e    = c1 - c2
      d    = e*0.5
      c3   = c2 + d
      kx   = max0(1,itr)
c
      do  10  k=1,kx
        call  raymrx(h,rho,vp,vs,l,w,c3,1,u,ek,y0,yij,ier)

        if(ier(0) .ne. 0) return
      
        f3   = y0(2)
c
c =================================================================
c
        if( f3*f2.le.0. )  go to  6
          ff   = c1
          c1   = c2
          c2   = ff
          ff   = f1
          f1   = f2
          f2   = ff
    6   if( abs(f3).le.abs(f2) )  go to  7
          ff   = c2
          c2   = c3
          c3   = ff
          ff   = f2
          f2   = f3
          f3   = ff
    7   e    = c2 - c3
        if( f3.eq.0. )  go to  12
        tolc = c3*tol1
        dd   = d
c
c       inverse quadratic interpolation
c
        f32  = f3/f2
        f31  = f3/f1
        f21  = f2/f1
        q    = f32*(e*(1. - f31) + f21*(f31 - f21)*(c1 - c3))
        s    = (f21 - 1.)*(f32 - 1.)*(f31 - 1.)
c
c       test range
c
        if( q.lt.0. )  s =-s
        q    = abs(q)
        if( q.lt.(e*s-abs(tolc*s)) )  go to  8
c
c       linear interpolation
c
        d    = e*f32/(f32 - 1.)
        go to  9
c
    8   d    = q/s
c
c       test convergence
c
    9   c1   = c2
        f1   = f2
        c2   = c3
        f2   = f3
        c3   = c2 + d
        if( abs(e).le.tolc )  go to  12
        if( abs(d).gt.tolc )  go to  10
        c3   = c2 + sign(tolc,d)
        if(abs(dd).gt.tolc )  go to  10
c
c       bisection
c
        d    = e*0.5
        c3   = c2 + d
c
   10 continue
c
c     slow convergence
c
c
c =================================================================
c
c     Removed by rua 17/03/95
c
csed      write(6,111)  kx
  111   format(20x,5('?'),3x,'(raydsp)   slow conv. after',
     *         i5,' iterations',3x,5('?'))
c
c =================================================================
c
   11 ier1 = 1
c
c     root is found
c
   12 call  raymrx(h,rho,vp,vs,l,w,c3,3,u,ek,y0,yij,ier)
      ier(0)  = ier1
      f3   = y0(2)

c
c =================================================================
c
      c    = c3
      return
      
c =================================================================
c     input error

cJMG   90 write(6,91)  l,w
   90 continue
   91 format(20x,5('?'),3x,'(raydsp)   input error   l =',
     *   i5,3x,'w =',1pe13.6,', should be gt 0.',3x,5('?'))
      ier(0)  =-1
      ier(1)  = l
      ier(2)  = w
      ier(3)  = c
      return
c
c =================================================================
c

cJMG   92 write(6,93)  c3
   92 continue
   93 format(20x,5('?'),3x,'(raydsp)   input error   c =',
     *   1pe13.6,', should be gt 0.',3x,5('?'))
	    return
c
c =================================================================
c
c     no root

csed   94 write(6,95)
   94 continue
   95 format(20x,5('?'),3x,'(raydsp)   root not found',3x,5('?'))
      return
c
c =================================================================
      end
      
      
      
      
      


      subroutine  raymrx(h,rho,vp,vs,l,w,c,ig,u,ek,y0,yij,ier)

      dimension  h(l),rho(l),vp(l),vs(l),y0(3),yij(15),ier(0:3)
      data  eps/1.e-30/
c
c     define sinh(x)/x and (cosh(x)-sinh(x)/x)/x**2
c
      sh0(x) = 0.9999997 + x*(0.1666667 + x*(0.0083361 + x*0.0001984))
      sh1(x) = 0.3333333 + x*(0.0333333 + x*(0.0011907 + x*0.0000220))
c
c     for double precision use
c
c     sh0(x) = 1.0d0 + x*(1.6666 66666 66666 7d-1
c    *       + x*(8.33 33333 33334 0d-3 + x*(1.9 84126 98412 7d-4
c    *       + x*(2.7557 31918 9d-6 + x*(2.50 12108 4d-8
c    *       + x*(1.60596 1d-10 + x*7.64 7d-13))))))
c     sh1(x) = 3.3333 33333 33333 3d-1 + x*(3.333 33333 33333 3d-2
c    *       + x*(1.19 04761 90476 2d-3 + x*(2.20458 55379 2d-5
c    *       + x*(2.505 21083 7d-7 + x*(1.9 27085 3d-9
c    *       + x*(1.0706 3d-11 + x*4.50d-14))))))
c
c     initial value
c
      if( l.le.0 .or. w.le.0. .or. c.le.0. )  go to  90
      i    = l
      if( rho(i).le.0. .or. c.ge.vp(i) )  go to  92
      ier(0)  = 0
c
c =================================================================
c
c     Inserted by rua 17/03/95
c
      ier(1)  = 0
      ier(2)  = 0
      ier(3)  = 0
c
c =================================================================
c
      wn   = w/c
      igg  = max0(1,min0(3,ig))
      ro   = rho(i)
      sv   = vs(i)
      cp   = c/vp(i)
      raa  = (1. + cp)*(1. - cp)
      ra   = sqrt(raa)
c
c     liquid bottom
c
      if( sv.gt.0. )  go to  1
      y1   = ra*eps
      y2   =-ro*eps
c
      if( igg.lt.2 )  go to  2
      y3   =-cp**2*eps/ra
      y4   = 0.
c
      if( igg.le.2 )  go to  2
      y5   = 0.
      y6   = 0.
      go to  2
c
c     solid bottom
c
    1 if( c.ge.sv )  go to  92
      cs   = c/sv
      rbb  = (1. + cs)*(1. - cs)
      rb   = sqrt(rbb)
      rg   = 2.*ro/cs**2
      y3   =-ra*eps
      y4   =-rb*eps
      y2   =-eps*(cp**2*rbb + cs**2)/(ro*(ra*rb + 1.))
      y1   = rg*y2 + eps
      y5   =-rg*(y1 + eps) + ro*eps
c
      if( igg.lt.2 )  go to  2
      y8   = eps*cp**2/ra
      y9   = eps*cs**2/rb
      y7   =-(rb*y8 + ra*y9)/ro
      y6   = rg*(y7 - y2-y2)
      ya   = rg*(y1+y1 + eps+eps - y6)
c
      if( igg.le.2 )  go to  2
      yb   = 0.
      yc   = 0.
      yd   = 0.
      ye   = 0.
      yf   = 0.
c
c     integrate upward
c
    2 if( l.le.1 )  go to  12
      do  11  ii=2,l
        i    = i - 1
        ro   = rho(i)
        pv   = vp(i)
        sv   = vs(i)
        if( pv.le.0. .or. ro.le.0. .or. h(i).le.0. )  go to  92
        if((sv.eq.0. .and. vs(i+1).eq.0.) .or.
     *     (sv.ne.0. .and. vs(i+1).ne.0.))  go to  4
c
c       solid to liquid boundary
c
        if( sv.ne.0. )  go to  3
        y0(3) =-y1/y3
        y1   = y3
        y2   = y5
c
        if( igg.lt.2 )  go to  4
        y3   = y8
        y4   = ya
c
        if( igg.le.2 )  go to  4
        y5   = yd
        y6   = yf
        go to  4
c
c       liquid to solid boundary
c
    3   y9   = y4
        y4   = y2
        y2   = y1
        y7   = y3
        yc   = y5
        ye   = y6
        y1   = 0.
        y3   = 0.
        y5   = 0.
c
        if( igg.lt.2 )  go to  4
        y6   = 0.
        y8   = 0.
        ya   = 0.
c
        if( igg.le.2 )  go to  4
        yb   = 0.
        yd   = 0.
        yf   = 0.
c
    4   r2   = 1./ro
        cp   = c/pv
        raa  = (1. + cp)*(1. - cp)
        cs   = 0.
        if( sv.gt.0. )  cs = c/sv
        rbb  = (1. + cs)*(1. - cs)
        hk   = h(i)*wn
        hkk  = hk**2
        xx   = raa*hkk
        one  = 1.
c
c       sinh(x)/x
c
        do  9  k=1,2
          cha  = chb
          sha  = shb
          dha  = dhb
          ax   = abs(xx)
          if( ax.le.1. )  go to  7
          ax   = sqrt(ax)
          if( xx.le.0. )  go to  5
          if( ax.gt.100. )  one = 0.
          if( ax.le.100. )  one = one/cosh(ax)
          chb  = 1.
          shb  = tanh(ax)/ax
          go to  6
    5     chb  = cos(ax)
          shb  = sin(ax)/ax
    6     if( igg.ge.2 )  dhb = (hkk/xx)*(chb - shb)
          go to  8
    7     shb  = sh0(xx)
          chb  = 1. + 0.5*xx*sh0(xx*0.25)**2
          if( igg.ge.2 )  dhb = sh1(xx)*hkk
    8     xx   = hkk*rbb
          shb  = hk*shb
    9   continue
c
c     layer matrices
c
c       liquid layer
c
        if( sv.gt.0. )  go to  10
        b11  = cha
        b12  =-raa*sha*r2
        b21  =-ro*sha
c
        z1   = y1
        z2   = y2
        y1   = b11*z1 + b12*z2
        y2   = b21*z1 + b11*z2
c
        if( igg.lt.2 )  go to  11
        c11  =-hk*sha
        c12  = (cp**2*sha + hk*cha)*r2
        c21  = ro*(sha + hk*dha)
c
        z3   = y3
        y3   = b11*z3 + b12*y4 + c11*z1 + c12*z2
        y4   = b21*z3 + b11*y4 + c21*z1 + c11*z2
c
        if( igg.le.2 )  go to  11
        w11  = hk*raa*sha
        w12  =-hk*raa*cha*r2
        w21  =-hk*ro*cha
c
        z5   = y5
        y5   = b11*z5 + b12*y6 + w11*z1 + w12*z2
        y6   = b21*z5 + b11*y6 + w21*z1 + w11*z2
        go to  11
c
c       solid layer
c
   10   g1   = 2./cs**2
        rg   = g1*ro
        r4   = rg - ro
        e1   = cha*chb
        e2   = e1 - one
        e3   = sha*shb
        e5   = sha*chb
        e6   = shb*cha
        f1   = e2 - e3
        f2   = r2*f1
        f3   = g1*f1 + e3
        b33  = e1
        b34  = raa*e3
        b43  = rbb*e3
        b25  =-r2*(f2 + r2*(e2 - raa*b43))
        b15  = rg*b25 + f2
        b16  =-rg*b15 - f3
        b22  = b16 + e1
        b12  = rg*b16 - r4*f3
        b52  =-rg*b12 + r4*(rg*f3 + r4*e3)
        b23  = r2*(e5 - rbb*e6)
        b13  = rg*b23 - e5
        b42  =-rg*b13 + r4*e5
        b24  = r2*(e6 - raa*e5)
        b14  = rg*b24 - e6
        b32  =-rg*b14 + r4*e6
        b11  = one - b16-b16
        b21  = b15+b15
        b31  = b14+b14
        b41  = b13+b13
        b51  = b12+b12
c
        z1   = y1
        z2   = y2
        z3   = y3
        z4   = y4
        z5   = y5
        y1   = b11*z1 + b12*z2 + b13*z3 + b14*z4 + b15*z5
        y2   = b21*z1 + b22*z2 + b23*z3 + b24*z4 + b25*z5
        y3   = b31*z1 + b32*z2 + b33*z3 + b34*z4 + b24*z5
        y4   = b41*z1 + b42*z2 + b43*z3 + b33*z4 + b23*z5
        y5   = b51*z1 + b52*z2 + b42*z3 + b32*z4 + b22*z5
c
        if( igg.lt.2 )  go to  11
        raac =-2.*cp*cp
        rbbc =-2.*cs*cs
        r3c  =-rg-rg
        e1c  =-hk*(e5 + e6)
        e3c  =-e3-e3 - hk*(dha*shb + dhb*sha)
        e5c  =-e5 - hk*(dha*chb + e3)
        e6c  =-e6 - hk*(dhb*cha + e3)
        f1c  = e1c - e3c
        f2c  = r2*f1c
        f3c  = g1*(f1c - f1-f1) + e3c
        c33  = e1c
        c34  = raa*e3c + raac*e3
        c43  = rbb*e3c + rbbc*e3
        c25  =-r2*(f2c + r2*(e1c - raa*c43 - raac*b43))
        c15  = rg*c25 + r3c*b25 + f2c
        c16  =-rg*c15 - r3c*b15 - f3c
        c22  = c16 + e1c
        c12  = rg*c16 + r3c*(b16 - f3) - r4*f3c
        c52  =-rg*c12 + r4*(rg*f3c + r4*e3c)
     *       + r3c*(-b12 + (rg + r4)*f3 + 2.*r4*e3)
        c23  = r2*(e5c - rbb*e6c - rbbc*e6)
        c13  = rg*c23 + r3c*b23 - e5c
        c42  =-rg*c13 + r4*e5c + r3c*(e5 - b13)
        c24  = r2*(e6c - raa*e5c - raac*e5)
        c14  = rg*c24 - e6c + r3c*b24
        c32  =-rg*c14 + r4*e6c + r3c*(e6 - b14)
        c11  =-c16-c16
        c21  = c15+c15
        c31  = c14+c14
        c41  = c13+c13
        c51  = c12+c12
c
        z6   = y6
        z7   = y7
        z8   = y8
        z9   = y9
        y6   = b11*z6 + b12*z7 + b13*z8 + b14*z9 + b15*ya
     *       + c11*z1 + c12*z2 + c13*z3 + c14*z4 + c15*z5
        y7   = b21*z6 + b22*z7 + b23*z8 + b24*z9 + b25*ya
     *       + c21*z1 + c22*z2 + c23*z3 + c24*z4 + c25*z5
        y8   = b31*z6 + b32*z7 + b33*z8 + b34*z9 + b24*ya
     *       + c31*z1 + c32*z2 + c33*z3 + c34*z4 + c24*z5
        y9   = b41*z6 + b42*z7 + b43*z8 + b33*z9 + b23*ya
     *       + c41*z1 + c42*z2 + c43*z3 + c33*z4 + c23*z5
        ya   = b51*z6 + b52*z7 + b42*z8 + b32*z9 + b22*ya
     *       + c51*z1 + c52*z2 + c42*z3 + c32*z4 + c22*z5
c
        if( igg.le.2 )  go to  11
        e1w  = hk*(raa*e5 + rbb*e6)
        e3w  = hk*(e5 + e6)
        e5w  = hk*(e1 + b43)
        e6w  = hk*(e1 + b34)
        f1w  = e1w - e3w
        f2w  = r2*f1w
        f3w  = g1*f1w + e3w
        w33  = e1w
        w34  = raa*e3w
        w43  = rbb*e3w
        w25  =-r2*(f2w + r2*(e1w - raa*w43))
        w15  = rg*w25 + f2w
        w16  =-rg*w15 - f3w
        w22  = w16 + e1w
        w12  = rg*w16 - r4*f3w
        w52  =-rg*w12 + r4*(rg*f3w + r4*e3w)
        w23  = r2*(e5w - rbb*e6w)
        w13  = rg*w23 - e5w
        w42  =-rg*w13 + r4*e5w
        w24  = r2*(e6w - raa*e5w)
        w14  = rg*w24 - e6w
        w32  =-rg*w14 + r4*e6w
        w11  =-w16-w16
        w21  = w15+w15
        w31  = w14+w14
        w41  = w13+w13
        w51  = w12+w12
c
        zb   = yb
        zc   = yc
        zd   = yd
        ze   = ye
        yb   = b11*zb + b12*zc + b13*zd + b14*ze + b15*yf
     *       + w11*z1 + w12*z2 + w13*z3 + w14*z4 + w15*z5
        yc   = b21*zb + b22*zc + b23*zd + b24*ze + b25*yf
     *       + w21*z1 + w22*z2 + w23*z3 + w24*z4 + w25*z5
        yd   = b31*zb + b32*zc + b33*zd + b34*ze + b24*yf
     *       + w31*z1 + w32*z2 + w33*z3 + w34*z4 + w24*z5
        ye   = b41*zb + b42*zc + b43*zd + b33*ze + b23*yf
     *       + w41*z1 + w42*z2 + w43*z3 + w33*z4 + w23*z5
        yf   = b51*zb + b52*zc + b42*zd + b32*ze + b22*yf
     *       + w51*z1 + w52*z2 + w42*z3 + w32*z4 + w22*z5
c
   11 continue
c
c     normal exit
c
   12 if( sv.ne.0. )  go to  13
      y0(1) = y1
      y0(2) = y2/abs(y1)
      yij(1) = y1
      yij(2) = y2
c
      if( igg.lt.2 )  return
      yij(3) = y3
      yij(4) = y4
c
      if( igg.le.2 )  return
      yij(5) = y5
      yij(6) = y6
      u    = c*y4/(y4 + y6)
      ek   =-(w*c*c/u)*(y4/y1)
      return
c
   13 y0(1) = y3
      y0(2) = y5/abs(y3)
      y0(3) =-y1/y3
      yij( 1) = y1
      yij( 2) = y2
      yij( 3) = y3
      yij( 4) = y4
      yij( 5) = y5
c
      if( igg.lt.2 )  return
      yij( 6) = y6
      yij( 7) = y7
      yij( 8) = y8
      yij( 9) = y9
      yij(10) = ya
c
      if( igg.le.2 )  return
      yij(11) = yb
      yij(12) = yc
      yij(13) = yd
      yij(14) = ye
      yij(15) = yf
      u    = c*ya/(ya + yf)
      ek   =-(w*c*c/u)*(ya/y3)
      return
c
c     input error
c
c =================================================================
c
c     Removed by rua 17/03/95
c
   90 continue
c
      ier(0)  =-1
c
c =================================================================
c
c     Inserted by rua 17/03/95
c
      ier(1)  = l
      ier(2)  = w
      ier(3)  = c
c
c =================================================================
c
      return
c
   92 if( i.le.0 )  i = 1
c
c =================================================================
c
      ier(0)  =-2
c
c =================================================================
c
c     Inserted by rua 17/03/95
c
      ier(1)  = i
      ier(2)  = rho(i)
      ier(3)  = vp(i)
c
c =================================================================
c
      return
      end

      

