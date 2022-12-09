
! Copyright 2022 by Jeremy M. Gosselin
! Department of Geoscience
! University of Calgary

! This program is free software: 
! you can redistribute it and/or modify it under the terms of the 
! GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) 
! any later version.

! This program is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
! for more details.

! You should have received a copy of the GNU General Public License along 
! with this program. 
! If not, see <https://www.gnu.org/licenses/>.



!c========================================================================

      SUBROUTINE SWD_calc_CPS_rayleigh(nlay,nfreq,HH,Vp,Vs,rho,&
                         freqs,vels,nmodes,group,ierr)

!c========================================================================
      implicit none
      integer nlay, nfreq, nmodes, group, ierr
      real HH(nlay),Vp(nlay),Vs(nlay),rho(nlay),freqs(nfreq)
      real vels(nfreq)
!c---------------------------------------------------------------- 
      integer idispl, idispr, nsph, mmaxin
      real(4) din(1000),ain(1000),bin(1000),rhoin(1000),rat(1000) 
      real(4) qbinv(1000)

      real(4) ddc, sone, h, t(512)
      real(8) c(1000), cb(1000), U(1000)

      integer iunit, ifunc, kmax, mode, igr, iq, is, ie
!c----------------------------------------------------------------
      integer ii
      real depth

      idispl   = 0 
      idispr   = nfreq
      nsph     = 1
      mmaxin   = nlay

      depth = 0.0
      do ii = 1,nlay
         din(ii) = HH(ii)
         ain(ii) = Vp(ii)
         bin(ii) = Vs(ii)
         rhoin(ii) = rho(ii)

         if (Vs(ii) .eq. 0.0) then
            rat(ii) = 0.5
         else
            rat(ii) = ( 2.*(Vp(ii)/Vs(ii))**2 - 1. )/ &
                     ( 4.*(Vp(ii)/Vs(ii))**2 - 1. )  
         endif

         depth = depth + HH(ii)
         qbinv(ii) = 0.0

         if (depth .lt. 80.0) then
            qbinv(ii) = 1./600.
         else if (depth .lt. 220.0) then
            qbinv(ii) = 1./80.
         else
            qbinv(ii) = 1./143.
         endif      
      enddo

      iunit = 0
!c!    ifunc - 1 Love,  2 Rayleigh
      ifunc = 2
      kmax = nfreq
      mode = 1
      ddc = 0.01
      sone = 2.0
!c!    igr - must by > 0 to calculate information for group velocities
      if (group .eq. 1) then
         igr = 1
      else
         igr = 0
      endif
      h = 0.001
      t(1:nfreq) = 1./freqs(1:nfreq)
      iq = 1
      is = 1
      ie = nfreq  

!!$OMP CRITICAL        

      call ph_vel(idispl,idispr,nsph, &
                  mmaxin,din,ain,bin,rhoin,rat,iunit,qbinv, &
                  ifunc, &
                  kmax,mode,ddc,sone,igr,h, &
                  t, &
                  iq,is,ie, &
                  c,cb)

!!$OMP END CRITICAL        


      if (sum(c(1:nfreq)) .lt. 0.0001) then
         ierr = 1
         return
      endif

      if (group .eq. 1) then
         do ii=1,nfreq
            vels(ii) =1./ &
           (1./c(ii) + ( 1./c(ii)-1./cb(ii)  )/2./h )
         enddo
      elseif (group .eq. 0) then
         do ii=1,nfreq
            vels=c(ii)
         enddo
      endif
         
      ierr = 0
      return

      END SUBROUTINE SWD_calc_CPS_rayleigh
!c========================================================================



!c========================================================================

      SUBROUTINE SWD_calc_CPS_love(nlay,nfreq,HH,Vp,Vs,rho, &
                                   freqs,vels,nmodes,group,ierr)

!c========================================================================
      implicit none
      integer nlay, nfreq, nmodes, group, ierr
      real HH(nlay),Vp(nlay),Vs(nlay),rho(nlay),freqs(nfreq)
      real vels(nfreq)
!c---------------------------------------------------------------- 
      integer idispl, idispr, nsph, mmaxin
      real(4) din(1000),ain(1000),bin(1000),rhoin(1000),rat(1000) 
      real(4) qbinv(1000)

      real(4) ddc, sone, h, t(512)
      real(8) c(1000), cb(1000), U(1000)

      integer iunit, ifunc, kmax, mode, igr, iq, is, ie
!c----------------------------------------------------------------
      integer ii
      real depth

      idispl   = nfreq 
      idispr   = 0
      nsph     = 1
      mmaxin   = nlay

      depth = 0.0
      do ii = 1,nlay
         din(ii) = HH(ii)
         ain(ii) = Vp(ii)
         bin(ii) = Vs(ii)
         rhoin(ii) = rho(ii)

         if (Vs(ii) .eq. 0.0) then
            rat(ii) = 0.5
         else
            rat(ii) = ( 2.*(Vp(ii)/Vs(ii))**2 - 1. )/ &
                      ( 4.*(Vp(ii)/Vs(ii))**2 - 1. )  
         endif

         depth = depth + HH(ii)
         qbinv(ii) = 0.0

         if (depth .lt. 80.0) then
            qbinv(ii) = 1./600.
         else if (depth .lt. 220.0) then
            qbinv(ii) = 1./80.
         else
            qbinv(ii) = 1./143.
         endif      
      enddo

      iunit = 0
!c!    ifunc - 1 Love,  2 Rayleigh
      ifunc = 1
      kmax = nfreq
      mode = 1
      ddc = 0.01
      sone = 2.0
!c!    igr - must by > 0 to calculate information for group velocities
      if (group .eq. 1) then
         igr = 1
      else
         igr = 0
      endif
      h = 0.001
      t(1:nfreq) = 1./freqs(1:nfreq)
      iq = 1
      is = 1
      ie = nfreq  

      call ph_vel(idispl,idispr,nsph, &
                  mmaxin,din,ain,bin,rhoin,rat,iunit,qbinv, &
                  ifunc, &
                  kmax,mode,ddc,sone,igr,h, &
                  t, &
                  iq,is,ie, &
                  c,cb)

      if (sum(c(1:nfreq)) .lt. 0.0001) then
         ierr = 1
         return
      endif

      if (group .eq. 1) then
         do ii=1,nfreq
            vels(ii) =1./ &
           ( 1./c(ii) + ( 1./c(ii)-1./cb(ii)  )/2./h )
         enddo
      elseif (group .eq. 0) then
         do ii=1,nfreq
            vels=c(ii)
         enddo
      endif
         
      ierr = 0
      return

      END SUBROUTINE SWD_calc_CPS_love
! c========================================================================
! c========================================================================
! c========================================================================
! c========================================================================
! c========================================================================











