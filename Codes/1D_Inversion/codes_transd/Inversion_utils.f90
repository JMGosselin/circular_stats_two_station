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



!========================================================================
   function Cauchy(idum,xo,b)
!  Generate random number from Cauchy distribution
!  xo is the location of the peak of the distribution (statistical median)
!  b is the half width at half maximum
   implicit none
   integer idum
   real xo,b
   real Cauchy
   real pi
   real ran1
   pi = 3.14159
   
   Cauchy = xo + (b*tan(pi*(ran1(idum)-0.5)))
   return
   end
!========================================================================

!========================================================================
   function Gasdev(idum)
!--Generate Gaussian-distributed random variable.
!  (C) Copr. 1986-92 Numerical Recipes Software 0H21.
!  USES ran1

   implicit none
   INTEGER idum
   REAL gasdev
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran1
   SAVE iset,gset     
   DATA iset/0/    
   external ran1
   
   if (iset.eq.0) then
1       v1=2.D0*ran1(idum)-1.D0
     v2=2.D0*ran1(idum)-1.D0
     rsq=v1**2+v2**2
     if(rsq.ge.1.D0.or.rsq.eq.0.D0) goto 1
     fac=sqrt(-2.D0*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
   else
     gasdev=gset
     iset=0
   endif      

   return
   end
!========================================================================

!========================================================================
   function Ran1(idum)
!--Generates a uniform random variable on [0,1].
!  (C) Copr. 1986-92 Numerical Recipes Software 0H21.

   implicit none
   INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
   REAL ran1,AM,EPS,RNMX
   PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
   PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   INTEGER j,k,iv(NTAB),iy
   SAVE iv,iy
   DATA iv /NTAB*0/, iy /0/
   
   if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do 11 j=NTAB+8,1,-1
       k=idum/IQ
       idum=IA*(idum-k*IQ)-IR*k
       if (idum.lt.0) idum=idum+IM
       if (j.le.NTAB) iv(j)=idum
11      continue
     iy=iv(1)
   endif
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   j=1+iy/NDIV
   iy=iv(j)
   iv(j)=idum
   ran1=min(AM*iy,RNMX)
   
   return
   end function
!========================================================================

!========================================================================
   subroutine linspace(d1,d2,n,grid)
      implicit none
      integer, intent(in) :: n
      real, intent(in) :: d1, d2
      real, dimension(n), intent(out) :: grid
      integer :: indxi

      grid(1) = d1
      do indxi= 0,n-2
         grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
      enddo
      grid(n) = d2

   end subroutine
!========================================================================

!========================================================================
   subroutine period_space(d1,d2,n,grid)
      implicit none
      integer, intent(in) :: n
      real, intent(in) :: d1, d2
      real, dimension(n), intent(out) :: grid
      real, dimension(n) :: grid2
      integer :: indxi

      call linspace(1./d1,1./d2,n,grid2)

      do indxi= 1,n
         grid(indxi) = 1./grid2(indxi)
      enddo

   end subroutine
!========================================================================

!========================================================================
   subroutine log10space(d1,d2,n,grid)
      implicit none
      integer, intent(in) :: n
      real, intent(in) :: d1, d2
      real, dimension(n), intent(out) :: grid
      real, dimension(n) :: grid2
      integer :: indxi

      call linspace(log10(d1),log10(d2),n,grid2)

      do indxi= 1,n
         grid(indxi) = 10**grid2(indxi)
      enddo

   end subroutine
!========================================================================

!========================================================================
   function factorial (n) result (res)
      implicit none
      integer, intent (in)  :: n
      real*8                :: res
      integer :: i
      res = 1.0
      do i = 1,n
         res = res*i
      enddo
   end function factorial
!========================================================================

!========================================================================
   function choose (n, k) result (res)
      implicit none
      integer, intent (in)  :: n
      integer, intent (in)  :: k
      real*8                :: factorial
      real*8                :: res
      res = factorial (n) / (factorial (k) * factorial (n - k))
   end function choose
!========================================================================

! !========================================================================
!    ! Returns the inverse of a matrix calculated by finding the LU
!    ! decomposition.  Depends on LAPACK.
!    function inv(A) result(Ainv)
!       real, dimension(:,:), intent(in) :: A
!       real, dimension(size(A,1),size(A,2)) :: Ainv

!       real, dimension(size(A,1)) :: work  ! work array for LAPACK
!       integer, dimension(size(A,1)) :: ipiv   ! pivot indices
!       integer :: n, info

!       ! External procedures defined in LAPACK
!       external DGETRF
!       external DGETRI

!       ! Store A in Ainv to prevent it from being overwritten by LAPACK
!       Ainv = A
!       n = size(A,1)

!       ! DGETRF computes an LU factorization of a general M-by-N matrix A
!       ! using partial pivoting with row interchanges.
!       call DGETRF(n, n, Ainv, n, ipiv, info)

!       if (info /= 0) then
!          write(6,*) "Matrix is numerically singular!"
!          stop 
!       end if

!       ! DGETRI computes the inverse of a matrix using the LU factorization
!       ! computed by DGETRF.
!       call DGETRI(n, Ainv, n, ipiv, work, n, info)

!       if (info /= 0) then
!          write(6,*) "Matrix inversion failed!"
!          stop 
!       end if
!    end function inv
! !========================================================================

!========================================================================
   subroutine get_n_depth(ndepth_int,depth_int_min,depth_int_max,depth_int_d,ndepth)
      implicit none
      integer, intent(in)  :: ndepth_int
      real,    intent(in)  :: depth_int_min(ndepth_int)
      real,    intent(in)  :: depth_int_max(ndepth_int)
      real,    intent(in)  :: depth_int_d(ndepth_int)
      integer, intent(out) :: ndepth
      
      integer ii

      ndepth = 0.0

      do ii = 1,ndepth_int
         ndepth = ndepth + int((depth_int_max(ii) - depth_int_min(ii))/depth_int_d(ii))
      enddo

      ndepth = ndepth + 1

   end subroutine
!========================================================================

!========================================================================
   subroutine get_depths(ndepth,ndepth_int,depth_int_min,depth_int_max,depth_int_d,&
                         depths,thickness)
      implicit none
      integer, intent(in)  :: ndepth, ndepth_int
      real,    intent(in)  :: depth_int_min(ndepth_int)
      real,    intent(in)  :: depth_int_max(ndepth_int)
      real,    intent(in)  :: depth_int_d(ndepth_int)
      real,    intent(out) :: depths(ndepth)
      real,    intent(out) :: thickness(ndepth-1)

      integer ii, jj, count

      count = 0
      do ii = 1,ndepth_int
         do jj = 1, int((depth_int_max(ii) - depth_int_min(ii))/depth_int_d(ii))

            depths(count + jj) = depth_int_min(ii) + depth_int_d(ii)*(jj-1)
            thickness(count + jj) = depth_int_d(ii)
         enddo
         count = count + int((depth_int_max(ii) - depth_int_min(ii))/depth_int_d(ii))
      enddo
      depths(ndepth) = depth_int_max(ndepth_int)

   end subroutine
!========================================================================

!========================================================================
   subroutine get_grid_size(topo_file,nx,ny)
      implicit none
      character(LEN=*), intent(in)  :: topo_file
      integer, intent(out) :: nx, ny
      integer dummy 

      open(unit=31,file=trim(adjustl(topo_file)) )
      read(31,*) dummy
      read(31,*) nx
      read(31,*) ny
      close(31)

   end subroutine
!========================================================================

!========================================================================
   subroutine get_grid(nx,ny,topo_file,ice,water,earth,igrid)
      implicit none
      integer, intent(in)           :: nx, ny 
      character(LEN=*), intent(in)  :: topo_file
      real, intent(out)             :: ice(nx,ny), water(nx,ny), earth(nx,ny) 
      integer, intent(out)          :: igrid(nx,ny,2)
      integer dummy, ix, iy, xt, yt
      real topot, icet, watert, eartht, dummyr
      

      open(unit=31,file=trim(adjustl(topo_file)) )
      read(31,*) dummy
      read(31,*) dummy
      read(31,*) dummy
      read(31,*) dummyr

      do ix = 1,nx
         do iy = 1,ny
            read(31,*) dummyr, dummyr, topot, icet, watert, eartht, xt, yt
            igrid(ix,iy,1) = xt 
            igrid(ix,iy,2) = yt

            ice(xt,yt)   = icet
            water(xt,yt) = watert
            ! earth(xt,yt) = eartht 
            earth(xt,yt) = topot 

         enddo
      enddo

      close(31)

   end subroutine
!========================================================================

!========================================================================
   subroutine draw_new_par(par,sig,min,max,iseed,Tstar,temp_par)
      implicit none
      real par, sig, min, max, Tstar, temp_par
      integer iseed
      real Gasdev, ran1

      do 
         ! Gaussian perturbation
         temp_par = par + sig*Gasdev(iseed)*Tstar
         
         if ( (temp_par .gt. min) .and. &
              (temp_par .lt. max) ) then
            exit
         endif
      enddo 

   end subroutine draw_new_par
!========================================================================

!========================================================================
   subroutine freq_match(nfreq,nfreq_o,freqs,observed_freqs,average_disp)
      implicit none

      integer nfreq,nfreq_o
      integer ifreq, ifreq_o, ilow

      real freqs(nfreq),observed_freqs(nfreq_o),average_disp(nfreq)
      real temp_disp(nfreq)

      temp_disp = 0.0
      ilow  = 1
      do ifreq_o = 1,nfreq_o
         do ifreq = ilow,nfreq
            if ( abs(freqs(ifreq) - observed_freqs(ifreq_o)) .lt. 0.00001) then
               temp_disp(ifreq_o) = average_disp(ifreq)
               ilow  = ifreq + 1
               exit
            endif
         enddo
      enddo

      average_disp = temp_disp

   end subroutine
!========================================================================

!========================================================================
   subroutine freq_match_index(nfreq,nfreq_o,freqs,observed_freqs,freq_index)
      implicit none

      integer nfreq,nfreq_o
      integer ifreq, ifreq_o, ilow
      integer freq_index(nfreq)

      real freqs(nfreq),observed_freqs(nfreq_o)
      
      freq_index = 0
      ilow  = 1
      do ifreq_o = 1,nfreq_o
         do ifreq = ilow,nfreq
            if ( abs(freqs(ifreq) - observed_freqs(ifreq_o)) .lt. 0.00001) then
               ilow  = ifreq + 1
               freq_index(ifreq_o) = ifreq
               exit
            endif
         enddo
      enddo

      ! write(6,*) nfreq_o
      ! write(6,*) 1./observed_freqs
      ! write(6,*) nfreq
      ! write(6,*) 1./freqs
      ! write(6,*) freq_index
      ! stop

   end subroutine
!========================================================================















