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



program test
   !use mpi
   use Mod_Env
   use omp_lib
   !use mpi
   implicit none


!=========== Define Variables ===================
   integer max_itr, ntemp, iseed, nfiles, ifile
   integer thread_id, i, ii, reason, ppos, cutoff, write_freq

   real temp_fact, r

   real, dimension(:), allocatable                    :: Tstar

   type(Mbounds)                                      :: bounds  

   character(LEN=150)                                 :: data_directory, filename
   character(LEN=200), dimension(:), allocatable      :: filenames
   character(LEN=200)                                 :: command
   
   logical :: file_exists




   open(31,file='inversion_parameters.in',action="read")

!=========== User Input ===================

   read(31,*) max_itr
   read(31,*) write_freq 
   read(31,*) cutoff  

   ! PT stuff, currently not used
   ntemp = 1
   temp_fact = 1.2
   allocate(Tstar(ntemp))
   do ii=1,ntemp
      Tstar(ii) = temp_fact**(ii-1)
   enddo

   read(31,*) iseed
   read(31,*) bounds%min_nodes
   read(31,*) bounds%max_nodes
   read(31,*) bounds%npar

   allocate(bounds%min_par(bounds%npar),&
            bounds%max_par(bounds%npar),&
            bounds%sig_par(bounds%npar))

   read(31,*) bounds%min_par(1)
   read(31,*) bounds%max_par(1) 
   read(31,*) bounds%sig_par(1)

   read(31,*) bounds%min_par(2)
   read(31,*) bounds%max_par(2) 
   read(31,*) bounds%sig_par(2)

   read(31,*) bounds%nthread

   read(31,'(a)') data_directory

   ! Stuff for CPS
   bounds%nmode = 1
   bounds%group = 0                            

!=========== User Input ===================

   close(31)






   ! Get data filenames
   write(6,'(a)') 'Reading data from ' // data_directory

   inquire(FILE=trim(adjustl(data_directory)) // "contents.txt", EXIST=file_exists)

   if (.not. file_exists) then
      command = "ls "  // trim(adjustl(data_directory)) // " > contents.txt"
      call system(trim(adjustl(command)))
   endif

   open(31,file='contents.txt',action="read")
   i = 0
   do
      read(31,FMT='(a)',iostat=reason) r
      if (reason/=0) EXIT
      i = i+1
   end do
   nfiles = i

   write(6,'(a,I0)') "Number of data files: " , nfiles

   allocate(filenames(nfiles))
   rewind(31)
   do ifile = 1,nfiles
      read(31,'(a)') filenames(ifile)
   enddo
   close(31)
   !!!!!!!!!!!!!!!!!!!!!


!    ! Parallel Loop over filenames
!    CALL OMP_SET_NUM_THREADS(bounds%nthread)
! !$OMP PARALLEL DEFAULT(NONE), SHARED(data_directory,filenames,Tstar,ntemp,iseed,max_itr,bounds,nfiles), &
! !$OMP PRIVATE(ifile, thread_id)

   ! thread_id = omp_get_thread_num()
   thread_id = 0

! $OMP DO

   do ifile = 1,nfiles

      ! if (ifile .lt. 10) then
      !    write(6,*) filenames(ifile)
      ! endif

      ! if ( "AK_BESE__AT_CRAG_phase.txt" .ne. &
      !    trim(adjustl(filenames(ifile))) )  then
      !    cycle 
      ! endif

      
      ! Example if you want to test on one station pair
      ! if ( "AT_CRAG__TA_R33M_phase.txt" .ne. &
      !    trim(adjustl(filenames(ifile))) )  then
      !    cycle 
      ! endif

      ! ppos = scan(trim(filenames(ifile)),".", BACK= .true.) - 1
      ! filename = filenames(ifile)
      ! command = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_misfit.out'
      ! inquire(FILE = trim(adjustl(command)), EXIST=file_exists)
      ! if (.not. file_exists) then
      !    cycle
      ! endif

      ! Meat of the MCMC code
      call transd_1d_phase_rjmcmc_grid(thread_id,ntemp,iseed,max_itr,filenames(ifile),data_directory,Tstar,bounds,cutoff)
      

   enddo

! $OMP END DO

! $OMP END PARALLEL

end program test
!========================================================================
!========================================================================
!========================================================================
!========================================================================

!========================================================================




!========================================================================
!========================================================================
!========================================================================
   subroutine transd_1d_phase_rjmcmc_grid(thread_id,ntemp,iseed,max_itr,filename,data_directory,Tstar,bounds,cutoff)    
      use Mod_Env
      implicit none

      integer ntemp,iseed,max_itr,thread_id,reason,ii

      integer nfreqs, ifreq, itemp, inode, ipar, ierr, itr

      integer ppos, cutoff

      real Tstar(ntemp), ss_dist, temp

      real ran1

      character(LEN=150) data_directory
      character(LEN=200) filename

      character(LEN=100) misfit_file, par_file, residual_file, model_file

      type(Mbounds)                                      :: bounds  

      type(dispersion)                                   :: observed

      type(model_1D), dimension(:), allocatable          :: current

      integer ndepth, nvs, ndisp

      integer, dimension(:,:), allocatable               :: profile, disp_grid
      real, dimension(:), allocatable                    :: depth_axis, vs_axis
      real, dimension(:), allocatable                    :: dispersion_axis



      integer nfreq_temp
      real,  dimension(:),     allocatable               :: freqs
      real,  dimension(:),     allocatable               :: phases
      real,  dimension(:),     allocatable               :: kappas

      integer ikeep_periods(10), ik

      ! ikeep_periods = (/ 10, 12, 15, 20, 24, 30, 40, 48, 60, 80 /) 

      ! Read phase observed spectrum parameters
      open(unit=30 + thread_id, file=trim(adjustl(data_directory)) // trim(adjustl(filename)) )
      ! read(31,*) observed%nfreqs
      read(30 + thread_id,*) nfreq_temp
      read(30 + thread_id,*) observed%ss_dist

      allocate( freqs(nfreq_temp), phases(nfreq_temp), kappas(nfreq_temp) )

      do ifreq = 1, nfreq_temp
         read(30 + thread_id,*) freqs(ifreq), phases(ifreq), kappas(ifreq)
         ! write(6,*) freqs(ifreq), phases(ifreq), kappas(ifreq)
      enddo
      close(30 + thread_id)


      observed%nfreqs = nfreq_temp

      allocate( observed%freqs(observed%nfreqs), &
                observed%kappas(observed%nfreqs), &
                observed%phases(observed%nfreqs) )

      ! write(6,*) ''

      do ik = 1, observed%nfreqs

         observed%freqs(ik) = 1./freqs(ik) ! Period to frequency
         observed%kappas(ik) = kappas(ik)  ! Can apply scaling factor here
         observed%phases(ik) = phases(ik)

         ! observed%freqs(ik) = 1./freqs(ikeep_periods(ik)) ! Period to frequency
         ! observed%kappas(ik) = 2.*kappas(ikeep_periods(ik))
         ! observed%phases(ik) = phases(ikeep_periods(ik))
         ! write(6,*) freqs(ikeep_periods(ik)), phases(ikeep_periods(ik)), kappas(ikeep_periods(ik))
      enddo
      
      ! stop












      ! Allocate model 
      ntemp = 1
      allocate(current(ntemp))
      do itemp=1,ntemp
         allocate( current(itemp)%par(bounds%npar,bounds%max_nodes),&
                   current(itemp)%predicted(observed%nfreqs), current(itemp)%vels(observed%nfreqs) )

      enddo


      ndepth = 400
      nvs = 100

      allocate(profile(ndepth,nvs),depth_axis(ndepth),vs_axis(nvs))
      profile = 0
      call linspace(bounds%min_par(1),bounds%max_par(1),ndepth,depth_axis)
      call linspace(1.5,5.5,nvs,vs_axis)

      ndisp = 100

      allocate(disp_grid(observed%nfreqs,ndisp),dispersion_axis(ndisp))
      disp_grid = 0
      call linspace(3.0,5.0,ndisp,dispersion_axis)



      !=========== Grid =================================

      ppos = scan(trim(filename),".", BACK= .true.) - 1

      misfit_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_misfit.out'
      par_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_par.out'
      residual_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_residual.out'


      open(unit=90 + thread_id,file=trim(misfit_file),form='formatted') 
      open(unit=60 + thread_id,file=trim(par_file),form='unformatted')
      open(unit=30 + thread_id,file=trim(residual_file),form='unformatted')


      itr = 0
      do 

         read(90 + thread_id,*,iostat=reason) current(1)%misfit, current(1)%nodes
         if (reason < 0) exit
         
         read(60 + thread_id,iostat=reason) current(1)%par 
         if (reason < 0) exit

         read(30 + thread_id,iostat=reason) current(1)%vels
         if (reason < 0) exit

         
         !! Specify cutoff for burn-in ...
         if (itr .gt. cutoff) then
            call profile_grid(current(1),bounds,ndepth,nvs,profile,depth_axis,vs_axis)  

            call dispersion_grid(current(1),bounds,observed,ndisp,disp_grid,dispersion_axis)  

         endif      

         itr = itr + 1

      enddo

      !=========== Terminate Program =================== 
      close(90 + thread_id)
      close(60 + thread_id)
      close(30 + thread_id)

      write(6,*) "completed Inversion for ", filename(1:ppos)


      model_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_profile_grid.out'
      open(unit=77,file=trim(model_file),status='replace',form='formatted') 
      do ii = 1,ndepth
         write(77,*) profile(ii,:)
      enddo
      close(77)

      model_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_dispersion_grid.out'
      open(unit=77,file=trim(model_file),status='replace',form='formatted') 
      do ii = 1,observed%nfreqs
         write(77,*) disp_grid(ii,:)
      enddo      
      close(77)

   end subroutine transd_1d_phase_rjmcmc_grid
!========================================================================
!========================================================================
!========================================================================



!========================================================================
   subroutine profile_grid(current,bounds,ndepth,nvs,profile,depth_axis,vs_axis)
      use Mod_Env
      implicit none
 
      type(model_1D)                      :: current
      type(Mbounds)                       :: bounds  

      !=================================================================

      integer ndepth, nvs, idepth, inode, imin, ivs

      integer profile(ndepth,nvs)
      real depth_axis(ndepth), vs_axis(nvs)
      real temp_diff, min_distance, depth, depth_diff, distance

      real ak135_vs

      depth_diff = depth_axis(2) - depth_axis(1)

      do idepth = 1,ndepth

         depth = (idepth-1)*depth_diff + depth_diff/2.

         min_distance = bounds%max_par(1)**2
         do inode = 1, current%nodes

            if (current%par(1,inode) .gt. depth) then
               distance = current%par(1,inode) - depth
               if (distance .lt. min_distance) then
                  imin = inode
                  min_distance = distance
               endif
            endif
         enddo

         call get_ak135_value(current%par(1,imin),current%par(2,imin),ak135_vs)

         ivs = MINLOC( ABS( vs_axis - ak135_vs ), 1)

         profile(idepth,ivs) = profile(idepth,ivs) + 1

      enddo

   end subroutine profile_grid
!========================================================================



!========================================================================
   subroutine dispersion_grid(current,bounds,observed,ndisp,disp_grid,dispersion_axis) 
      use Mod_Env
      implicit none
 
      type(model_1D)                      :: current
      type(Mbounds)                       :: bounds  
      type(dispersion)                    :: observed

      !=================================================================

      integer ndisp, ifreq, imin

      integer disp_grid(observed%nfreqs,ndisp)
      real dispersion_axis(ndisp)

      ! write(6,*) current%vels(1:5)

      do ifreq = 1,observed%nfreqs

         imin = MINLOC( ABS( dispersion_axis - current%vels(ifreq) ), 1)

         disp_grid(ifreq,imin) = disp_grid(ifreq,imin) + 1

      enddo

   end subroutine dispersion_grid
!========================================================================

!========================================================================
!========================================================================
!========================================================================
   subroutine get_ak135_value(depth,perturb,value)
      implicit none
      real depth, perturb, value

      real ak135_z(13), ak135_vs(13), ak135_shift(13)

      integer imin

      real temp_vs, temp_shift

      ! ak135_z     = (/ 10.,20.,35.,77.,100. /)
      ! ak135_vs    = (/ 3.0,3.46,3.85,4.49,4.50 /)
      ! ak135_shift = (/ 0.5,0.4,0.35,0.25,0.2 /)

      ! imin = minloc( abs(ak135_z - depth), 1)
      ! value = perturb*ak135_shift(imin)*ak135_vs(imin) + ak135_vs(imin)


      ak135_z     = (/ 0.,  19.9, 20.1, 34.9, 35.1, 77.0, 120., 165., 210., 260., 310., 360., 410. /)
      ak135_vs    = (/ 3.0, 3.46, 3.85, 3.85, 4.48, 4.49, 4.50, 4.51, 4.52, 4.61, 4.70, 4.78, 4.87 /)
      ak135_shift = (/ 0.5, 0.4,  0.35, 0.25, 0.2,  0.2,  0.2,  0.2,  0.15,  0.15,  0.15,  0.1,  0.1  /)

      call lin_interp1(13, ak135_z, ak135_vs, depth, temp_vs)
      call lin_interp1(13, ak135_z, ak135_shift, depth, temp_shift)

      value = perturb*temp_shift*temp_vs + temp_vs

   end subroutine get_ak135_value
!========================================================================
!========================================================================
!========================================================================

!========================================================================
!========================================================================
!========================================================================
   subroutine lin_interp1(ndat, xdat, ydat, xval, yval)
      implicit none

      integer ndat

      real xdat(ndat), ydat(ndat)
      real xval, yVal, slope, intercept

      integer ilay

      ilay = 1
      do while (xdat(ilay) .lt. xval)
         ilay = ilay + 1
      enddo

      slope = (ydat(ilay) - ydat(ilay-1)) / (xdat(ilay) - xdat(ilay-1)) 
      intercept = ydat(ilay) - slope*xdat(ilay)

      yval = slope*xval + intercept

   end subroutine lin_interp1
!========================================================================
!========================================================================
!========================================================================