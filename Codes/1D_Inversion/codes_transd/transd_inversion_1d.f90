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



program inversion_main
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
   write(6,'(a)') 'Reading data from ' // trim(adjustl(data_directory))

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


   ! Parallel Loop over filenames
   CALL OMP_SET_NUM_THREADS(bounds%nthread)
!$OMP PARALLEL DEFAULT(NONE), SHARED(data_directory,filenames,Tstar,ntemp,iseed,max_itr,bounds,nfiles,write_freq), &
!$OMP PRIVATE(ifile, thread_id, ppos, command, filename, file_exists)

   thread_id = omp_get_thread_num()

!$OMP DO

   do ifile = 1,nfiles

      ! if (ifile .lt. 10) then
      !    write(6,*) filenames(ifile)
      ! endif

      ! if ( "AK_BESE__AT_CRAG_phase.txt" .ne. &
      !    trim(adjustl(filenames(ifile))) )  then
      !    cycle 
      ! endif


      ! EXAMPLE FOR RUNNING ON SPECIFIC STATION PAIR
      ! if ( "AT_CRAG__TA_R33M_phase.txt" .ne. &
      !    trim(adjustl(filenames(ifile))) )  then
      !    cycle 
      ! endif




      ! ppos = scan(trim(filenames(ifile)),".", BACK= .true.) - 1
      ! filename = filenames(ifile)
      ! command = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_misfit.out'
      ! inquire(FILE = trim(adjustl(command)), EXIST=file_exists)
      ! if (file_exists) then
      !    write(6,*) 'Station pair already inverted: ' // trim(adjustl(filename(1:ppos)))
      !    cycle
      ! endif


      ! Meat of the MCMC code
      call transd_1d_phase_rjmcmc(thread_id,ntemp,iseed,max_itr,filenames(ifile),data_directory,Tstar,bounds,write_freq)
      

   enddo

!$OMP END DO

!$OMP END PARALLEL

end program inversion_main
!========================================================================
!========================================================================
!========================================================================
!========================================================================

!========================================================================





!========================================================================
!========================================================================
!========================================================================
   subroutine transd_1d_phase_rjmcmc(thread_id,ntemp,iseed,max_itr,filename,data_directory,Tstar,bounds,write_freq)    
      use Mod_Env
      implicit none

      integer ntemp,iseed,max_itr,thread_id

      integer nfreqs, ifreq, itemp, inode, ipar, ierr, itr

      integer ppos, write_freq

      real Tstar(ntemp), ss_dist, temp

      real ran1

      character(LEN=150) data_directory
      character(LEN=200) filename

      character(LEN=100) misfit_file, par_file, residual_file, wrapped_file

      type(Mbounds)                                      :: bounds  

      type(dispersion)                                   :: observed

      type(model_1D), dimension(:), allocatable          :: current

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
      allocate(current(ntemp))
      do itemp=1,ntemp
         allocate( current(itemp)%par(bounds%npar,bounds%max_nodes),&
                   current(itemp)%predicted(observed%nfreqs), current(itemp)%vels(observed%nfreqs) )

      enddo


      ! ! Initialize Geopsy forward model (for dispersion curves) 
      ! call dispersion_curve_init(0)


      !=========== Random Starting Model ==============
      do itemp=1,ntemp

         do    
               
            current(itemp)%par = 0.0
            current(itemp)%nodes = 0
            current(itemp)%predicted = 0.0

            current(itemp)%attempt = 0
            current(itemp)%success = 0


            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            current(itemp)%nodes = ceiling(ran1(iseed))*10

            do inode = 1, current(itemp)%nodes

               do ipar = 1, bounds%npar

                  current(itemp)%par(ipar,inode) = ran1(iseed)*( bounds%max_par(ipar)-bounds%min_par(ipar) ) +  &
                                                   bounds%min_par(ipar)

               enddo

            enddo
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! current(itemp)%nodes = 3

            ! current(itemp)%par(1,1) = 10.
            ! current(itemp)%par(1,2) = 40.
            ! current(itemp)%par(1,3) = 60.

            ! current(itemp)%par(2,1) = 3.0
            ! current(itemp)%par(2,2) = 3.5
            ! current(itemp)%par(2,3) = 4.5

            ierr = 0
            call forward_calc_1d(current(itemp),bounds,observed,ierr)

            if ((ierr == 0)) then
               ! write(6,*) current(itemp)%misfit
               exit
            endif

         enddo
      enddo


      !=========== MHS =================================

      ppos = scan(trim(filename),".", BACK= .true.) - 1

      misfit_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_misfit.out'
      par_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_par.out'
      residual_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_residual.out'

      wrapped_file = 'results_1d/' // trim(adjustl(filename(1:ppos))) // '_wrapped.out'


      open(unit=90 + thread_id,file=trim(misfit_file),status='replace',form='formatted') 
      open(unit=60 + thread_id,file=trim(par_file),status='replace',form='unformatted')
      open(unit=30 + thread_id,file=trim(residual_file),status='replace',form='unformatted')

      open(unit=8 + thread_id,file=trim(wrapped_file),status='replace',form='formatted') 


      itr = 1
      do while (itr .lt. max_itr)


         do itemp=1,ntemp 

            if (mod(itr,5) .eq. 0) then
                        
               select case(ceiling(ran1(iseed)*2))

                  case(1)
                     call birth_1d(current(itemp),bounds,observed,Tstar(itemp),iseed)   

                  case(2)
                     call death_1d(current(itemp),bounds,observed,Tstar(itemp),iseed)

               end select

            else

               call perturb_one_1d(current(itemp),bounds,observed,Tstar(itemp),iseed)  

            endif

         enddo


         ! ! Write update to screen
         ! if (mod(itr,100) == 0) then
         !    write(6,*) " "
         !    write(6,*) "Iteration  |   Temp  |    Misfit   |  Nodes  |",&
         !    " Perturb | Birth |  Death  |   Swap   |   Fail"
         !    do itemp=1,ntemp 

         !       write(6,"(I10,F10.3,F13.3,I10,5F10.3)") &
         !                                 itr, &
         !                                 Tstar(itemp), &
         !                                 current(itemp)%misfit, &
         !                                 current(itemp)%nodes, &
         !                                (current(itemp)%success(1)/current(itemp)%attempt(1)), &
         !                                (current(itemp)%success(2)/current(itemp)%attempt(2)), &
         !                                (current(itemp)%success(3)/current(itemp)%attempt(3)), &
         !                                (current(itemp)%success(6)/current(itemp)%attempt(6)),&
         !                                (current(itemp)%success(5)/current(itemp)%attempt(5))

         !    enddo
         !    ! write(6,*) 2.0, 2.5, 3.0, 3.2, 3.3, 3.4, 3.5, 4.5, 4.5, 4.5, 5.5, 40.0
         !    ! write(6,*) current%par

         ! endif


         ! Write model to file
         if (mod(itr,write_freq) == 0) then
            ! write(6,*) itr, ' / ', max_itr

            write(90 + thread_id,*) current(1)%misfit, current(1)%nodes
            write(60 + thread_id) current(1)%par 
            write(30 + thread_id) current(1)%vels

            write(8 + thread_id,*) current(1)%predicted

            ! write(6,*) 1./observed%freqs(9), current(1)%vels(9)

         endif

         ! ! Attempt t-swap
         ! if (mod(itr,20) == 0) then
         !    call tswap(current,bounds,ntemp,iseed,Tstar)
         ! endif


         itr = itr + 1
      enddo



      !=========== Terminate Program =================== 
      close(90 + thread_id)
      close(60 + thread_id)
      close(30 + thread_id)
      close(8 + thread_id)

      write(6,*) "completed Inversion for ", filename(1:ppos)


   end subroutine transd_1d_phase_rjmcmc
!========================================================================
!========================================================================
!========================================================================






!========================================================================
!========================================================================
!========================================================================
!========================================================================
   subroutine forward_calc_1d(current,bounds,observed,errout)
      use Mod_Env
      implicit none
      integer errout
      type(model_1D)                      :: current
      type(Mbounds)                       :: bounds
      type(dispersion)                    :: observed

      integer ii ,ix,iy,isamp,itest,index
      real lppd
      integer*8 count1, count2, count_rate, count_max

      errout = 0
      call predict_disp(current,bounds,observed,errout)

      if (errout .eq. 1) then
         lppd = huge(lppd)
         current%misfit = lppd
         return
      endif

      ! L2 misfit
      ! lppd = 0.0     
      ! do ii=1,bounds%nfreq(1)

      !    lppd = lppd + ((observed(index)%vels(ii) - current%predicted(ii) )**2.) /&
      !                   (2.*(observed(index)%err(ii)**2.) )

      ! enddo

         ! write(6,*) ''
         ! write(6,*) ''

      ! Von Mises Misfit 
      lppd = 0.0     
      do ii=1,observed%nfreqs

         lppd = lppd - observed%kappas(ii)*cos(current%predicted(ii) - observed%phases(ii))

         ! write(6,*) 1./observed%freqs(ii), current%predicted(ii), observed%phases(ii)
         ! write(6,*) ''
      enddo


      current%misfit = lppd
         
   end subroutine forward_calc_1d
!========================================================================
!========================================================================
!========================================================================


!========================================================================
!========================================================================
!========================================================================
   subroutine predict_disp(current,bounds,observed,errout)
      use Mod_Env
      implicit none
      integer errout
 
      type(model_1D)                      :: current
      type(Mbounds)                       :: bounds  
      type(dispersion)                    :: observed

      !=================================================================

      real h(bounds%max_nodes+10), vs(bounds%max_nodes+10)
      real vp(bounds%max_nodes+10), rho(bounds%max_nodes+10)

      integer nlay, ifreq
      real vels(observed%nfreqs), slow(observed%nfreqs)
      ! integer index

      integer inode !, inode2, imin, iprev
      real depth_prev, temp

      ! integer nfreq
      ! real freqs(3), vels_temp(3), distance


      logical mk(current%nodes)
      integer isort


      ! real ak135_hh(7), ak135_vs(7), ak135_vp(7), ak135_rho(7)

      !=================================================================
      interface   
        subroutine SWD_calc_CPS_rayleigh(nlay,nfreq,&
                            H,Vp,Vs,rho,&
                            freqs,vels,&
                            nmodes,group,ierr)
           integer :: nlay,nfreq,nmodes,group,ierr
           real    :: H(nlay),Vp(nlay),Vs(nlay),rho(nlay)
           real    :: freqs(nfreq),vels(nfreq)
        end subroutine 
      end interface
      !=================================================================

      !=================================================================
      interface   
        subroutine SWD_calc(nlay,nfreq,&
                            H,Vp,Vs,rho,&
                            freqs,vels,&
                            nmodes,group,ierr)
           integer :: nlay,nfreq,nmodes,group,ierr
           real    :: H(nlay),Vp(nlay),Vs(nlay),rho(nlay)
           real    :: freqs(nfreq),vels(nfreq)
        end subroutine 
      end interface
      !=================================================================

      ! ak135_hh  = (/ 20., 45., 45., 50., 50., 50., 50. /) 
      ! ak135_rho = (/ 3.4268, 3.3711, 3.3243, 3.3663, 3.4110, 3.4577, 3.5068 /) 
      ! ak135_vp  = (/ 8.0505, 8.1750, 8.3007, 8.4822, 8.6650, 8.8476, 9.0302 /) 
      ! ak135_vs  = (/ 4.5000, 4.5090, 4.5184, 4.6094, 4.6964, 4.7832, 4.8702 /) 






      mk = .true.

      errout = 0


      ! write(6,*) ' '
      ! write(6,*) ' '
      ! write(6,*) current%nodes
      ! write(6,*) ' '

      nlay = 0
      depth_prev = 0.0
      ! Find node depth sort order
      do inode = 1, current%nodes
         isort = MINLOC(current%par(1,1:current%nodes),1,mk)
         mk(isort) = .false.
 
         ! if (inode .eq. 1) then 
         !    vs_min = current%par(2,isort)
         ! endif

         nlay      = nlay + 1
         h(nlay)   = current%par(1,isort) - depth_prev
         ! vs(nlay)  = current%par(2,isort)
         call get_ak135_value(current%par(1,isort),current%par(2,isort),temp)
         vs(nlay)  = temp

         vp(nlay)  = vs(nlay)*1.75
         rho(nlay) = 2.35 + 0.036*((vp(nlay)-3.0)**2.0)

         depth_prev = current%par(1,isort)

         ! write(6,*) nlay, h(nlay), vs(nlay), vp(nlay), rho(nlay)

      enddo


      ! stop

      ! if (nlay .ne. maxloc(vs(1:nlay),1) ) then
      !    errout = 1 
      !    return
      ! endif

      ! if (1 .ne. minloc(vs(1:nlay),1) ) then
      !    errout = 1 
      !    return
      ! endif

      ! nlay      = nlay + 1
      ! h(nlay)   = bounds%max_par(1) - depth_prev
      ! vs(nlay)  = ( current%par(2,isort) + ak135_vs(1) )/2.
      ! vp(nlay)  = vs(nlay)*1.75
      ! rho(nlay) = 2.35 + 0.036*((vp(nlay)-3.0)**2.0)

      ! write(6,*) nlay, h(nlay), vs(nlay), vp(nlay), rho(nlay)
      ! write(6,*) ' '



      ! ! Add extra layers based on AK135
      ! do inode = 1,7

      !    nlay      = nlay + 1
      !    h(nlay)   = ak135_hh(inode)
      !    vs(nlay)  = ak135_vs(inode)
      !    vp(nlay)  = ak135_vp(inode)
      !    rho(nlay) = ak135_rho(inode)

      ! enddo



      ! Try to get more accurate derivatives for group velocity
      if (bounds%group .eq. 1) then

         write(6,*) 'Fix code for group velocity'
 
         ! nfreq = 3

         ! do ifreq = 1, bounds%nfreq(idataset)

         !    freqs(1) = 1./((1./bounds%freqs(idataset,ifreq)) + 1.0)  
         !    freqs(2) = bounds%freqs(idataset,ifreq)  
         !    freqs(3) = 1./((1./bounds%freqs(idataset,ifreq)) - 1.0)  

         !    !===========================================
         !    call SWD_calc_CPS_rayleigh(nlay,nfreq,&
         !                         h(1:nlay),&
         !                         vp(1:nlay),&
         !                         vs(1:nlay),&
         !                         rho(1:nlay),&
         !                         freqs,&
         !                         vels_temp,&
         !                         1,&
         !                         1,errout)
         !    !===========================================

         !    if (errout .eq. 1) then
         !       return
         !    else
         !       vels(ifreq) = vels_temp(2)
         !       ! write(6,*) vels_temp(2)
         !    endif

         ! enddo

      
         ! ! stop


      else

         ! write(6,*) nlay

         !===========================================
         call SWD_calc(nlay,observed%nfreqs,&
                              h(1:nlay),&
                              vp(1:nlay),&
                              vs(1:nlay),&
                              rho(1:nlay),&
                              observed%freqs,&
                              vels,&
                              1,&
                              0,errout)
         !===========================================
         
         ! write(6,*) 'test 2', errout

         if (errout .eq. 1) return


         ! call dispersion_curve_rayleigh(nlay,&
         !                                h,&
         !                                vp,&
         !                                vs,&
         !                                rho,&
         !                                observed%nfreqs,&
         !                                2.*3.14159*observed%freqs,&
         !                                1,slow,0)
         ! vels = 1./slow

      endif

      
      ! Unwrapped dispersion
      current%vels = vels

      call wrap_phase(observed%nfreqs,observed%freqs,observed%ss_dist,vels)
      ! stop
      ! Wrapped phase
      current%predicted = vels

      ! write(6,*)
      ! do ifreq = 1, observed%nfreqs
      !    write(6,*) 1./observed%freqs(ifreq), vels(ifreq)
      ! enddo
      ! write(6,*)

      ! stop

   end subroutine predict_disp
!========================================================================
!========================================================================
!========================================================================



!========================================================================
!========================================================================
!========================================================================
   subroutine wrap_phase(nfreqs,freqs,ss_dist,vels)
      use Mod_Env
      implicit none

      integer nfreqs, ifreq, revolutions

      real freqs(nfreqs), ss_dist, vels(nfreqs), temp

      !=================================================================

      do ifreq = 1,nfreqs

         ! write(6,*) ss_dist, freqs(ifreq), vels(ifreq)

         temp = -1.*ss_dist*2.*3.14159*freqs(ifreq)/vels(ifreq)


         ! write(6,*) ''
         ! write(6,*) temp

         if (temp .gt. 0.0) then

            revolutions = 0 
            do while (temp > 3.14159)
               temp = temp - 2.*3.14159
               revolutions = revolutions + 1
               ! write(6,*) temp

            enddo

         else

            revolutions = 0 
            do while (temp < -3.14159)
               temp = temp + 2.*3.14159
               revolutions = revolutions - 1
            enddo

         endif

         ! write(6,*) temp
         ! stop

         vels(ifreq) = temp

      enddo

   end subroutine wrap_phase
!========================================================================
!========================================================================
!========================================================================


!========================================================================
!========================================================================
!========================================================================
   subroutine perturb_one_1d(current,bounds,observed,Tstar,iseed)      
      use Mod_Env
      implicit none

      real Tstar, Gasdev,ran1, ratio, temp_par

      integer iseed, inode, ipar, errout

      type(model_1D)                          :: current
      type(model_1D)                          :: proposed1

      type(Mbounds)                           :: bounds
      type(dispersion)                        :: observed




      current%attempt(1) = current%attempt(1) + 1
      
      inode = ceiling(ran1(iseed)*current%nodes)

      ! Set proposed model to current model 
      proposed1 = current 


      ! Perturb node
      ipar = ceiling(ran1(iseed)*bounds%npar)
      
      call draw_new_par(current%par(ipar,inode),&
                        bounds%sig_par(ipar),&
                        bounds%min_par(ipar),&
                        bounds%max_par(ipar),iseed,1.0,&
                        temp_par)

      proposed1%par(ipar,inode) = temp_par


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Recompute predicted data
      errout = 0
      call forward_calc_1d(proposed1,bounds,observed,errout)
      current%attempt(5) = current%attempt(5) + 1       
      if (errout .eq. 1) then
         current%success(5) = current%success(5) + 1
         return
      endif  

      ratio = exp( (current%misfit - proposed1%misfit)*(1./Tstar) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ran1(iseed) .lt. ratio) then
         current%par                  = proposed1%par
         current%misfit               = proposed1%misfit
         current%predicted            = proposed1%predicted
         current%vels                 = proposed1%vels

         current%success(1)           = current%success(1) + 1
      endif

   end subroutine
!========================================================================

!========================================================================
   subroutine birth_1d(current,bounds,observed,Tstar,iseed)      
      use Mod_Env
      implicit none


      real Tstar, Gasdev,ran1, ratio

      integer iseed, inode, ipar, errout

      type(model_1D)                          :: current
      type(model_1D)                          :: proposed1

      type(Mbounds)                           :: bounds
      type(dispersion)                        :: observed



      current%attempt(2) = current%attempt(2) + 1

      !!! Not the proper way to return... fix this
      if (current%nodes .eq. bounds%max_nodes) return 

      ! Set proposed model to current model 
      proposed1 = current 

      proposed1%nodes = proposed1%nodes + 1

      ! Draw parameter value and new node location
      do ipar = 1, bounds%npar
         proposed1%par(ipar,proposed1%nodes) = (bounds%max_par(ipar) - bounds%min_par(ipar))*ran1(iseed) + bounds%min_par(ipar)
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Recompute predicted data
      errout = 0
      call forward_calc_1d(proposed1,bounds,observed,errout)
      current%attempt(5) = current%attempt(5) + 1       
      if (errout .eq. 1) then
         current%success(5) = current%success(5) + 1
         return
      endif  

      ratio = exp( (current%misfit - proposed1%misfit)*(1./Tstar) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Poisson prior on nodes
      ! ratio = ratio*real(current%nodes)/real(proposed1%nodes)
      ratio = ratio*real(bounds%min_nodes+2)/real(proposed1%nodes)


      if (ran1(iseed) .lt. ratio) then

         current%par                  = proposed1%par
         current%misfit               = proposed1%misfit
         current%predicted            = proposed1%predicted
         current%vels                 = proposed1%vels

         current%nodes                = proposed1%nodes

         current%success(2)           = current%success(2) + 1

      endif

   end subroutine
!========================================================================

!========================================================================
   subroutine death_1d(current,bounds,observed,Tstar,iseed)      
      use Mod_Env
      implicit none


      real Tstar, Gasdev,ran1, ratio

      integer iseed, inode, ipar, errout

      type(model_1D)                          :: current
      type(model_1D)                          :: proposed1

      type(Mbounds)                           :: bounds
      type(dispersion)                        :: observed

      !!! Not the proper way to return... fix this
      if (current%nodes .eq. bounds%min_nodes) return 

      current%attempt(3) = current%attempt(3) + 1

      ! Set proposed model to current model 
      proposed1 = current 

      proposed1%nodes = proposed1%nodes - 1

      ! Choose random node
      inode = ceiling(ran1(iseed)*current%nodes)

      do ipar = 1,bounds%npar
         proposed1%par(ipar,inode) = proposed1%par(ipar,current%nodes)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Recompute predicted data
      errout = 0
      call forward_calc_1d(proposed1,bounds,observed,errout)
      current%attempt(5) = current%attempt(5) + 1       
      if (errout .eq. 1) then
         current%success(5) = current%success(5) + 1
         return
      endif  

      ratio = exp( (current%misfit - proposed1%misfit)*(1./Tstar) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Poisson prior on nodes
      ! ratio = ratio*real(current%nodes)/real(proposed1%nodes)
      ratio = ratio*real(current%nodes)/real(bounds%min_nodes+2)

      if (ran1(iseed) .lt. ratio) then
         
         current%par                  = proposed1%par
         current%misfit               = proposed1%misfit
         current%predicted            = proposed1%predicted
         current%vels                 = proposed1%vels

         current%nodes                = proposed1%nodes

         current%success(3)           = current%success(3) + 1

      endif

   end subroutine
!========================================================================




!=======================================================================
! Probabilistic temperature swapping between parallel sampling chains
!=======================================================================
   subroutine tswap(current,bounds,nT,iseed,Tstar)
      use Mod_Env
      implicit none

      integer ii, first, second, group
      integer iseed, nT
      real ran1, ratio
      real Tstar(nT)

      type(model_1D), dimension(nT)           :: current
      type(model_1D)                          :: holder
      type(Mbounds)                           :: bounds

      first  = ceiling(ran1(iseed)*nT) 
      second = ceiling(ran1(iseed)*nT) 
	      
      if (first .ne. second) then
            
            ratio = exp( (current(first)%misfit - current(second)%misfit)* &
                          ( (1./Tstar(first)) - (1./Tstar(second)) ) )
      
            if (ran1(iseed) .lt. ratio) then
                 holder          = current(first)
                 current(first)  = current(second)
                 current(second) = holder

                 holder%attempt          = current(first)%attempt
                 holder%success          = current(first)%success
                 current(first)%attempt  = current(second)%attempt
                 current(first)%success  = current(second)%success
                 current(second)%attempt = holder%attempt
                 current(second)%success = holder%success

                 current(first)%success(6)  = current(first)%success(6) + 1
                 current(second)%success(6) = current(second)%success(6) + 1
            endif
            current(first)%attempt(6) = current(first)%attempt(6) + 1
            current(second)%attempt(6) = current(second)%attempt(6) + 1
      endif

   end subroutine    
!========================================================================
!========================================================================
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

      ! Define a custom (wide) prior bound that is depth dependent 
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
