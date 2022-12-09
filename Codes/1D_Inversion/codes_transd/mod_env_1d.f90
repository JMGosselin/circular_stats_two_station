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



!=========================================================================
!module containing fortran derived data-types describing 
!(Mod) model and parameter (Env) environment variables
module Mod_Env 
      implicit none 

      type :: model_1D 

            real                                   :: misfit

            integer                                :: nodes 
            real, dimension(:,:), allocatable      :: par

            real*8, dimension(10)                  :: attempt
            real*8, dimension(10)                  :: success
            
            real, dimension(:), allocatable        :: predicted
            real, dimension(:), allocatable        :: vels

      end type model_1D

      type :: dispersion

            integer                                :: nfreqs

            real                                   :: ss_dist

            real,  dimension(:),     allocatable   :: freqs
            real,  dimension(:),     allocatable   :: phases
            real,  dimension(:),     allocatable   :: kappas

      end type dispersion 
 

      type :: Mbounds
            integer                              :: npar  
            
            integer                              :: nthread
            
            ! For SWD calculation
            integer                              :: geopsy
            integer                              :: nmode
            integer                              :: group
            
            ! Model parameter priors
            integer                              :: min_nodes
            integer                              :: max_nodes
            real, dimension(:),   allocatable    :: min_par
            real, dimension(:),   allocatable    :: max_par
            real, dimension(:),   allocatable    :: sig_par

      end type Mbounds   

end module Mod_Env
!=========================================================================
