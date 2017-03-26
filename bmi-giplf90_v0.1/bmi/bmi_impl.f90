      module bmif
      use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
      
      implicit none
      
      integer, parameter :: BMI_VAR_TYPE_DOUBLE = 8

      integer, parameter :: BMI_GRID_TYPE_UNKNOWN = 0
      integer, parameter :: BMI_GRID_TYPE_UNIFORM = 1
      integer, parameter :: BMI_GRID_TYPE_RECTILINEAR = 2
      integer, parameter :: BMI_GRID_TYPE_STRUCTURED = 3
      integer, parameter :: BMI_GRID_TYPE_UNSTRUCTURED = 4
      integer, parameter :: BMI_GRID_TYPE_NUMBER = 5

      type :: BMI_Model
        private
        
        ! Parameter read from cmd file
		integer :: restart                                     ! 0/1 start from previous time step / start from the begining
		real*8 :: time_step                                    ! step is the timestep in the example it is 1 yr
		real*8 :: TAUM                                         ! taum is the convergence parameter used by the stefan subroutine
		real*8 :: TMIN                                         ! tmin minimal timestep used in the Stefan subroutine
		real*8 :: time_beg,time_end                            ! inbegin time, end time
		integer :: itmax                                       ! maximum number of iterations in Stefan subroutine
		integer :: n_time                                      ! number of time steps that temp will be averaged over
		integer :: n_frz_max                                   ! maximum number of freezing fronts
		real*8 :: smooth_coef                                  ! smoothing factor
		real*8 :: unf_water_coef                               ! unfrozen water coefficient
		real*8 :: n_sec_day                                    ! number of second in a day
		real*8 :: frz_frn_max,frz_frn_min                      ! freezing front min and max depth [meters]
		real*8 :: sat_coef                                     ! saturation coefficient [dimensionless, fraction of 1]
		real*8 :: sea_level                                    ! how many meter above the sea level the borehole is
  
		character*64 :: file_sites,file_bound,file_snow,file_rsnow,file_init
		character*64 :: file_grid,file_organic,file_mineral
		character*64 :: restart_file,result_file,aver_res_file
	
		integer :: n_site                                       ! number of sites
		integer, allocatable :: snow_code(:),veg_code(:)        ! (not necccessary) required for runing in parallel
		integer, allocatable :: geo_code(:),gt_zone_code(:)     ! (not necccessary) required for runing in parallel
		real*8, allocatable  :: temp_grd(:)                     ! temprature gradient at the lower boundary

		integer :: n_temp                                       ! number of upper boundary points for temperature (input)
		real*8,allocatable::  utemp_time(:), utemp(:,:)         ! upper boundary time and temperature (input)
		real*8,allocatable::  utemp_time_i(:), utemp_i(:,:)     ! time and upper boundary temprature (interpolated)

		integer :: n_snow                                       ! number of upper boundary points for snow (input)
		real*8 ,allocatable:: snd_time(:),snd(:,:)              ! upper boundary snow time and snow depth (input)
		integer :: n_stcon
		real*8 ,allocatable:: stcon_time(:),stcon(:,:)          ! snow thermal conductivity time and itself (input)
		real*8 ,allocatable:: snd_i (:,:), stcon_i (:,:)        ! snow depth and thermal conductivity (interpolated)

		integer :: n_ini                                        ! number of vertical grid cells in init file
		real*8, allocatable :: zdepth_ini(:),ztemp_ini(:,:)     ! depth and correspoding initial temperature (time=0) 'zdepth_ini(n_ini)'

		real*8 :: time_restart                                  ! restart time in restart file

		integer :: n_grd                                        ! total number of grid points with depth (grid.txt)
		real*8,allocatable:: zdepth(:),dz(:)                    ! vertical grid and distance between grid point 'zdepth(n_grd)'
		integer :: m_grd                                        ! number of grid points to store in res file
		integer,allocatable:: zdepth_id(:)                      ! index vector of stored grid points 'zdepth_id(m_grid)'

		! thermo physical parameters of soil for each soil layer
		real*8,allocatable:: vwc(:,:)                           ! volumetric water content
		real*8,allocatable:: a_coef(:,:),b_coef(:,:)            ! a and b unfrozen water curve coefficients
		real*8,allocatable:: temp_frz(:,:)                      ! temperature freezing depression
		real*8,allocatable:: EE(:,:)
		real*8,allocatable:: hcap_frz(:,:),hcap_thw(:,:)        ! soil layer heat capacity thawed/frozen
		real*8,allocatable:: tcon_frz(:,:),tcon_thw(:,:)        ! soil layer thermal conductivity thawed/frozen

		integer,allocatable:: n_lay_cur(:)                      ! current number of soil layers <= n_lay
		real, allocatable:: n_bnd_lay(:,:)                      ! number of boundaries between layer in soil

		integer :: n_lay                          				! total allowed number of soil layer

		real*8, allocatable :: temp(:,:)                        ! soil temperature  
		integer,allocatable:: lay_id(:,:)                       ! layer index
		integer,allocatable:: i_time(:)                         ! internal time step with the the main loop

		integer,allocatable:: n_frz_frn(:,:)                    ! number of freezing front (e.g. when freezup is about to happened)
		real*8 ,allocatable:: z_frz_frn(:,:,:)                  ! depth of the freezing front

		real*8 ,allocatable:: RES(:,:)                          ! unified variable for the writing results into the file
		real*8 :: hcap_s                                        ! heat capacity of snow (constant) nondimentional
		real*8 :: L_fus                                         ! Latent heat of fusion [W/mK]

		character(210) :: FMT1,FMT2                             ! results formating type
		real*8 :: TINIR

!the below is from the old script, has to cleaned
          real :: dt
          real :: t
          real :: t_end

          integer :: n_x
          integer :: n_y

          real :: dx
          real :: dy

          real, pointer :: z(:,:)
          real, pointer :: z_temp(:,:)
!from the old script, ends here
      
      end type BMI_Model


      integer, parameter :: component_name_length = 60
      character (len=component_name_length), target :: &
        component_name = "gipl model PermaModel component"

      ! start exchange item list
      integer, parameter :: input_item_count = 1
      integer, parameter :: output_item_count = 1
      integer, parameter :: item_name_length = 22
      character (len=item_name_length), target, &
        dimension (input_item_count) :: &
        input_items = (/'air_temp    '/)

      character (len=item_name_length), target, &
        dimension (output_item_count) :: &
        output_items = (/'ground_temps   '/)
      ! end exchange item list

      contains
          subroutine BMI_Initialize (self, config_file)
            implicit none
            type (BMI_Model), intent (out) :: self
            character (len=*), intent (in) :: config_file
            ! end declaration section
            
			call set_config(self,config_file)
            call set_ids(self)
            call set_upper_bc(self)
            call set_snow_properties(self)
            call set_ic(self)
            call set_grid(self)
            call set_soil_layer_props(self)
            call set_nondim_params(self)
            call assign_layer_id(self)
            call interpolate_ic(self)
            call interpolate_temp_snd(self)
            call open_res_files(self)

          end subroutine BMI_Initialize
          
          
          subroutine set_config(self,cfg_file)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: cfg_file
            character*64 stdummy
		
			if (len (cfg_file)>0) then
			  !Read the configuration file: config_file				  
              open(60,file=cfg_file)
              !read input files
				read(60,'(A)')stdummy
				read(60,'(A)')self%file_sites
				read(60,'(A)')self%file_bound
				read(60,'(A)')self%file_snow
				read(60,'(A)')self%file_rsnow
				read(60,'(A)')self%file_init
				read(60,'(A)')self%file_grid
				read(60,'(A)')self%file_organic
				read(60,'(A)')self%file_mineral

			  !set output files
				read(60,'(A)')stdummy
				read(60,'(A)')stdummy
				read(60,'(A)')self%aver_res_file
				read(60,'(A)')self%result_file
				read(60,'(A)')self%restart_file
              
		      !read input specs
		      	read(60,'(A)')stdummy
		      	read(60,'(A)')stdummy
				read(60,*)self%restart
				read(60,'(A)')stdummy
				read(60,*)self%time_step,self%TAUM,self%TMIN
				read(60,'(A)')stdummy
				read(60,*) self%time_beg,self%time_end
				read(60,'(A)')stdummy
				read(60,*) self%smooth_coef,self%unf_water_coef,self%itmax  
				read(60,'(A)')stdummy
				read(60,*) self%n_sec_day,self%n_time 
				read(60,'(A)')stdummy
				read(60,*) self%sea_level,self%n_frz_max
				read(60,'(A)')stdummy
				read(60,*) self%frz_frn_min,self%frz_frn_max
				read(60,'(A)')stdummy
				read(60,*) self%sat_coef
    		  close(60)
    		  print*, 'Configuration file is read!'
            else            
              self%file_sites="../../examples/sites.txt"
			  self%file_bound="../../examples/bound.txt"
			  self%file_snow="../../examples/snow.txt"
			  self%file_rsnow="../../examples/rsnow.txt"
			  self%file_init="../../examples/initial.txt"

		      self%file_grid="../../examples/grid.txt"
			  self%file_organic="../../examples/organic.txt"
			  self%file_mineral="../../examples/mineral.txt"
			  
			  self%aver_res_file="../../examples/out/mean.txt"
			  self%result_file="../../examples/out/result.txt"
			  self%restart_file="../../examples/out/start.txt"

              self%restart=1
              self%time_step=1.0
              self%TAUM=0.1
              self%TMIN=0.001
              self%time_beg=0
              self%time_end=1
              self%smooth_coef=0.01
              self%unf_water_coef=0.01
              self%itmax=5
              self%n_sec_day=86400.0
              self%n_time=365
              self%sea_level=0
              self%n_frz_max=4
              self%frz_frn_min=0.05
              self%frz_frn_max=10.0
              self%sat_coef=0.95
              print*, 'Configuration file is NOT found!'
              print*, 'Setting up default parameters ...'
            end if            
            
            call filexist(self%file_sites)
			call filexist(self%file_bound)
			call filexist(self%file_snow)
			call filexist(self%file_rsnow)
			call filexist(self%file_grid)
			call filexist(self%file_init)
			call filexist(self%file_mineral)
			call filexist(self%file_organic)
            
          end subroutine set_config
          
          
          subroutine set_ids(self)
            implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: IREAD,ierr,i_site
            character*64 :: stdummy
            
            open(60,FILE=self%file_sites)
            read(60,'(A)')stdummy
		    read(60,*)self%n_site
				allocate(self%snow_code(self%n_site),STAT=IERR)
				allocate(self%veg_code(self%n_site),STAT=IERR)
				allocate(self%geo_code(self%n_site),STAT=IERR)
				allocate(self%gt_zone_code(self%n_site),STAT=IERR)
				allocate(self%temp_grd(self%n_site),STAT=IERR)
				read(60,'(A)')stdummy
				do i_site=1,self%n_site
				  read(60,*) IREAD,self%snow_code(i_site),self%veg_code(i_site),&
				  			self%geo_code(i_site),self%gt_zone_code(i_site),&
				  			self%temp_grd(i_site)
				enddo
    		close(60)
    		print*, trim(self%file_sites),' has been read.'
 		  end subroutine set_ids

          
          subroutine set_upper_bc(self)
            implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: i,ierr,i_site
            character*64 :: stdummy
            
			open(60,file=self%file_bound)
				read(60,'(A)')stdummy
				read(60,*)self%n_temp
				allocate(self%utemp_time(self%n_temp),STAT=IERR)
				allocate(self%utemp(self%n_temp,self%n_site),STAT=IERR)
				read(60,'(A)')stdummy
				do i=1,self%n_temp
					read(60,*) self%utemp_time(i),(self%utemp(i,i_site),i_site=1,self%n_site)
				enddo
			close(60)
		    print*, trim(self%file_bound),' has been read.'
 		  end subroutine set_upper_bc


          subroutine set_snow_properties(self)
            implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: i,ierr,i_site
            character*64 :: stdummy
            
			open(60,file=self%file_rsnow)
			read(60,'(A)')stdummy
			read(60,*)self%n_stcon
			allocate(self%stcon_time(self%n_stcon),STAT=IERR)
			allocate(self%stcon(self%n_stcon,self%n_site),STAT=IERR)
			read(60,'(A)')stdummy
			do i=1,self%n_stcon
				read(60,*) self%stcon_time(i),(self%stcon(i,i_site),i_site=1,self%n_site)
			enddo
			close(60)
		    print*,trim(self%file_rsnow),' has been read.'

			open(60,file=self%file_snow)
			read(60,'(A)')stdummy
			read(60,*)self%n_snow
			allocate(self%snd_time(self%n_snow),STAT=IERR)
			allocate(self%snd(self%n_snow,self%n_site),STAT=IERR)
			read(60,'(A)')stdummy
			do i=1,self%n_snow
			   read(60,*) self%snd_time(i),(self%snd(i,i_site),i_site=1,self%n_site)
			enddo
			  close(60)
		    print*,trim(self%file_snow),' has been read.'
 		  end subroutine set_snow_properties


		  subroutine set_ic(self)
		    implicit none
            type (BMI_Model), intent (inout) :: self
            real*8 ,allocatable :: gtzone(:,:)
            integer :: i,ierr,j,z_num,k
            character*64 :: stdummy
			
			open(60,file=self%file_init,action='read')
			read(60,*)stdummy
			read(60,*)z_num,self%n_ini!,time_restart
			allocate(self%zdepth_ini(self%n_ini),STAT=IERR)
			allocate(self%ztemp_ini(self%n_ini,self%n_site),STAT=IERR)
			allocate(gtzone(self%n_ini,z_num+1),STAT=IERR)
			read(60,*)stdummy
				do i=1,self%n_ini
					read(60,*) (gtzone(i,j),j=1,z_num+1)
				enddo
			close(60)
		    print*,trim(self%file_init),' has been read.'
		    
		    self%time_restart=self%utemp_time(1)
			self%zdepth_ini(:)=gtzone(:,1)
			do i=1,self%n_site
				k=self%gt_zone_code(i)
				self%ztemp_ini(:,i)=gtzone(:,k+1)
			enddo
			deallocate(gtzone)

		  end subroutine set_ic
		  
		  
		  subroutine set_grid(self)
		    implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: i,ierr,j
            character*64 :: stdummy
            
			open(60,file=self%file_grid)
			read(60,'(A)')stdummy
			read(60,*)self%n_grd
			read(60,'(A)')stdummy
			allocate(self%zdepth(self%n_grd),STAT=IERR)
				do i=1,self%n_grd
					read(60,*) self%zdepth(i)
				enddo
				read(60,'(A)')stdummy
				read(60,*)self%m_grd
				read(60,'(A)')stdummy
				allocate(self%zdepth_id(self%m_grd),STAT=IERR)
				do j=1,self%m_grd
					read(60,*)self%zdepth_id(j)
				enddo
			close(60)
		    print*,trim(self%file_grid),' has been read.'
		  end subroutine set_grid
		  
		  	  
		  subroutine set_soil_layer_props(self)
		    implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: i,ierr,j,k,kk
            character*64 :: stdummy
			
			real*8,allocatable:: A1(:,:),A2(:,:),A3(:,:),A4(:,:),A5(:,:)
			real*8,allocatable:: A6(:,:),A7(:,:),A8(:,:),A9(:,:),A10(:,:)
			integer, allocatable :: veg_class(:), num_vl(:)
			integer :: vln
			
			real*8,allocatable:: B1(:,:),B2(:,:),B3(:,:),B4(:,:),B5(:,:)
			real*8,allocatable:: B6(:,:),B7(:,:),B8(:,:)
			integer, allocatable :: geo_class(:), num_gl(:)
			real*8 :: layer_thick
			integer :: gln

			! note: that all max n_lay_cur layers has to be read or it will a give segmantation error
			! n_lay=10!MAXVAL(n_lay_cur)
			self%n_lay=10       
			
			! Read organic and mineral layer props and then combine then together
			! for the entire column
			open (60, file=self%file_organic)
				read(60,'(A)')stdummy
				read(60,*) vln ! reads numbers of  classes
				allocate(A1(self%n_lay,vln),STAT=IERR) ! vwc
				allocate(A2(self%n_lay,vln),STAT=IERR) ! a_coef
				allocate(A3(self%n_lay,vln),STAT=IERR) ! b_coef
				allocate(A4(self%n_lay,vln),STAT=IERR) ! hcap_frz
				allocate(A5(self%n_lay,vln),STAT=IERR)  !hcap_thw
				allocate(A6(self%n_lay,vln),STAT=IERR)  !tcon_frz
				allocate(A7(self%n_lay,vln),STAT=IERR)  !tcon_thw
				allocate(A8(vln,self%n_lay),STAT=IERR)  !bot_cond
				allocate(veg_class(vln),STAT=IERR) !veg_class
				allocate(num_vl(vln),STAT=IERR)  !num_vl number of vegetation layers
				read(60,'(A)')stdummy
				do I = 1,vln
					read(60,*)veg_class(i),num_vl(i)
					read(60,'(A)')stdummy
					do j=1,num_vl(i)
						read(60,*)A1(J,I),A2(J,I),A3(J,I), &
							A4(J,I),A5(J,I),A6(J,I),A7(J,I),A8(I,J)
					enddo
				enddo
			close(60)
		    print*,trim(self%file_organic),' has been read.'
		    
		    open (60, file=self%file_mineral)
		    	read(60,'(A)')stdummy
				read(60,*) gln ! reads numbers of  classes
				allocate(B1(self%n_lay,gln),STAT=IERR) ! vwc
				allocate(B2(self%n_lay,gln),STAT=IERR) ! a_coef
				allocate(B3(self%n_lay,gln),STAT=IERR) ! b_coef
				allocate(B4(self%n_lay,gln),STAT=IERR) ! hcap_frz
				allocate(B5(self%n_lay,gln),STAT=IERR)  !hcap_thw
				allocate(B6(self%n_lay,gln),STAT=IERR)  !tcon_frz
				allocate(B7(self%n_lay,gln),STAT=IERR)  !tcon_thw
				allocate(B8(gln,self%n_lay),STAT=IERR)  !bot_cond
				allocate(geo_class(gln),STAT=IERR) !geo_class
				allocate(num_gl(gln),STAT=IERR)  !num_vl number of lithologic layers
				read(60,'(A)')stdummy
				do i = 1,gln
					read(60,*)geo_class(i),num_gl(i)
					read(60,'(A)')stdummy
					do j=1,num_gl(i)
					read(60,*)B1(J,I),B2(J,I),B3(J,I), &
							B4(J,I),B5(J,I),B6(J,I),B7(J,I),B8(I,J)
					enddo
				enddo
			close(60)
			print*,trim(self%file_mineral),' has been read.'

			allocate(self%vwc(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%a_coef(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%b_coef(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%EE(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%hcap_frz(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%hcap_thw(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%tcon_frz(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%tcon_thw(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%n_lay_cur(self%n_site),STAT=IERR)
			allocate(self%n_bnd_lay(self%n_site,self%n_lay+1),STAT=IERR)
			
			do i = 1,self%n_site
				layer_thick=0
				self%n_bnd_lay(i,1)=layer_thick
				k=self%veg_code(i)
				do j=1,num_vl(k)
					self%vwc(J,I)=A1(j,k);
					self%a_coef(J,I)=A2(j,k);
					self%b_coef(J,I)=A3(j,k);
					self%hcap_thw(J,I)=A4(j,k);
					self%hcap_frz(J,I)=A5(j,k);
					self%tcon_thw(J,I)=A6(j,k);
					self%tcon_frz(J,I)=A7(j,k);
					if (j.eq.1) then 
						layer_thick=A8(k,j)
					else
						layer_thick=layer_thick+A8(k,j);
					endif
					self%n_bnd_lay(i,j+1)=layer_thick
					self%EE(J,I)=0
				!	     write(*,'(3(f8.3),2(f12.1),3(f8.3))') vwc(J,I),a_coef(J,I),b_coef(J,I), &
				!		  hcap_thw(J,I),hcap_frz(J,I),tcon_thw(J,I),tcon_frz(J,I),n_bnd_lay(i,j+1)
				enddo
				k=1
				kk=self%geo_code(i)
				self%n_lay_cur(I)=num_vl(self%veg_code(i))+num_gl(kk) ! maximum number of soil layer = organic layers + mineral layers
				do j=num_vl(self%veg_code(i))+1,self%n_lay_cur(I)
				   self%vwc(J,I)=B1(k,kk);
				   self%a_coef(J,I)=B2(k,kk);
				   self%b_coef(J,I)=B3(k,kk);
				   self%hcap_thw(J,I)=B4(k,kk);
				   self%hcap_frz(J,I)=B5(k,kk);
				   self%tcon_thw(J,I)=B6(k,kk);
				   self%tcon_frz(J,I)=B7(k,kk);
				   self%EE(J,I)=0
				   layer_thick=layer_thick+B8(kk,k);
				   self%n_bnd_lay(i,j+1)=layer_thick!B8(geo_code(i),j)
				   k=k+1
				enddo
				self%n_bnd_lay(i,self%n_lay_cur(I)+1)=self%zdepth(self%n_grd)
			enddo
		
			deallocate(A1,STAT=IERR) ! vwc
			deallocate(A2,STAT=IERR) ! a_coef
			deallocate(A3,STAT=IERR) ! b_coef
			deallocate(A4,STAT=IERR) ! hcap_frz
			deallocate(A5,STAT=IERR)  !hcap_thw
			deallocate(A6,STAT=IERR)  !tcon_frz
			deallocate(A7,STAT=IERR)  !tcon_thw
			deallocate(A8,STAT=IERR)  !bot_cond
			deallocate(veg_class,STAT=IERR) !veg_class
			
			deallocate(B1,STAT=IERR) ! vwc
			deallocate(B2,STAT=IERR) ! a_coef
			deallocate(B3,STAT=IERR) ! b_coef
			deallocate(B4,STAT=IERR) ! hcap_frz
			deallocate(B5,STAT=IERR)  !hcap_thw
			deallocate(B6,STAT=IERR)  !tcon_frz
			deallocate(B7,STAT=IERR)  !tcon_thw
			deallocate(B8,STAT=IERR)  !bot_cond
			deallocate(geo_class,STAT=IERR) !geo_class
			deallocate(num_gl,STAT=IERR)  !num_vl number of lithologic layers
		  end subroutine set_soil_layer_props


		  subroutine set_nondim_params(self)
		    implicit none
            type (BMI_Model), intent (inout) :: self
            integer :: i,ierr,j,i_grd
            real*8, allocatable :: z(:) ! vertical grid
            real*8 :: hcscale
            real*8, parameter  :: hcap_snow=840000.0                ! heat capacity of snow (constant)
			real*8, parameter  :: Lf=333.2*1.D+6					! Latent of water fusion

            
			allocate(z(self%n_grd),STAT=IERR)
			allocate(self%dz(self%n_grd),STAT=IERR)
			allocate(self%temp(self%n_site,self%n_grd),STAT=IERR)
			allocate(self%lay_id(self%n_site,self%n_grd),STAT=IERR)
			allocate(self%i_time(self%n_site),STAT=IERR)
			allocate(self%z_frz_frn(self%n_time,self%n_frz_max,self%n_site),STAT=IERR)
			allocate(self%n_frz_frn(self%n_time,self%n_site),STAT=IERR)
			allocate(self%temp_frz(self%n_lay,self%n_site),STAT=IERR)
			allocate(self%RES(self%n_time,self%m_grd+3),STAT=IERR)
			self%i_time=1  ! active_layer uses it below, needs to be initialized here
			
			z=self%zdepth/self%zdepth(self%n_grd)
			do i_grd=2,self%n_grd
				self%dz(i_grd)=z(i_grd)-z(i_grd-1)
			enddo
			
			hcscale=self%zdepth(self%n_grd)*self%zdepth(self%n_grd)/self%n_sec_day
			self%hcap_frz=self%hcap_frz*hcscale
			self%hcap_thw=self%hcap_thw*hcscale
			self%hcap_s=hcap_snow*hcscale
			self%L_fus=hcscale*Lf
			
			deallocate(z)
		  end subroutine set_nondim_params		


		  subroutine assign_layer_id(self)
			!assigns correspond layer id to the grid point
			!starting from surface to the bottom
			implicit none
			type (BMI_Model), intent (inout) :: self
			integer :: isite,igrd,ilay

			do isite=1,self%n_site
				do 6 igrd=1,self%n_grd
					self%lay_id(isite,igrd)=self%n_lay_cur(isite)
					do ilay=1,self%n_lay_cur(isite)-1
						if (self%n_bnd_lay(isite,ilay).LE.self%zdepth(igrd).AND. &
							self%zdepth(igrd).LT.self%n_bnd_lay(isite,ilay+1))then
							self%lay_id(isite,igrd)=ilay
							GOTO 6
						endif
					enddo
				6 continue
			enddo
			return
		  end subroutine assign_layer_id

		subroutine interpolate_ic(self)
		   implicit none
		   type (BMI_Model), intent (inout) :: self
		   integer :: i,j
		   character*64 :: file_init

			if(self%restart.EQ.1)then !restart=1 means reading initial data from
				do i=1,self%n_site
					call interpolate(self%zdepth_ini,self%ztemp_ini(:,i),&
									self%n_ini,self%zdepth,self%temp(i,:),self%n_grd)
				enddo
			elseif(self%restart.EQ.0)then  			!restart=0 enbales spinup
				write(file_init,'(A14)') 'out/start.txt'
				open(60,file=file_init,action='READ')
					read(60,*)self%time_restart              ! day number in restart file
					do j=1,self%n_grd
						read (60,* ) ( self%temp(i,j),i=1,self%n_site)
					enddo
				close(60)
			endif
		end subroutine interpolate_ic
		

		subroutine interpolate_temp_snd(self)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8 :: vw,a,b
			integer :: i,ierr,j,i_site,i_lay
			integer, parameter :: lbound=2                          ! 1 const temp, 2 heat flux condition at the bottom boundary

			allocate(self%utemp_time_i(self%n_time+2),STAT=IERR)                  ! allocating interval varialbe after interation
			allocate(self%utemp_i(self%n_time+2,self%n_site),STAT=IERR)
			allocate(self%snd_i(self%n_time+2,self%n_site),STAT=IERR)
			allocate(self%stcon_i(self%n_time+2,self%n_site),STAT=IERR)

			do j=1,self%n_time+2
				self%utemp_time_i(j)=self%time_restart+DBLE(j-1)*self%time_step
			enddo
	
			do i_site=1,self%n_site
				if (lbound.EQ.2)self%temp_grd(i_site)=self%temp_grd(i_site)*self%zdepth(self%n_grd)
				do i_lay=1,self%n_lay_cur(i_site)
					vw=self%vwc(i_lay,i_site)
					a=self%a_coef(i_lay,i_site)
					b=self%b_coef(i_lay,i_site)
					self%temp_frz(i_lay,i_site)=-(vw/a)**(1.d0/b)
				enddo
				call interpolate(self%utemp_time,self%utemp(:,i_site),self%n_temp,&
								self%utemp_time_i,self%utemp_i(:,i_site),self%n_time+2)
				call interpolate(self%snd_time,self%snd(:,i_site),self%n_snow,&
								self%utemp_time_i,self%snd_i(:,i_site),self%n_time+2)
				call snowfix(self%utemp_i(:,i_site),self%snd_i(:,i_site),self%n_time+2)
				call interpolate(self%stcon_time,self%stcon(:,i_site),self%n_stcon,&
								self%utemp_time_i,self%stcon_i(:,i_site),self%n_time+2)
				call active_layer(self,i_site)
			enddo

		end subroutine interpolate_temp_snd
		
		subroutine open_res_files(self)
			implicit none
			type (BMI_Model), intent (inout) :: self

			open(1,file=self%result_file,STATUS='unknown')
			open(2,file=self%aver_res_file,STATUS='unknown')
			open(3,file=self%restart_file,STATUS='unknown')
			write(self%FMT1,'(A30,I0,A12)')'(1x,I10,1x,F12.3,2(1x,F16.12),',self%m_grd,'(1x,F16.12))'
			write(FMT11,'(A34,I0,A12)')'(7x,A4,1x,A10,2(3x,A10),',self%m_grd,'(1x,F16.4))'
			write(1,FMT11) 'id','day','bnd_temp','snd', (self%zdepth(self%zdepth_id(i)),i=1,self%m_grd)
			
			write(self%FMT2,'(A28,I0,A40)')'(1x,I10,1x,F12.3,2(1x,F8.3),',self%m_grd,'(1x,F8.3),(1x,F8.3,1x,F12.3),(1x,F12.3))'
			write(FMT22,'(A34,I0,A20,A12)')'(7x,A4,1x,A10,2(3x,A10),',self%m_grd,'(1x,F8.3),3(A12))'
			write(2,FMT22) 'id','day','bnd_temp','snd', (self%zdepth(self%zdepth_id(i)),i=1,self%m_grd),&
							'alt','frz_up_cur','frz_up_tot'
		end subroutine open_res_files

		subroutine active_layer(self,k)
			implicit none
			type (BMI_Model), intent (inout) :: self
			integer :: k,j,jj
			real*8 GA,GB,YFRON,GX,GY
			real*8 fsat

			self%z_frz_frn(self%i_time(k),:,k)=self%sea_level
			self%n_frz_frn(self%i_time(k),k)=0
			do 1329 JJ=1,self%n_grd-1
				J=self%n_grd-JJ
				if (self%zdepth(J).GE.self%sea_level.AND.self%zdepth(J+1).LE.self%frz_frn_max)then
					call fsat_unf_water(self,self%temp(k,J),self%lay_id(k,J),k,GA)
					call fsat_unf_water(self,self%temp(k,J+1),self%lay_id(k,J+1),k,GB)
					if((GA-self%sat_coef)*(GB-self%sat_coef).LE.0.D0) then
						GY=(GA-GB)/(self%zdepth(J)-self%zdepth(J+1))
						GX=(GA+GB-GY*(self%zdepth(J)+self%zdepth(J+1)))/2.D0
						if(GY.EQ.0.D0) then
							YFRON=(self%zdepth(J)+self%zdepth(J+1))/2.D0
						else
							YFRON=(self%sat_coef-GX)/GY
						endif
					else
					GOTO 1329
				endif
				if(self%n_frz_frn(self%i_time(k),k).LT.self%n_frz_max)then
					self%n_frz_frn(self%i_time(k),k)=self%n_frz_frn(self%i_time(k),k)+1
					self%z_frz_frn(self%i_time(k),self%n_frz_frn(self%i_time(k),k),k)=YFRON
					endif
				endif
			1329 CONTINUE

		end subroutine active_layer


		subroutine fsat_unf_water(self,T,NNN,I,sat)!Saturated unforzen water
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent(in) :: T
			integer, intent(in) :: NNN, I
			real*8, intent(out) :: sat
			real*8 :: temp_dep
			real*8 :: a,b,e
			real*8 :: theta
	  
			temp_dep=self%temp_frz(NNN,I) ! freezing temprature depression
			e=self%EE(NNN,I)
			theta=self%vwc(NNN,I)
			a=self%a_coef(NNN,I)
			b=self%b_coef(NNN,I)
			IF(T.LE.temp_dep-e)THEN
				sat=a*((DABS(T))**b)
			ELSEIF(T.GT.temp_dep)THEN
				sat=theta
			ELSE
				sat=a*((DABS(temp_dep-e))**b)
				sat=sat+(theta-sat)*(T+e-temp_dep)/e
			ENDIF
			sat=sat/theta
		end subroutine fsat_unf_water


		subroutine interpolate(XIN,YIN,NIN,XOUT,YOUT,n_itime)
		! Linear interpolation
			implicit none
			real*8, intent(in) :: XIN(NIN),YIN(NIN)
			real*8, intent(out) :: XOUT(n_itime),YOUT(n_itime)
			integer :: i,j,NIN,n_itime
			do I=1,n_itime
				if(XOUT(I).LE.XIN(1))THEN
					YOUT(I)=YIN(1)
					GOTO 1
				elseif(XOUT(I).GT.XIN(NIN))THEN
					YOUT(I)=YIN(NIN)
					GOTO 1
				else
					do J=1,NIN-1
						if (XIN(J).LT.XOUT(I).AND.XOUT(I).LE.XIN(J+1))THEN
							YOUT(I)=YIN(J)+(XOUT(I)-XIN(J))*(YIN(J+1)-YIN(J))/(XIN(J+1)-XIN(J))
							GOTO 1
						endif
					enddo
				endif
				1 continue
			enddo
			return
		end subroutine interpolate


		subroutine snowfix(air_temp,stcon,n)
			implicit none
			real*8, intent (in)  :: air_temp(n)
			real*8, intent (out) :: stcon(n)
			integer :: i,n

			   if(air_temp(1).gt.0.and.stcon(1).gt.0)stcon(1)=0 
			   do i=2,n 
				  if(air_temp(i).gt.0.and.stcon(i).gt.0)then 
				if (stcon(i-1).eq.0)stcon(i)=0 ! puts zeros only at the begining of the year
				  endif
			   enddo

			return
		end subroutine snowfix


          subroutine BMI_Finalize (self)
            implicit none
            type (BMI_Model), intent (inout) :: self
            ! end declaration section

			close(1);close(2);close(3)
			
			deallocate(self%snow_code)
			deallocate(self%veg_code)
			deallocate(self%geo_code)
			deallocate(self%gt_zone_code)
			deallocate(self%temp_grd)
		
			deallocate(self%utemp_time)
			deallocate(self%utemp)
			deallocate(self%stcon_time)
			deallocate(self%stcon)
			deallocate(self%snd_time)
			deallocate(self%snd)
			deallocate(self%zdepth_ini)
			deallocate(self%ztemp_ini)
			deallocate(self%zdepth)
			deallocate(self%zdepth_id)
			
			deallocate(self%vwc)
			deallocate(self%a_coef)
			deallocate(self%b_coef)
			deallocate(self%EE)
			deallocate(self%hcap_frz)
			deallocate(self%hcap_thw)
			deallocate(self%tcon_frz)
			deallocate(self%tcon_thw)
			deallocate(self%n_lay_cur)
			deallocate(self%n_bnd_lay)
			
			deallocate(self%dz)
			deallocate(self%temp)
			deallocate(self%lay_id)
			deallocate(self%i_time)
			deallocate(self%z_frz_frn)
			deallocate(self%n_frz_frn)
			deallocate(self%temp_frz)
			deallocate(self%RES)
			
			deallocate(self%utemp_time_i)       
			deallocate(self%utemp_i)
			deallocate(self%snd_i)
			deallocate(self%stcon_i)

          end subroutine BMI_Finalize

          subroutine BMI_Update (self)
            implicit none
            type (BMI_Model), intent (inout) :: self
            ! variables
			real*8 :: res_save(self%m_grd+3,self%n_site)               ! save results into 2D array
			real*8 :: dfrz_frn(self%n_time)                       ! depth of the freezing front
			real :: frz_up_time_cur                          ! freezeup time current (within a year)
			real :: frz_up_time_tot                          ! freezeup time global
			! counters (time,steps)
			real*8 :: time_s,time_e                          ! internal start and end times
			real*8 :: time_loop                      		 ! main looping time
			real*8 :: time_cur                     			 ! current time (e.g. current day)
			integer :: i_site,j_time,i_grd,i_lay
    
			time_s=self%time_step*DBLE(self%n_time*self%time_beg)
			time_e=self%time_step*DBLE(self%n_time*self%time_end)
			self%i_time=1
			self%TINIR=0.0D0 
			time_loop=0.0D0
			
			do while (time_loop.LT.time_e)
				do i_site=1,self%n_site
					time_cur=time_loop+self%time_restart
					call save_results(self,i_site,time_cur,time_loop)
					6666  continue
					call stefan1D(self,i_site,time_loop)
					time_loop=time_loop+self%time_step
					time_cur=time_loop+self%time_restart
				
					if(self%i_time(i_site).LT.self%n_time)  then
						self%i_time(i_site)=self%i_time(i_site)+1
						call save_results(self,i_site,time_cur,time_loop)
						call active_layer(self,i_site)
						GOTO 6666
					endif
				
					if(time_s.LT.time_e.AND.time_loop.GT.time_s)then
						do j_time=1,self%n_time			! WRITTING RESULTS
							write(1,self%FMT1) i_site, (self%RES(j_time,i_grd),i_grd=1,self%m_grd+3)
						enddo
					endif
					do i_grd=1,self%m_grd+3
						res_save(i_grd,i_site)=sum((self%RES(:,i_grd)))
					enddo
				enddo

				self%i_time=1
				do i_site=1,self%n_site
					frz_up_time_cur=-7777.D0
					frz_up_time_tot=frz_up_time_cur
					do j_time=2,self%n_time
						if((self%n_frz_frn(j_time,i_site)-self%n_frz_frn(j_time-1,i_site)).EQ.-2)then
							if(self%z_frz_frn(j_time-1,self%n_frz_frn(j_time-1,i_site),i_site).GE.self%frz_frn_min) &
								frz_up_time_cur=SNGL(self%RES(j_time,1))
						endif
					enddo

					if(frz_up_time_cur.GT.0.0)then
						frz_up_time_tot=AMOD(frz_up_time_cur,REAL(self%n_time))
						if(frz_up_time_tot.EQ.0.0)frz_up_time_tot=REAL(self%n_time)
					endif
					dfrz_frn=self%z_frz_frn(:,1,i_site)
		
					call save_results(self,i_site,time_cur,time_loop)
					call active_layer(self,i_site)

					!____WRITTING MEAN
					write(2,self%FMT2) i_site,(res_save(i_grd,i_site)/DBLE(self%n_time),i_grd=1,self%m_grd+3), &
								   dfrz_frn(self%n_time),frz_up_time_cur,frz_up_time_tot
					do j_time=1,self%n_time+2
						self%utemp_time_i(j_time)=time_cur+DBLE(j_time-1)*self%time_step
					enddo
					call interpolate(self%utemp_time,self%utemp(:,i_site),self%n_temp,&
									self%utemp_time_i,self%utemp_i(:,i_site),self%n_time+2)
					call interpolate(self%snd_time,self%snd(:,i_site),self%n_snow, &
									self%utemp_time_i,self%snd_i(:,i_site),self%n_time+2)
					call snowfix(self%utemp_i(:,i_site),self%snd_i(:,i_site),self%n_time+2)
					call interpolate(self%stcon_time,self%stcon(:,i_site),self%n_stcon,&
											self%utemp_time_i,self%stcon_i(:,i_site),self%n_time+2)
				enddo
				call save_restart(self)
				self%TINIR=time_loop
			enddo
		
          end subroutine BMI_Update



		subroutine save_restart(self)
			implicit none
			type (BMI_Model), intent (inout) :: self
			integer :: i_site,i_grd
	
			rewind(3) 
			write(3, * ) self%time_restart
			do i_grd=1,self%n_grd
				write (3,* ) ( self%temp(i_site,i_grd),i_site=1,self%n_site)
			enddo    
	
		end subroutine save_restart


		subroutine stefan1D(self,isite,time_loop)
			implicit none
			type (BMI_Model), intent (inout) :: self

			real*8 :: temps(self%n_grd)
			integer:: lay_idx(self%n_grd)
			real*8, intent(inout) :: time_loop
			real*8 :: flux,C_app

			integer :: isite,i_grd,IT

		! tridiagonal variables
			real*8 :: RAB1,RAB2,AKAPA2,AMU2,Q2
			real*8 :: A,B,C,D
			real*8 :: ALF(self%n_grd),BET(self%n_grd)
			real*8 :: EEY,EEY1,abs1,abs2

			real*8 :: temp_o(self%n_grd)             ! old temperature before tridiagonal method
			real*8 :: temp_n(self%n_grd)             ! new temperature after tridiagonal method

		! time counter internal to this subroutine
			real*8 :: time_l                    ! loop time in a subroutine
			real*8 :: time_p                    ! present time in a subroutine
			real*8 :: timei                     ! main subroutine timer
			real :: time_swith                  ! for timei
			real*8 :: lambda
			integer, parameter :: lbound=2   

			lay_idx=self%lay_id(isite,:)
			flux=self%temp_grd(isite)

			time_l=time_loop
			time_swith=-1.0
			timei=self%TAUM
			temps=self%temp(isite,:) 
		64  continue
			time_p=time_l+timei
			temp_o=temps
			IT=1
			ALF(2)=0.D0
			call futemp(self,time_p,isite,BET(2))
			 
		22  continue
			if(IT.GT.self%itmax) then
				timei=timei/2.D0
				time_swith=-1.0
				GOTO 64
			endif

			do i_grd=2,self%n_grd-1
				call fapp_hcap(self,temp_o,isite,i_grd,C_app)
				D=C_app/timei
				call ftcon(self,temp_o(i_grd),isite,i_grd,time_p,lambda)
				A=2.D0*lambda/(self%dz(i_grd)*(self%dz(i_grd)+self%dz(i_grd+1)))
				call ftcon(self,temp_o(i_grd+1),isite,i_grd+1,time_p,lambda)
				B=2.D0*lambda/(self%dz(i_grd+1)*(self%dz(i_grd)+self%dz(i_grd+1)))
				C=A+B+D
				ALF(i_grd+1)=B/(C-A*ALF(i_grd))
				BET(i_grd+1)=(A*BET(i_grd)+D*temps(i_grd))/(C-A*ALF(i_grd))
			enddo
	
			call ftcon(self,temp_o(self%n_grd),isite,self%n_grd,time_p,RAB1)
			call fapp_hcap(self,temp_o,isite,self%n_grd,RAB2)
			
			AKAPA2=2.D0*RAB1/(((RAB2*self%dz(self%n_grd)*self%dz(self%n_grd))/timei+2.D0*RAB1))
			Q2=RAB1*flux
			AMU2=(temps(self%n_grd)*RAB2/timei+2.D0*Q2/self%dz(self%n_grd))/ &
					(RAB2/timei+2.D0*RAB1/self%dz(self%n_grd)**2.D0)
			if(DABS(AKAPA2)>1.D0) then
				print*,'Tridiagonal method is failed - chang you time step tau'
				print*,rab1,rab2,akapa2
				STOP
			endif

		! assigns boundary condition check
			if (lbound.EQ.2)then
				temp_n(self%n_grd)=(AMU2+AKAPA2*BET(self%n_grd))/(1.D0-ALF(self%n_grd)*AKAPA2)
			else
				temp_n(self%n_grd)=flux
			endif
		! calculates new tempratures
			do i_grd=1,self%n_grd-1
				temp_n(self%n_grd-i_grd)=ALF(self%n_grd-i_grd+1)*temp_n(self%n_grd-i_grd+1)+ &
											BET(self%n_grd-i_grd+1)
			enddo
			if(timei>self%tmin) then
				do i_grd=1,self%n_grd
					call fsat_unf_water(self,temp_n(i_grd),lay_idx(i_grd),isite,EEY)
					call fsat_unf_water(self,temp_o(i_grd),lay_idx(i_grd),isite,EEY1)
					abs1=DABS(EEY-EEY1)
					abs2=DABS(temp_o(i_grd)-temp_n(i_grd))
					if((abs1.GT.self%unf_water_coef).or.(abs2.GT.self%smooth_coef)) then
						temp_o=temp_n
						IT=IT+1	
						GOTO 22
					endif
				enddo
			endif
			if(time_p.LT.time_loop+self%time_step-1.D-12)then
				time_l=time_p
				temps=temp_n
				if(time_swith>0) then
					if(timei.LT.self%TAUM) then
						timei=timei*2.D0
						time_swith=-1.0
					endif
				else
					time_swith=1.0
				endif
				GOTO 64
				elseif(time_p.GT.time_loop+self%time_step+1.D-12)then
					timei=(time_loop+self%time_step-time_l)
					goto 64
				else
					temps=temp_n
			endif
		self%temp(isite,:)=temps

		end subroutine stefan1D


		subroutine ftcon(self,T,id,j,time_cur,lambda)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent (inout) :: lambda
			real*8, intent (in) :: T ! scalar here
			real*8, intent (in) :: time_cur
			integer, intent (in) :: id,j
			real*8 :: gr_sur,dsnow,snd,unf_water,WC
			integer :: II, NS

			gr_sur=self%sea_level
			call fsnow_level(self,id,time_cur,snd)
			dsnow=self%sea_level-snd
			NS=self%lay_id(id,j)
			if(self%zdepth(j).le.dsnow)then							!atmosphere
				lambda=1.d4
			elseif (self%zdepth(j).Lt.gr_sur)then					!snow
				II=1+IDINT((time_cur-self%tinir)/self%time_step)
					lambda=self%stcon_i(II,id)+(time_cur+self%time_restart-self%utemp_time_i(II))* &
					(self%stcon_i(II+1,id)-self%stcon_i(II,id))/ &
					(self%utemp_time_i(II+1)-self%utemp_time_i(II))
			else   
				call funf_water(self,T,NS,id,unf_water)				!ground
				WC=unf_water/self%vwc(NS,id)
				lambda=(self%tcon_thw(NS,id)**WC)*(self%tcon_frz(NS,id)**(1.0-WC))
			endif

		end subroutine ftcon


		subroutine fapp_hcap(self,T,I,J,Capp)       ! Apparent heat capacity
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8 :: gr_sur 
			real*8, intent (inout) :: Capp
			integer, intent (in) :: I,J
			real*8, intent (in) :: T(self%n_grd)
			real*8  :: WW(2)
			real*8  :: WC, unf_water,hcap
			integer :: NN(2)
			integer :: li 
	
			li=self%lay_id(I,J)                          ! layer index
			gr_sur=self%sea_level                        ! ground surface
			if(self%zdepth(J).lE.gr_sur)then
				Capp=self%hcap_s                  ! heat capacity for snow
			else
				call funf_water(self,T(J),li,I,unf_water)
				WC=unf_water/self%vwc(li,I)
				Capp=self%hcap_thw(li,I)*WC+self%hcap_frz(li,I)*(1.0-WC)
				if(J.GT.(1).AND.J.LT.self%n_grd)then
					WW(1)=(T(J-1)+T(J))/2.D0
					NN(1)=self%lay_id(I,J-1)
					WW(2)=T(J)
					NN(2)=self%lay_id(I,J)
					call fhcap(self,WW,NN,I,hcap)
					Capp=Capp+hcap*self%dz(J)/(self%dz(J+1)+self%dz(J))
					WW(1)=T(J)
					NN(1)=self%lay_id(I,J)
					WW(2)=(T(J+1)+T(J))/2.D0
					NN(2)=self%lay_id(I,J+1)
					call fhcap(self,WW,NN,I,hcap)
					Capp=Capp+hcap*self%dz(J+1)/(self%dz(J+1)+self%dz(J))
				elseif(J.EQ.1)then
					WW(1)=T(J)
					NN(1)=self%lay_id(I,J)
					WW(2)=(T(J+1)+T(J))/2.D0
					NN(2)=self%lay_id(I,J+1)
					call fhcap(self,WW,NN,I,hcap)
					Capp=Capp+hcap
				elseif(J.EQ.self%n_grd)then
					WW(1)=(T(J-1)+T(J))/2.D0
					NN(1)=self%lay_id(I,J-1)
					WW(2)=T(J)
					NN(2)=self%lay_id(I,J)
					call fhcap(self,WW,NN,I,hcap)
					Capp=Capp+hcap
				endif
			endif

		end subroutine fapp_hcap

		subroutine  fhcap(self,T,nnus,I,hcap)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent(inout) :: hcap
			real*8, intent(in) :: T(2) ! temperature
			integer, intent(in) :: nnus(2)
			real*8  :: H
			real*8  :: dunf_water1,dunf_water2
			real*8  :: fuw11,fuw21,fuw12,fuw22
			integer :: I
		
			H=1/(T(1)-T(2))
			if(DABS(T(1)-T(2)).LT.1.D-6) THEN
				call fdunf_water(self,T(1),NNUS(1),I,dunf_water1)
				call fdunf_water(self,T(2),NNUS(2),I,dunf_water2)
				hcap=0.5d0*(dunf_water1+dunf_water2)
			else
				call funf_water(self,T(1),NNUS(1),I,fuw11)
				call funf_water(self,T(2),NNUS(2),I,fuw22)
				if (nnus(1).ne.nnus(2)) THEN   	
					call funf_water(self,T(2),NNUS(1),I,fuw21)
					call funf_water(self,T(1),NNUS(2),I,fuw12)
					hcap=0.5D0*( H*(fuw11-fuw21)+ H*(fuw12-fuw22) )
				else
					hcap=H*(fuw11-fuw22)
				endif
			endif
			hcap=self%L_fus*DABS(hcap)
		end subroutine  fhcap

		subroutine fdunf_water(self,T,NNN,I,dunf_water)

			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent(in) :: T ! temprature
			integer, intent(in) :: NNN, I
			real*8 :: temp_dep
			real*8 :: a,b,e
			real*8 :: theta
			real*8, intent(inout) :: dunf_water

		temp_dep=self%temp_frz(NNN,I)
		e=self%EE(NNN,I)
		theta=self%vwc(NNN,I)
		a=self%a_coef(NNN,I)
		b=self%b_coef(NNN,I)

		if(T.LE.temp_dep-e)THEN
			dunf_water=-b*a*((DABS(T))**(b-1.0D0))
		elseif(T.GT.temp_dep)THEN
			dunf_water=0.0D0
		else
			dunf_water=a*((DABS(temp_dep-e))**b)
			dunf_water=(b-dunf_water)/e
		endif

		end subroutine fdunf_water


		subroutine funf_water(self,T,NNN,I,fuw)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent(in) :: T ! temperature
			integer, intent(in) :: NNN, I
			real*8 :: temp_dep
			real*8 :: a,b,e
			real*8 :: theta
			real*8, intent(inout) :: fuw

			temp_dep=self%temp_frz(NNN,I) ! change I to k0 everywhere except temp_dep
			e=self%EE(NNN,I)
			theta=self%vwc(NNN,I)
			a=self%a_coef(NNN,I)
			b=self%b_coef(NNN,I)
	
			IF(T.LE.temp_dep-e)THEN
				fuw=a*((DABS(T))**b)
			ELSEIF(T.GT.temp_dep)THEN
				fuw=theta
			ELSE
				fuw=a*((DABS(temp_dep-e))**b)
				fuw=fuw+(theta-fuw)*(T+e-temp_dep)/e
			endif
		end subroutine funf_water

          
		subroutine save_results(self,k,time1,time2)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent (in) :: time1,time2
			integer :: k,j

			self%RES(self%i_time(k),1)=time1
			call futemp(self,time2,k,self%RES(self%i_time(k),2))
			call fsnow_level(self,k,time2,self%RES(self%i_time(k),3))
			do  J=1,self%m_grd
			    self%RES(self%i_time(k),J+3)=self%temp(k,self%zdepth_id(J))
			enddo

		end subroutine save_results
	  
		subroutine futemp(self,T,I,temp)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8, intent (in) :: T
			integer :: I,II
			real*8, intent (out) :: temp

			II=1+IDINT((T-self%TINIR)/self%time_step)
			temp=self%utemp_i(II,I)+(T+self%time_restart-self%utemp_time_i(II)) &
				*(self%utemp_i(II+1,I)-self%utemp_i(II,I))/&
				(self%utemp_time_i(II+1)-self%utemp_time_i(II))

		end subroutine futemp         
          

		subroutine fsnow_level(self,site_id,time,snd)
			implicit none
			type (BMI_Model), intent (inout) :: self
			real*8 :: time
			integer :: site_id,II
			real*8, intent (out) :: snd

			II=1+IDINT((time-self%TINIR)/self%time_step)
			snd=self%snd_i(II,site_id)+(time+self%time_restart-self%utemp_time_i(II))* &
						(self%snd_i(II+1,site_id)-self%snd_i(II,site_id))/ &
						(self%utemp_time_i(II+1)-self%utemp_time_i(II))

		end subroutine fsnow_level 
		
		subroutine filexist(filename)
			implicit none
			character*64 filename
			logical chf

			inquire(file=filename,exist=chf)
			if (.not.chf) then 
				write(*,'(/'' FILE '',a, '' DOESNT EXIST'')')trim(filename)
				stop
			endif
		end subroutine filexist
		!-----------------------------------------------

        subroutine BMI_Update_until (self, t)
            implicit none
            type (BMI_Model), intent (inout) :: self
            real, intent (in) :: t
            ! end declaration section

			print*, 'Nothing here'
			
        end subroutine BMI_Update_until
        

        subroutine BMI_Get_start_time (self, t_start)
            implicit none
            type (BMI_Model), intent (inout) :: self
            real, intent (out) :: t_start
            ! end declaration BMI_Get_start_time

            t_start = self%time_step*DBLE(self%n_time*self%time_beg)     
        end subroutine BMI_Get_start_time

        subroutine BMI_Get_end_time (self, t_end)
            implicit none
            type (BMI_Model), intent (inout) :: self
            real, intent (out) :: t_end
            ! end declaration BMI_Get_end_time

            t_end = self%time_step*DBLE(self%n_time*self%time_end)
        end subroutine BMI_Get_end_time

        subroutine BMI_Get_current_time (self, time)
            implicit none
            type (BMI_Model), intent (inout) :: self
            real, intent (out) :: time
            ! end declaration BMI_Get_current_time

            time = self%t
        end subroutine BMI_Get_current_time

          subroutine BMI_Get_time_step (self, dt)
            implicit none
            type (BMI_Model), intent (inout) :: self
            real, intent (out) :: dt
            ! end declaration section

            dt = self%dt
          end subroutine BMI_Get_time_step

          subroutine BMI_Get_time_units (self, units)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (out) :: units
            ! end declaration section

            units = "-"
          end subroutine BMI_Get_time_units

          subroutine BMI_Get_var_type (self, var_name, type)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            integer, intent (out) :: type
            ! end declaration section

            type = BMI_VAR_TYPE_DOUBLE
          end subroutine BMI_Get_var_type

          subroutine BMI_Get_var_units (self, var_name, units)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            character (len=*), intent (out) :: units
            ! end declaration section

            units = "meter"
          end subroutine BMI_Get_var_units

          subroutine BMI_Get_var_rank (self, var_name, rank)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            integer, intent (out) :: rank
            ! end declaration section

            rank = 2
          end subroutine BMI_Get_var_rank

          subroutine BMI_Get_grid_type (self, var_name, type)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            integer, intent (out) :: type
            ! end declaration section

            type = BMI_GRID_TYPE_UNSTRUCTURED
          end subroutine BMI_Get_grid_type

          subroutine BMI_Get_grid_shape (self, var_name, shape)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            integer, dimension (:), intent (out) :: shape
            ! end declaration section

            shape(1) = self%n_grd
          end subroutine BMI_Get_grid_shape

          subroutine BMI_Get_grid_spacing (self, var_name, spacing)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            real, dimension (:), intent (out) :: spacing
            ! end declaration section

            !spacing(1) = self%dz
          end subroutine BMI_Get_grid_spacing

          subroutine BMI_Get_grid_origin (self, var_name, origin)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            real, dimension (:), intent (out) :: origin
            ! end declaration section

            origin(1) = self%zdepth(1)
          end subroutine BMI_Get_grid_origin

          subroutine BMI_Get_double (self, var_name, dest)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (len=*), intent (in) :: var_name
            real, pointer, intent (inout) :: dest(:)
            ! end declaration section

            real, pointer :: src_as_1d (:)
            type (c_ptr) :: src

            select case (var_name)
              case ('surface_elevation')
                src = c_loc (self%z(1,1))
            end select

            if (associated (dest)) then
              call C_F_POINTER (src, src_as_1d, [self%n_x*self%n_y])
              dest = src_as_1d
            else
              call C_F_POINTER (src, dest, [self%n_x*self%n_y])
            endif

          end subroutine BMI_Get_double

          subroutine BMI_Get_double_copy (self, var_name, dest)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (len=*), intent (in) :: var_name
            real, intent (out) :: dest (*)
            ! end declaration section

            select case (var_name)
              case ('surface_elevation')
                call copy_array (dest, self%z, self%n_x * self%n_y)
            end select

          end subroutine BMI_Get_double_copy

          subroutine BMI_Get_double_at_indices (self, var_name, dest, inds)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (len=*), intent (in) :: var_name
            real, pointer, intent (inout) :: dest(:)
            integer, intent (in) :: inds(:)
            ! end declaration section

            real, pointer :: src_as_1d (:)
            type (c_ptr) :: src
            integer :: i

            select case (var_name)
              case ('surface_elevation')
                src = c_loc (self%z(1,1))
            end select

            call C_F_POINTER (src, src_as_1d, [self%n_x*self%n_y])

            do i = 1, size (inds)
              dest(i) = src_as_1d(inds(i))
            end do

          end subroutine BMI_Get_double_at_indices

          subroutine BMI_Set_double (self, var_name, src)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            real, intent (in) :: src (*)
            ! end declaration section

            call copy_array (self%z, src, self%n_x * self%n_y)

          end subroutine BMI_Set_double

          subroutine BMI_Set_double_at_indices (self, var_name, inds, src)
            implicit none
            type (BMI_Model), intent (inout) :: self
            character (len=*), intent (in) :: var_name
            integer, intent (in) :: inds(:)
            real, intent (in) :: src (*)
            ! end declaration section

            real, pointer :: dest_as_1d (:)
            type (c_ptr) :: dest
            integer :: i

            select case (var_name)
              case ('surface_elevation')
                dest = c_loc (self%z(1,1))
            end select

            call C_F_POINTER (dest, dest_as_1d, [self%n_x*self%n_y])

            do i = 1, size (inds)
              dest_as_1d(inds(i)) = src(i)
            end do

          end subroutine BMI_Set_double_at_indices

          subroutine copy_array (dest, src, n)
            real, intent (out) :: dest(n)
            real, intent (in) :: src(n)
            integer, intent (in) :: n
            dest = src
          end subroutine copy_array

          subroutine BMI_Get_input_var_names (self, names)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (*), pointer, intent (out) :: names(:)
            ! end declaration section

            names => input_items
          end subroutine BMI_Get_input_var_names

          subroutine BMI_Get_output_var_names (self, names)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (*), pointer, intent (out) :: names(:)
            ! end declaration section

            names => output_items
          end subroutine BMI_Get_output_var_names

          subroutine BMI_Get_component_name (self, name)
            implicit none
            type (BMI_Model), intent (in) :: self
            character (len=*), pointer, intent (out) :: name
            ! end declaration section

            name => component_name
          end subroutine BMI_Get_component_name

      end module

