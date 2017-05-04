      program test
      use bmif
      implicit none

      type (BMI_Model) :: m
      integer :: rank
      integer, allocatable, dimension (:) :: shape
      real, allocatable, dimension (:) :: spacing
      real, allocatable, dimension (:) :: origin
      integer :: i
      real :: time
      character (len=22), pointer :: names(:)
      character (len=component_name_length), pointer :: name

      write (*,"(A)") "Initializing..."
      call BMI_Initialize (m, "../../examples/gipl_cfg.cfg")

      write (*,"(A)") "Component info:"
      call BMI_Get_component_name (m, name)
      write (*,"(a30, a30)") "Component name: ", name

      call BMI_Get_start_time (m, time)
      write (*,"(A30, f8.2)") "Start time: ", time
      
      call BMI_Get_end_time (m, time)
      write (*,"(A30, f8.2)") "End time: ", time


      write (*,"(A)") "Running..."
      call BMI_Update (m)
      write (*,"(A)") "Done."


      write (*,"(A)", advance="no") "Finalizing..."
      call BMI_Finalize (m)
      write (*,*) "Done"

      end program test

