module variables


  integer,parameter :: MAX_CONS_TYPE=15
  real*8,dimension(:),allocatable::shellcharge
  real*8,dimension(:,:),allocatable:: xyz_pts,xyz_molecule,p2p_response_static
  real*8,dimension(:,:,:),allocatable:: xyz_drude_eq
  character(3),dimension(:),allocatable::atomtype,atom_label,atmtype
  integer,parameter::rank=3
  character(3) :: hard_constraints="yes"  ! if we are plugging in definite values for pol tensors 
  real*8,dimension(:),allocatable :: static_pols
  character(3),dimension(MAX_CONS_TYPE) :: hard_cons_type
  real*8,dimension(MAX_CONS_TYPE) :: hard_cons_values
  real*8,parameter :: springconstant=.1 ! au
  real*8,parameter :: thole = 2.0,deltashell=1D-8
  real*8, parameter :: shell_min=0.0
  

end module variables
