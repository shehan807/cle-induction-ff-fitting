!*****************************************************************
! this program is for fitting Drude oscillators to the
! point to point response calculations performed by camcasp
!*****************************************************************

program p2p
  use variables
  use routines
  use eq_drude
  use amoeba_routines

  character(50) :: name
  character(70) :: consfile, molecule_file
  character(8)::date
  character(10)::time
  integer:: i,j,k,n_param,atoms,atms,h_atms,iter
  real*8,dimension(:),allocatable::p0,p1,df0,delta
  real*8 :: f0,f1,fdel,rms,max_d,min_d
  real*8,parameter :: small=1d-10
  real*8:: r_tol=5D-4
  real*8,dimension(:),allocatable:: y_amoeba
  real*8,dimension(:,:),allocatable::p_amoeba

  call getarg(1,name); 

  write(*,*) "*******************************************"
  write(*,*) "IMPORTANT:  All units are assumed to be in a.u."
  write(*,*) "please make sure geometry in .clt file is in"
  write(*,*) "atomic units"
  write(*,*) "*******************************************"
  write(*,*) ""

  consfile=trim(name)//".constraint"
  write(*,*) ""
  write(*,*) "reading constraint file ", consfile


  molecule_file=trim(name)//".clt"
  write(*,*) ""
  write(*,*) "reading clt file ", molecule_file
  call get_molecule_data(molecule_file)
  ! atom types will be fit if initial values are non zero, and will not be fit if initial values are zero
  atoms=size(atomtype)
  allocate(static_pols(atoms))

  write(*,*) "atomtypes"
  do i=1,atoms
    write(*,*) atomtype(i)
  enddo

  ! get hard_constraints
  if ( hard_constraints .eq. "yes" ) then
     call get_hard_constraints(consfile)
  endif

  ! need to fill in hard constraints
  if ( hard_constraints .eq. "yes" ) then
     call fill_in_hard_constraints
  endif

  call numindatoms(atms)
  allocate(atmtype(atms))
  call collapseatomtype(h_atms)

 write(*,*) "number of atomtypes that will be fit", atms
 do i=1,atms
   write(*,*) i, atmtype(i)
 enddo

  ! get static point to point response
  call get_p2p(name)
  write(*,*) "finished reading p2p"


  write(*,*) "initializing drude oscillators..."

  ! initialize drude oscillators
  call initialize_drude_oscillators
  call rms_p2p(rms,max_d,min_d)

  write(*,*) "initializing amoeba data..."
  ! initialize amoeba data
  allocate( p_amoeba(atms+1,atms), y_amoeba(atms+1) )
  do i=1,atms
     do j=1,atoms
        if ( atomtype(j) .eq. atmtype(i) ) then
           p_amoeba(:,i) = shellcharge(j)
        endif
     enddo
  enddo
  y_amoeba(1) = rms



  do i= 2, atms+1
     p_amoeba(i,i-1) = 1.5d0 * p_amoeba(i,i-1)
     call fill_charges(p_amoeba(i,:))

     call equilibrate_drude_positions
     call rms_p2p(rms,max_d,min_d)
     y_amoeba(i)=rms

  enddo


  call amoeba(p_amoeba,y_amoeba,r_tol,iter) 


  call rms_p2p(rms,max_d,min_d)

  write(*,*) "converged Drude charges"
  do i=1,atms
     do j=1,atoms
        if ( atomtype(j) .eq. atmtype(i) ) then  
           write(*,*) atmtype(i), shellcharge(j)
        endif
     enddo
  enddo
  write(*,*) "rms, max_d,min_d"
  write(*,*) rms, max_d,min_d
  stop



end program p2p
