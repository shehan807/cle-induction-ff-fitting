module routines

contains

  !***************************************************
  ! this subroutine reads a .clt file and obtains the molecule coordinates and atom types
  !**********************************************
  subroutine get_molecule_data(molecule_file)
    use variables

    character(*),intent(in) :: molecule_file
    character(50) :: line
    integer::inputstatus,ind,count,i,j
    real*8 :: junkr
    character(4) :: junkc


    open(unit=7,file=molecule_file,status='old')

    ! get number of atom types

    do
       read(7,'(A)',Iostat=inputstatus) line
       if(inputstatus < 0 ) then
        write(*,*) "Can't find Charge section in clt file!"
        stop
       endif
       ind=index(line,'Charge')
       if(ind .ne.0) exit
    enddo



    count=0
    do
       read(7,'(A)',Iostat=inputstatus) line
       if(inputstatus < 0 ) exit
       ind=index(line,'End')
       if(ind .ne.0) exit
       count=count+1
    enddo

    allocate(atomtype(count),atom_label(count),xyz_molecule(count,3))

    do i=1,count+1
       backspace(7)
    enddo

    ! fill in atom types
    do i=1,count
       read(7,*) atom_label(i),junkr,(xyz_molecule(i,j),j=1,3),junkc,atomtype(i)
    enddo

    close(7)

  end subroutine get_molecule_data

  !*****************************************************
  ! this subroutine reads the positions and response values for the
  ! point to point response sites from the camcasp .p2p binary files
  !******************************************************
  subroutine get_p2p(p2pbase)
    use variables
    character(*),intent(in)::p2pbase
    character(70) :: response_file
    character(1)  :: temp1
    character(2)  :: temp2
    integer*4 :: n_points
    integer:: i,j,k

    ! first open the static file
    response_file=trim(p2pbase)//"_000.p2p"

    write(*,*) ""
    write(*,*) "reading p2p response file ", response_file


    open (unit=7,file=response_file,form='unformatted',access='sequential')

    read(7) n_points

    allocate(xyz_pts(n_points,3),p2p_response_static(n_points,n_points))

    do i=1,n_points
       read(7) ( xyz_pts(i,j),j=1,3)
    enddo

    do i=1,n_points
       do j=1, i
          read(7) p2p_response_static(i,j)
       enddo
    enddo

    close(7)

  end subroutine get_p2p


  !****************************************************
  ! this subroutine reads in frequency dependent polarizability
  ! values that will be hard constraints
  !****************************************************
  subroutine get_hard_constraints(ifile)
    use variables
    character(*),intent(in) :: ifile
    character(50) :: line
    integer :: count,i,j,rank_l,flag


    open(unit=7,file=ifile,status='old')

    read(7,'(A)') line
    read(7,'(A)') line
    read(7,'(A)') line

    count=1
    ! keep reading hard constraint parameters until end of file
    do
       read(7,*,Iostat=inputstatus) hard_cons_type(count)
       if(inputstatus < 0 ) exit
       i=1
       read(7,*) rank_l
       read(7,*) hard_cons_values(count)
       count=count+1
       if (count .ge. MAX_CONS_TYPE ) then
          stop "more atom types in hard constraint file than MAX_CONS_TYPE"
       endif
    enddo

    close(7)

    ! now check that all atom types have been given initial values
    do i=1, size(atomtype)
       flag=0
       do j=1,count
          if ( atomtype(i) .eq. hard_cons_type(j) ) then
             flag=1
          endif
       enddo
       if ( flag == 0 ) then
          write(*,*) "Please input an initial polarizability for atomtype ", atomtype(i)
          stop
       endif
    enddo

  end subroutine get_hard_constraints


  !******************************************************
  ! this subroutine fills in hard values for the distributed 
  ! polarizability tensors of atoms that match labels in
  ! the hard_constraint input file
  !******************************************************
  subroutine fill_in_hard_constraints
    use variables
    integer :: i,j

    do i=1,size(atomtype)
       do j=1,size(hard_cons_type)
          if (atomtype(i) .eq. hard_cons_type(j) ) then
             ! remember, parameters can be negative if we want to constrain them
             ! so take absolute value here
                static_pols(i) = abs(hard_cons_values(j))
          endif
       enddo
    enddo

  end subroutine fill_in_hard_constraints


  !******************************************************
  ! this subroutine determines the number of independent atomtypes
  ! as well as the number of hydrogen atoms
  !
  ! if parameter is initially set to zero or negative value, don't fit it
  ! if it's negative use absolute value of that number as a hard constraint
  !******************************************************
  SUBROUTINE numindatoms(atms)
    use variables
    INTEGER,INTENT(out)::atms
    CHARACTER(2)::value
    INTEGER::count1,count2,atoms,count,i,j,ind
    integer,dimension(MAX_CONS_TYPE):: hard_cons_zero
    CHARACTER(3),DIMENSION(:),ALLOCATABLE::typemon1

    atoms=SIZE(atomtype)
    hard_cons_zero=1

    ALLOCATE(typemon1(atoms))

    DO i=1,atoms
       WRITE(typemon1(i),'(A3)') atomtype(i)
    ENDDO

    count1=0
    count2=0
    DO i=1,atoms
       if ( hard_constraints .eq. "yes" ) then
!!!!! fit atom type if it is initially greater than zero
          do j=1,size(hard_cons_type)
             if(atomtype(i).eq.hard_cons_type(j)) then
                if ( hard_cons_values(j) > 0d0 ) then
                   hard_cons_zero(j)=0
                endif
             endif
          enddo
       endif
    ENDDO


    ! number of atom types fit is number of non zero initial values
    if ( hard_constraints .eq. "yes" ) then
       do j=1,size(hard_cons_type)
          if ( hard_cons_zero(j) .eq. 0 ) then
             count2=count2 + 1
          endif
       enddo
    endif

    atms=count2


  END SUBROUTINE numindatoms


  subroutine collapseatomtype(h_atms)
    use variables
    INTEGER,INTENT(out)::h_atms
    INTEGER::count,countt,count1,count2,atoms,atms,i,j,k,ind

    atoms=SIZE(atomtype)
    atms=SIZE(atmtype)

    count=0
    countt=0
    count1=0

    DO i=1,atms
       count=countt
       DO j=i+count,atoms
!!! check if we are fitting this atom type
          if ( hard_constraints .eq. "yes" ) then
             do k=1,size(hard_cons_type)
                if(atomtype(j).eq.hard_cons_type(k)) then
                   ! fitting if non zero initial value
                   if ( hard_cons_values(k) > 0d0 ) then
                      ! nothing
                   else
                      countt=countt+1
                      go to 202 
                   endif
                endif
             enddo
          endif
          DO k=1,j-1
             IF(atomtype(k).EQ.atomtype(j)) THEN
                countt=countt+1
                go to 202
             ENDIF
          ENDDO
          atmtype(i)=atomtype(j)
          !see if this is a hydrogen
          ind=index(atomtype(j),"H")
          if (ind .ne. 0 ) count1=count1+1

          go to 203
202       CONTINUE
       ENDDO
203    CONTINUE
    ENDDO

    h_atms = count1 


  END SUBROUTINE collapseatomtype




  subroutine fill_charges(param)
    use variables
    real*8,dimension(:),intent(in) :: param

    integer :: i,j,atms,atoms

    atms=size(atmtype)
    atoms=size(atomtype)


    do i=1, atoms
       do j=1, atms
          if ( atmtype(j) .eq. atomtype(i) ) then
             shellcharge(i)= -abs(param(j))
          endif
       enddo
    enddo

  end subroutine fill_charges



end module routines
