program numerov_main

    implicit none
    real(8), parameter :: etol = 1d-9
    integer, parameter :: mxprop = 200
    real(8) :: amass, xmin, xmax, dx, enrl, enrh, deltae, energy, anorm
    ! These are arrays with variable size which will be allocated in run time
    real(8), allocatable :: x(:), psi(:), highe, lowe
    integer, allocatable :: nodes(:)
    integer :: i, ngrid, nodes, nodes_h, nodes_l, ilevel, nlevel, iprop
    logical :: lnodecnt

    ! Namelist variables to be read from the input file 
    namelist/input/ amass, xmin, xmax, dx, enrl, enrh

    ! Default values
    amass = 1d0   ! amu
    xmin  = -10d0 ! ang
    xmax  = 10d0  ! ang
    dx    = 0.001d0 ! and
    enrl  = 0d0   ! cm-1
    enrh  = 10d0  ! cm-1

    ! Now read parameters of input file
    read(*,input)

    ! Sanity check
    if (enrl>enrh) stop "ENRL should be less than ENRH"

    ngrid = int((xmax-xmin)/dx)+1

    ! Allocating grid and wavefunction arrays
    allocate(x(ngrid))
    allocate(psi(ngrid))

    ! Defining x grid
    do i = 1, ngrid
       x(i) = xmin+(i-1)*dx
    enddo
    
    ! Calculate node count for enrh by calling this subroutine with lnodecnt as .true.
    lnodecnt = .true. 
    call numerov(ngrid, x, dx, amass, enrh, psi, nodes_h, lnodecnt)
    write(*,*)"Number of nodes for ENRH is:", nodes_h

    ! Calculate node count for enrl by calling this subroutine with lnodecnt as .true.
    lnodecnt = .true. 
    call numerov(ngrid, x, dx, amass, enrl, psi, nodes_l, lnodecnt)
    write(*,*)"Number of nodes for ENRL is:", nodes_l
    
    ! The number of levels must be the difference between the 2 node counts nodes_h and nodes_l
    nlevel = nodes_h - nodes_l

    if (nlevel==0) stop "No bound states available in the given energy window"

    write(*,*)"This means number of levels within this energy window is:", nlevel
    
    ! nodes array will contain the node count information for each level. So nodes is allocated with size nlevel
    allocate(nodes(nlevel))
    
    !The low energy boundary will change for each level calculation
    lowe  = enrl

    ! This loop does the energy and wavefunction calculation for each level
    do ilevel = 1, nlevel 

       !Store the node count in nodes for each level
       nodes(ilevel) = nodes_l+ilevel-1

       !The high energy boundary should be same for each level calculation
       highe = enrh

       !The energy difference deltae will be minimized to etol(etolerance)
       deltae = highe-lowe
       write(*,*)
       write(*,*) 'Solving for bound state number:', ilevel

       write(*,*) 'Propagation number:'

       !do while(abs(deltae)>etol) ! bisection loop
       do iprop = 1, mxprop
       
         !Guess energy is the average value
         energy = (highe+lowe)/2d0

         !Calling this subroutine bcoz to calculate node count for this energy
         lnodecnt = .true.
         call numerov(ngrid, x, dx, amass, energy, psi, nodes, lnodecnt)
         
         !If the node count is greater than the actual node which is stored in
         !nodes(ilevel), then shift highe to this energy, otherwise shift the 
         !lowe to this energy
         if (nodes > nodes(ilevel)) then
           highe = energy
         else
           lowe = energy
         endif

         deltae = highe-lowe
         
         ! Write some information
         write(*,*) iprop, ", Nodes:", nodes, ", Energy:", energy

         !Since we need the difference of energy to a perticular small value (etol),
         !for which we dont need to go for mxprop (all iteration), so we exit from the loop 
         if (abs(deltae) <= etol) exit
       enddo ! bisection loop ends

       ! Energy is now converged to the eigenvalue for this level
       ! Now call numerov again to calculate the wave function for this energy
       lnodecnt = .false.
       call numerov(ngrid, x, dx, amass, energy, psi, nodes, lnodecnt)

       ! Calculate normalization constant
       anorm = 0d0
       do i = 1, ngrid
          anorm = anorm+psi(i)*psi(i)*dx
       enddo

       write(20+ilevel,'(" # Energy (cm-1):", f12.6)') energy
       write(20+ilevel,*)"# No. of nodes:", nodes(ilevel)
       write(20+ilevel,*)"# Wave function:"
       do i = 1, ngrid
          write(20+ilevel,'(f18.6,es18.6)') x(i), psi(i)/sqrt(anorm)
       enddo

     enddo ! levels loop ends

     print*, ""
     print*, "Calculation complete"

end program numerov_main





    subroutine numerov(ngrid, x, dx, amass, energy, psi, nodes, lnodecnt)
    implicit none
    real(8), parameter :: bfct = 16.85762919164018d0  !cm^-1 amu ang^2
    real(8), parameter :: psiinit = 1d-4
    integer, intent(in) :: ngrid
    !lnodecnt(logical) can take true or false 
    logical, intent(in) :: lnodecnt
    real(8), intent(in) :: amass, energy, x(ngrid), dx
    real(8), intent(out) :: psi(ngrid)
    integer, intent(out) :: nodes 
    real(8) :: f(ngrid), psi_r(ngrid), psi_l(ngrid), dx2, prefac, signl, xmatch
    integer :: i
    real(8), external :: potential

    dx2 = dx*dx
    prefac = amass/bfct
    !Define the array f
    do i = 1, ngrid
      f(i) = prefac*(energy-potential(x(i))) 
    enddo

    ! The following quantity should not be less than -1 in order to avoid any
    ! spurious node in the wave function 
    if (dx2*f(1)/12d0 <= -1d0) stop "XMIN too small or dx too large in numerov"
    if (dx2*f(ngrid)/12d0 <= -1d0) stop "XMAX too large or dx too large in numerov"

    ! Propagation should always start or end at classically forbidden region
    if (f(1) > 0d0) stop "XMIN should not be in classically allowed region"
    if (f(ngrid) > 0d0) stop "XMAX should not be in classically allowed region"

    ! Count the number of nodes
    nodes = 0
    if (lnodecnt) then
       ! Initialization 
       psi(1) = 0d0
       psi(2) = psiinit
        
       ! Numerov Propagation Loop
       do i = 2, ngrid-1
           ! See Eq. (7) of numero_algo.pdf in Dropbox/BMGroup/Codes/Numerov/
           psi(i+1) = (2.0d0 * (1.0d0 - 5.0d0/12.0d0 * dx2 * f(i)) * psi(i) &
                      - (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i-1)) * psi(i-1)) &
                      / (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i+1))
           if (psi(i+1)*psi(i) < 0d0) nodes = nodes+1
       end do
       return  
       ! Subroutine should end here if only node count was wanted
    end if

    ! Subroutine continues here if wave function is needed

    !Initialize right wavefunction for inward propagation
    psi_r = 0d0
    psi_r(ngrid-1) = psiinit
   
    ! Numerov inward propagation Loop (right to left)
    do i = ngrid-1, 2, -1
      psi_r(i-1) = (2.0d0 * (1.0d0 - 5.0d0/12.0d0 * dx2 * f(i)) * psi_r(i) &
                 - (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i+1)) * psi_r(i+1)) &
                 / (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i-1))
      xmatch = x(i)
      if (f(i) > 0d0) exit  
    end do

    !Initialize left wavefunction for outward propagation
    psi_l = 0d0
    psi_l(2) = psiinit
   
    ! Numerov outward propagation Loop (left to right)
    do i = 2, ngrid-1
      if (x(i) > xmatch) exit
      psi_l(i+1) = (2.0d0 * (1.0d0 - 5.0d0/12.0d0 * dx2 * f(i)) * psi_l(i) &
                 - (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i-1)) * psi_l(i-1)) &
                 / (1.0d0 + 1.0d0/12.0d0 * dx2 * f(i+1))
      signl = sign(1d0,psi_l(i))     
    end do

    do i = 1, ngrid
      if (x(i) <= xmatch) psi(i) = psi_l(i) 
      if (x(i) > xmatch) psi(i) = signl*psi_r(i) 
    enddo
    return

    end subroutine numerov

  
