program TopoStresses_run

    implicit none

    integer nx,ny,nz
    double precision xl,yl,zl
    double precision rhog,ym,pr,ks,exx,eyy
    integer ibc
    double precision, target, dimension(:), allocatable :: h
    double precision, target, dimension(:), allocatable :: stress
    double precision, dimension(:,:), pointer, contiguous :: h2
    double precision, dimension(:,:,:,:), pointer, contiguous :: stress2

    nx=40 ! model resolution in the x-direction
    ny=40 ! model resolution in the y-direction
    nz=5 ! model resolution in the z-direction

    allocate (h(nx*ny),stress(6*nx*ny*nz)) ! allocates memory to store the surface topography and the stress tensor in a 1D array
    h2(1:nx,1:ny)=>h ! creates a 2D array equivalenet to h
    stress2(1:6,1:nz,1:nx,1:ny)=>stress ! creates a 4D array equivalent to stress
    ! stress2(1,:,:,:) contains the xx-component of the stress tensor at all points of the nx x ny x nz grid
    ! stress2(2,:,:,:) contains the yy-component of the stress tensor at all points of the nx x ny x nz grid
    ! stress2(3,:,:,:) contains the zz-component of the stress tensor at all points of the nx x ny x nz grid
    ! stress2(4,:,:,:) contains the xy-component of the stress tensor at all points of the nx x ny x nz grid
    ! stress2(5,:,:,:) contains the zx-component of the stress tensor at all points of the nx x ny x nz grid
    ! stress2(6,:,:,:) contains the yz-component of the stress tensor at all points of the nx x ny x nz grid

    open (101,file='Topo_test.txt',status='old') ! open the topo-test file
    read (101,*) h ! reads the test topography
    close (101) ! closes the topo-test file

    xl=10.d3 ! size of the model in the x-directoin (in m)
    yl=10.d3 ! size of the model in the y-direction (in m)
    zl=10.d3 ! size of the model in the z-direction (in m)
    rhog=2800.d0*9.81d0 ! product of rho (density in kg/m^3) by g (in m/s^2)
    ym=1.d11 ! Young's modulus (in Pa)
    pr=0.25d0 ! Poisson's ratio
    ks=-1.d0 ! yield stress (in Pa)
    exx=0.d0 ! imposed strain in the x-direction
    eyy=0.d0 ! imposed strain in the y-direction
    ibc=2 ! side bondary conditions
    ! ibc=0 is free displacements 
    ! ibc=1 is no-slip
    ! ibc=2 is free-slip
    ! ibc=3 is no-tilt

    call TopoStresses_Setup(nx,ny,nz) ! initializes and sends the model resolution through the interface
    call TopoStresses_Set_xlylzl(xl,yl,zl) ! sends the model sizes through the interface
    call TopoStresses_Set_Mechanical_Properties(rhog,exx,eyy,ym,pr,ks) ! sends the mechanical properties through the interface
    call TopoStresses_Set_BC(ibc) ! setes the boundary conditions through the interface
    call TopoStresses_Set_H(h) ! sends the surface topography through the interface
    call TopoStresses_Execute() ! executes TopoStress
    call TopoStresses_Copy_stress(stress) ! recover the computed stresses through the interface
    print*,minval(stress2(4,:,:,:)),maxval(stress2(4,:,:,:)) ! prints the minimum and maximum values of sigma_xy 
    call TopoStresses_Destroy() ! releases memory through the interface

end program TopoStresses_run