!---------------

subroutine vtk (xl,yl,zl,topo,nx,ny,nz,ks,stress)

    ! subroutine to store the stresses in a vtk file for display
    ! use Paraview or similar software to open Stress.vtk and display the computed stresses
    ! in 3D
    
    implicit none
    
    integer iunit,nn,nx,ny,nz,i,j,k,ne,ie,ijk
    double precision xl,yl,zl,ks
    double precision stress(6,(nx-1)*(ny-1)*(nz-1)),topo(nx,ny)
    double precision r,s,t,h(8),xint,yint,zint
    double precision, dimension(:), allocatable :: x,y,z
    double precision, dimension(:), allocatable :: p,s2d,fail
    integer, dimension(:,:), allocatable :: icon
    
    nn=nx*ny*nZ
    ne=(nx-1)*(ny-1)*(nz-1)
    
    allocate (x(nn),y(nn),z(nn))
    allocate (p(ne),s2d(ne),fail(ne))
    
    do k=1,nz
      do j=1,ny
        do i=1,nx
        ijk=i+(j-1)*nx+(k-1)*nx*ny
        x(ijk)=xl*float(i-1)/(nx-1)
        y(ijk)=yl*float(j-1)/(ny-1)
        z(ijk)=(zl+topo(i,j))*float(k-1)/(nz-1)
        enddo
      enddo
    enddo
    
    allocate (icon(8,ne))
    
    ie=0
      do k=1,nz-1
        do j=1,ny-1
          do i=1,nx-1
          ie=ie+1
          icon(1,ie)=i+(j-1)*nx+(k-1)*nx*ny
          icon(2,ie)=icon(1,ie)+1
          icon(3,ie)=icon(1,ie)+nx
          icon(4,ie)=icon(1,ie)+nx+1
          icon(5,ie)=icon(1,ie)+nx*ny
          icon(6,ie)=icon(2,ie)+nx*ny
          icon(7,ie)=icon(3,ie)+nx*ny
          icon(8,ie)=icon(4,ie)+nx*ny
          enddo
        enddo
      enddo
    
    p=(stress(1,:)+stress(2,:)+stress(3,:))/3.d0
    s2d=sqrt(0.5d0*((stress(1,:)-p)**2+(stress(2,:)-p)**2+(stress(3,:)-p)**2)+ &
              stress(4,:)**2+stress(5,:)**2+stress(6,:)**2)
    fail=s2d**2-ks**2
    fail=-min(0.,fail)
    
    iunit=33
    
    open(unit=iunit,file='Stress.vtk')
    write(iunit,'(a)')'# vtk DataFile Version 3.0'
    write(iunit,'(a)')'Stress'
    write(iunit,'(a)')'ASCII'
    write(iunit,'(a)')'DATASET STRUCTURED_GRID'
    write(iunit,'(a11,i4,i4,i4)')'DIMENSIONS ',(nx-1),(ny-1),(nz-1)
    write(iunit,'(a7,i10,a6)')'POINTS ',ne,' float'
    do k=1,nz-1
      do j=1,ny-1
        do i=1,nx-1
        ie=(k-1)*(nx-1)*(ny-1)+(j-1)*(nx-1)+i
        r=0.d0
        s=0.d0
        t=0.d0
        h(1)=(1.d0-r)*(1.d0-s)*(1.d0-t)/8.d0
        h(2)=(1.d0+r)*(1.d0-s)*(1.d0-t)/8.d0
        h(3)=(1.d0-r)*(1.d0+s)*(1.d0-t)/8.d0
        h(4)=(1.d0+r)*(1.d0+s)*(1.d0-t)/8.d0
        h(5)=(1.d0-r)*(1.d0-s)*(1.d0+t)/8.d0
        h(6)=(1.d0+r)*(1.d0-s)*(1.d0+t)/8.d0
        h(7)=(1.d0-r)*(1.d0+s)*(1.d0+t)/8.d0
        h(8)=(1.d0+r)*(1.d0+s)*(1.d0+t)/8.d0
        xint=sum(h*x(icon(:,ie)))
        yint=sum(h*y(icon(:,ie)))
        zint=sum(h*z(icon(:,ie)))
        write(iunit,'(3E14.6)') xint,yint,zint
        enddo
      enddo
    enddo
    write(iunit,'(a11,i10)')'POINT_DATA ',ne
    write(iunit,'(a)')'SCALARS sxx float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(1,ie)
      enddo
    write(iunit,'(a)')'SCALARS syy float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(2,ie)
      enddo
    write(iunit,'(a)')'SCALARS szz float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(3,ie)
      enddo
    write(iunit,'(a)')'SCALARS sxy float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(4,ie)
      enddo
    write(iunit,'(a)')'SCALARS sxz float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(5,ie)
      enddo
    write(iunit,'(a)')'SCALARS syz float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') stress(6,ie)
      enddo
    write(iunit,'(a)')'SCALARS p float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') p(ie)
      enddo
    write(iunit,'(a)')'SCALARS s2d float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') s2d(ie)
      enddo
    write(iunit,'(a)')'SCALARS fail float'
    write(iunit,'(a)')'LOOKUP_TABLE default'
      do ie=1,ne
      write(iunit,'(E14.6)') fail(ie)
      enddo
    
    deallocate (x,y,z)
    deallocate (p,s2d,fail)
    deallocate (icon)
    
    return
    end
    