module TopoStressesContext

    ! Context module for TopoStresses api
    ! should not be accessed or changed
    ! see API for name of routines and externally accessible variables
    
    implicit none
    
    integer :: nx,ny,nz,nn,nns,ne
    integer :: ibc
    double precision :: xl,yl,zl,rhog,exx,eyy,ym,pr,ks
    logical :: setup_has_been_run
    double precision, target, dimension(:), allocatable :: h
    double precision, target, dimension(:,:), allocatable :: stress
    double precision, dimension(:,:), pointer, contiguous :: h2
    
    contains
        
    !---------------------------------------------------------------
    
        subroutine SetUp(nxp,nyp,nzp)
    
        implicit none
    
        integer :: nxp,nyp,nzp

        nx=nxp
        ny=nyp
        nz=nzp

        nn=nx*ny*nz
        nns=nx*ny*nz
        ne=(nx-1)*(ny-1)*(nz-1)
    
        allocate (h(nns),stress(6,nn))
    
        h2(1:nx,1:ny)=>h
        
        setup_has_been_run = .true.
    
        return
    
        end subroutine SetUp
    
    !---------------------------------------------------------------
    
        subroutine Destroy()
    
        if (allocated(h)) deallocate(h)
        if (allocated(stress)) deallocate(stress)
    
        return
    
        end subroutine Destroy
    
    !---------------------------------------------------------------
    
        subroutine InitH (hp)
    
        double precision, intent(in), dimension(*) :: hp
    
        if (.not.setup_has_been_run) stop 'InitH - You need to run SetUp first'
    
        h = hp(1:nn)
    
        return
    
        end subroutine InitH
    
    !---------------------------------------------------------------
    
        subroutine View()
    
        write (*,*) 'TopoStressContext:'
        write (*,*) 'nx,ny,nz',nx,ny,nz
        write (*,*) 'nn,nnn,ne',nn,nns,ne
        write (*,*) 'xl,yl,zl',xl,yl,zl
        write (*,*) 'rhog,exx,eyy,ym,pr,ks',rhog,exx,eyy,ym,pr,ks
        write (*,*) 'ibc',ibc
        write (*,*) 'h',minval(h),sum(h)/nn,maxval(h)
        write (*,*) 'stress',minval(stress),sum(stress)/nn/6,maxval(stress)
    
        return
    
        end subroutine View
    
    !---------------------------------------------------------------
    
        subroutine SetXLYLZL (xxl,yyl,zzl)
    
        double precision, intent(in) :: xxl,yyl,zzl
    
        xl = xxl
        yl = yyl
        zl = zzl

        return
    
        end
    
    !---------------------------------------------------------------
    
        subroutine SetMechanicalParam (rrhog,eexx,eeyy,yym,ppr,kks)
    
        double precision, intent(in) :: rrhog,eexx,eeyy,yym,ppr,kks
        
        rhog=rrhog
        exx=eexx
        eyy=eeyy
        ym=yym
        pr=ppr
        ks=kks

        return
    
        end
        
    !---------------------------------------------------------------
    
        subroutine SetBC (jbc)
    
        integer, intent(in) :: jbc
    
        ibc = jbc
    
        return
    
        end subroutine SetBC
    
   !---------------------------------------------------------------
                 
        subroutine Copy_stress (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
        integer k,i

        do k=1,6
            do i=1,nn
            stressp((k-1)*nn+i)=stress(k,i)
            enddo
        enddo

        return
            
        end subroutine Copy_stress
            
   !---------------------------------------------------------------
                 
        subroutine Copy_sxx (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                        
        stressp(1:nn)=stress(1,1:nn)
                
        return
                
        end subroutine Copy_sxx
                
    !---------------------------------------------------------------
                 
        subroutine Copy_syy (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=stress(2,1:nn)
            
        return
            
        end subroutine Copy_syy

   !---------------------------------------------------------------
                 
        subroutine Copy_szz (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=stress(3,1:nn)
            
        return
            
        end subroutine Copy_szz


   !---------------------------------------------------------------
                 
        subroutine Copy_sxy (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=stress(4,1:nn)
            
        return
            
        end subroutine Copy_sxy

   !---------------------------------------------------------------
                 
        subroutine Copy_szx (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=stress(5,1:nn)
            
        return
            
        end subroutine Copy_szx

   !---------------------------------------------------------------
                 
        subroutine Copy_syz (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=stress(6,1:nn)
            
        return
            
        end subroutine Copy_syz

   !---------------------------------------------------------------
                 
        subroutine Copy_j2d (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
        
        stressp(1:nn)=(stress(1,1:nn)+stress(2,1:nn)+stress(3,1:nn))/3.d0
        stressp(1:nn)=sqrt(0.5d0*((stress(1,1:nn)-stressp(1:nn))**2+ &
                                    (stress(2,1:nn)-stressp(1:nn))**2+ &
                                    (stress(3,1:nn)-stressp(1:nn))**2)+ &
                                stress(4,1:nn)**2+stress(5,1:nn)**2+stress(6,1:nn)**2)
        print*,'j2d',minval(stressp(1:nn)),maxval(stressp(1:nn))
        return
            
        end subroutine Copy_j2d

   !---------------------------------------------------------------
                 
        subroutine Copy_press (stressp)
    
        double precision, intent(out), dimension(*) :: stressp
                    
        stressp(1:nn)=(stress(1,1:nn)+stress(2,1:nn)+stress(3,1:nn))/3.d0
            
        return
            
        end subroutine Copy_press

    !---------------------------------------------------------------

end module
    