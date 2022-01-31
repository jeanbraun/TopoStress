! TopoStressAPI

! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Setup(nnx,nny,nnz)

use TopoStressesContext
        
implicit none
        
integer, intent(in) :: nnx,nny,nnz
        
call SetUp (nnx,nny,nnz)
        
return
        
end subroutine TopoStresses_Setup
    
! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Set_xlylzl(xlp,ylp,zlp)

use TopoStressesContext
            
implicit none
            
double precision, intent(in) :: xlp,ylp,zlp
            
call SetXLYLZL (xlp,ylp,zlp)
            
return
            
end subroutine TopoStresses_Set_xlylzl
        
! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Set_Mechanical_Properties(rrhog,eexx,eeyy,yym,ppr,kks)

use TopoStressesContext
                
implicit none
                
double precision, intent(in) :: rrhog,eexx,eeyy,yym,ppr,kks
                
call SetMechanicalParam (rrhog,eexx,eeyy,yym,ppr,kks)
                
return
                
end subroutine TopoStresses_Set_Mechanical_Properties

! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Set_BC(jbc)

use TopoStressesContext
                    
implicit none
                    
integer, intent(in) :: jbc
                    
call SetBC (jbc)
                    
return
                    
end subroutine TopoStresses_Set_BC
                
            
! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Set_H(hp)

use TopoStressesContext
            
implicit none
            
double precision, intent(inout), dimension(*) :: hp
            
call InitH (hp)
            
return
            
end subroutine TopoStresses_Set_H

! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Execute()

use TopoStressesContext
                
implicit none
                                
call TopoStresses ()!h,nx,ny,nz,xl,yl,zl,rhog,exx,eyy,ym,pr,ks,ibc,stress)

return
                
end subroutine TopoStresses_Execute
        
! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_stress(stressp)

use TopoStressesContext
                    
implicit none
                    
double precision, intent(inout), dimension(*) :: stressp

!print*,'b',minval(stressp(1:6,1:nn)),maxval(stressp(1:6,1:nn))
call Copy_stress (stressp)
!print*,'a',minval(stressp(1:6,1:nn)),maxval(stressp(1:6,1:nn))

return
                    
end subroutine TopoStresses_Copy_stress

! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_sxx(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_sxx (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_sxx
    
! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_syy(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_syy (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_syy
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_szz(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_szz (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_szz
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_sxy(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_sxy (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_sxy
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_syz(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_syz (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_syz
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_szx(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_szx (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_szx
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_j2d(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_j2d (stressp)
    print*,'copy',minval(stressp(1:nn)),maxval(stressp(1:nn))
                        
    return
                        
    end subroutine TopoStresses_Copy_j2d
    
    ! -----------------------------------------------------------------------------------------

subroutine TopoStresses_Copy_press(stressp)

    use TopoStressesContext
                        
    implicit none
                        
    double precision, intent(inout), dimension(*) :: stressp
                        
    call Copy_press (stressp)
                        
    return
                        
    end subroutine TopoStresses_Copy_press
    
    !--------------------------------------------------------------------------

subroutine TopoStresses_Destroy()

use TopoStressesContext
    
implicit none
    
call Destroy()
    
return
    
end subroutine TopoStresses_Destroy
    
       
