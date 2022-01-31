      subroutine solve_skyline_renumber_n &
                    (ndof,n,ael,f_true,x,icon_true,nnode,nelem,iop)

! subroutine to solve by Cholesky factorization a symmetrical positive
! definite matrix. The matrix is provided as elemental ndofxndof matrices
! ael(k,j,ielem). The first step is to store it in a skyline fashion.
! The second step is the Cholesky factorization.
! The third step is the double back substitution.

! in input ndof: number of d.o.f. per node
! n is the number of nodes per element
! ael: elemental matrices (nelem of them of dimension ndof*n, ndof*n)
!      note that only the upper half of the ael's are used because
!      the ael's are SPD (ael(i,j,ie) are defined for i>=j)
! f: right hand-side vector
! x: solution vector
! icon: connectivity array (n nodes per element)
! nnode: number of nodes
! nelem: number of elements
! a: global matrix sotred in skyline format
! l, ll, m: working integer arrays of dimension ndof*nnode
! namax: maximum dimension of a in calling routine
! iop: operation to be performed (-1 means assemblage only
!                                 0 means factorization only
!                                 1 means back-subtitution only
!                                 2 means ordering, formation, factorization and
!                                         back-substitution
!                                 3 means formation, factorization and
!                                         back-substitution

      implicit none!real*8 (a-h,o-z)

      integer ndof,n,nnode,nelem
      double precision  ael(n*ndof,n*ndof,nelem)
      double precision  x(ndof*nnode),f_true(ndof*nnode)
      integer icon_true(n,nelem)

      double precision,dimension(:),allocatable::a,f
      integer,dimension(:,:),allocatable::icon,ixyz,mixyz
      integer,dimension(:),allocatable::l,ll,m,sort,jcon,adj,xadj,iw

      integer p,pmin,namax
      integer iop,ireorder,i,j,k,ie,ii,id,ka,ia,kb,ib
      integer idf,iip,jd,mi,ind,mkk,mk,k1,klk,mjj,j1,mkj
      integer ijd,mj,jlj,ip,ij,ijp,kllk,mpk
      double precision akj,akk,xk

      allocate (f(ndof*nnode))
      allocate (icon(n,nelem))
      allocate (l(nnode*ndof),ll(nnode*ndof),m(nnode*ndof))
      allocate (ixyz(ndof,n))
      allocate (mixyz(ndof,n))
      allocate (sort(nnode),jcon(nelem+1),adj(nelem*96))
      allocate (xadj(nnode+1),iw(3*nnode+1))

      if (iop.eq.0) goto 111
      if (iop.eq.1) goto 222
      if (iop.eq.3) goto 555

! first we have to find the two arrays l and m
!    m: where the diagonal elements of the matrix will be stored in the
!       reduced skyline format
!    l: width of each row in the skyline format (note l(1)=0)

! re-ordering

      ireorder=1

        if (ireorder.eq.0) then
          do i=1,nnode
          sort(i)=i
          enddo
        else
          do ie=1,nelem
            do k=1,n
            icon(k,ie)=icon_true(k,ie)
            enddo
          enddo
        call re_number_sloane (icon,nelem,n, &
                               nnode,jcon,adj,xadj,iw,sort)
        endif

        do ie=1,nelem
          do k=1,n
          icon(k,ie)=sort(icon_true(k,ie))
          enddo
        enddo

! initialization
!    note that the x-degrees of freedom can have l=0
!    whereas the y-degrees of freedom have l=1 (as a minimum) 
!       and the z-degrees of freedom have l=2 (as aminimum) because
!       every x-dof is connected to a y- and z-dof

      do j=1,ndof

        do i=j,nnode*ndof,ndof
        l(i)=j-1
        enddo

      enddo

! we first compute the l's and ll's using a first pass through the icon array
! we only have to store the a(k,j) for which k>j

        do ie=1,nelem
          do k=1,n
          ii=icon(k,ie)
            do id=1,ndof
            ixyz(id,k)=(ii-1)*ndof+id
            enddo
          enddo
          do ka=1,n
          ia=icon(ka,ie)
            do kb=ka+1,n
            ib=icon(kb,ie)
            idf=(ia-ib)*ndof
              if (idf.gt.0) then
                do id=1,ndof
                l(ixyz(id,ka))=max0(l(ixyz(id,ka)),idf+id-1)
                enddo
              else
                do id=1,ndof
                l(ixyz(id,kb))=max0(l(ixyz(id,kb)),-idf+id-1)
                enddo
              endif
            enddo
          enddo
        enddo

      ll=0
        do i=nnode*ndof,1,-1
          do j=i-l(i),i-1
          ll(j)=max0(ll(j),i-j)
          enddo
        enddo

! from the l's we can easily compute the m's

      m(1)=1
        do i=2,nnode*ndof
        m(i)=m(i-1)+l(i)+1
        enddo

! check if we have reserved enough storage for the vector a

      namax=m(nnode*ndof)
      allocate (a(namax))

555   continue

! initialize vector a

        do i=1,namax
        a(i)=0.
        enddo

! we now distribute the elemental nxn matrices into the big a vector

        do ie=1,nelem
          do k=1,n
          ii=icon(k,ie)
            do id=1,ndof
            iip=(ii-1)*ndof+id
            ixyz(id,k)=iip
            mixyz(id,k)=m(iip)
            enddo
          enddo
          do id=1,ndof
            do jd=1,id
            ijd=jd-id
              do k=1,n
              mi=mixyz(id,k)+ijd
              a(mi)=a(mi)+ael(ndof*(k-1)+id,ndof*(k-1)+jd,ie)
              enddo
            enddo
          enddo
          do ka=1,n
          ia=icon(ka,ie)
            do kb=ka+1,n
            ib=icon(kb,ie)
              if (ib.gt.ia) then
                do id=1,ndof
                  do jd=1,ndof
                  ind=mixyz(id,kb)-ixyz(id,kb)+ixyz(jd,ka)
                  a(ind)=a(ind)+ael(ndof*(kb-1)+id,ndof*(ka-1)+jd,ie)
                  enddo
                enddo
              else
                do id=1,ndof
                  do jd=1,ndof
                  ind=mixyz(jd,ka)-ixyz(jd,ka)+ixyz(id,kb)
                  a(ind)=a(ind)+ael(ndof*(kb-1)+id,ndof*(ka-1)+jd,ie)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo

      if (iop.eq.-1) then
      deallocate (f)
      deallocate (icon)
      deallocate (l,ll,m)
      deallocate (ixyz)
      deallocate (mixyz)
      deallocate (sort,jcon,adj)
      deallocate (xadj,iw)
      deallocate (a)
      return
      endif

111   continue

! Cholesky factorization

        do k=1,nnode*ndof
        mkk=m(k)
        mk=mkk-k
        k1=k-1
        klk=k-l(k)
          do j=klk,k1
          mjj=m(j)
          mj=mjj-j
          j1=j-1
          jlj=j-l(j)
          pmin=max0(klk,jlj)
          mkj=mk+j
          akj=a(mkj)
            do p=pmin,j1
            akj=akj-a(mk+p)*a(mj+p)
            enddo
          a(mkj)=akj/a(mjj)
          enddo
        akk=a(mkk)
          do p=klk,k1
          akk=akk-a(mk+p)**2
          enddo
          if (akk.le.0.) then
          print*,k,akk
          print*,'negative pivot at ',k
          stop
          endif
        a(mkk)=sqrt(akk)
        enddo

! end of factorization

      if (iop.eq.0) then
      deallocate (f)
      deallocate (icon)
      deallocate (l,ll,m)
      deallocate (ixyz)
      deallocate (mixyz)
      deallocate (sort,jcon,adj)
      deallocate (xadj,iw)
      deallocate (a)
      return
      endif

222   continue

! re-ordering

        do i=1,nnode
        ip=sort(i)
          do j=1,ndof
          ij=(i-1)*ndof+j
          ijp=(ip-1)*ndof+j
          f(ijp)=f_true(ij)
          enddo
        enddo

! first backsubstitution

        do k=1,nnode*ndof
        xk=f(k)
        mkk=m(k)
        mk=mkk-k
        k1=k-1
        klk=k-l(k)
          do p=klk,k1
          xk=xk-a(mk+p)*x(p)
          enddo
        x(k)=xk/a(mkk)
        enddo

! second backsubstitution

        do k=nnode*ndof,1,-1
        xk=x(k)
        mkk=m(k)
        kllk=k+ll(k)
          do p=k+1,kllk
          if (p-k.le.l(p)) then
          mpk=m(p)-p+k
          xk=xk-a(mpk)*x(p)
          endif
          enddo
        x(k)=xk/a(mkk)
        enddo

! re-re-ordering

        do i=1,nnode
        ip=sort(i)
          do j=1,ndof
          ij=(i-1)*ndof+j
          ijp=(ip-1)*ndof+j
          f(ij)=x(ijp)
          enddo
        enddo

        do i=1,nnode*ndof
        x(i)=f(i)
        enddo

      deallocate (f)
      deallocate (icon)
      deallocate (l,ll,m)
      deallocate (ixyz)
      deallocate (mixyz)
      deallocate (sort,jcon,adj)
      deallocate (xadj,iw)
      deallocate (a)

      return
      end

!---
      subroutine re_number_sloane (icon,nelem,n, &
                                   nnode,jcon,adj,xadj,iw,sort)

! routine to renumber nodes to minimize bandwisth based on Sloane's
! subroutines

! in input icon, nelem, nnode

! in output icon_new,nelem_new

! working arrays : jcon(nelem+1),adj(nelem*96),xadj(nnode+1),iw(3*nnode+1)

      implicit none

      integer   nelem,n,nnode,ie
      integer   icon(n,nelem),jcon(nelem+1)
      integer   adj(nelem*96),xadj(nnode+1),iw(3*nnode+1)
      integer   sort(nnode),oldpro,newpro

        do ie=1,nelem
        jcon(ie)=1+(ie-1)*n
        enddo
      jcon(nelem+1)=jcon(nelem)+1

      call graph_sloan (nnode,nelem,n*nelem,icon,jcon, &
                        nelem*96,adj,xadj)
      if (xadj(nnode+1)-1.gt.nelem*96) stop 'must increase size of adj'

      call label_sloan (nnode,xadj(nnode+1)-1,adj,xadj, &
                        sort,iw,oldpro,newpro)

      return
      end


