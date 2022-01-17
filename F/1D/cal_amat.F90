      subroutine cal_amat(n,ne,kd,cny,num,coor,younglst,u,r)
      implicit none
      ! input
      integer n,ne,kd
      integer cny(2,ne),num(ne)
      real*4 coor(n),younglst(kd)
      real*4 u(n)
      ! output
      real*4 r(n)
      ! temporal
      integer i,ie,in,cny1,cny2
      real*4 ue1,ue2,young,detjt
      real*4 bdbu1,bdbu2,ds

      ! clear left hand side vector
      do i=1,n
        r(i)=0.0e0
      enddo

      do ie=1,ne
        ! load right hand side vector
        cny1=cny(1,ie)
        cny2=cny(2,ie)
        ue1=u(cny1)
        ue2=u(cny2)

        ! compute components of B
        ds=coor(cny2)-coor(cny1)

        ! compute components of D
        in=num(ie)
        young=younglst(in)

        ! compute determinant
        detjt=1.0e0/ds

        ! compute BDBu
        bdbu1=young*(ue1-ue2)*detjt
        bdbu2=young*(ue2-ue1)*detjt

        ! add BDBu into left hand side vector
        r(cny1)=r(cny1)+bdbu1
        r(cny2)=r(cny2)+bdbu2
      enddo
        
      return
      end

