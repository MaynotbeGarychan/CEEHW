#define NL 16 
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
      integer i,ie,ieo,in,cny1,cny2
      real*4 ue1(NL),ue2(NL),young(NL),detjt
      real*4 bdbu1(NL),bdbu2(NL),ds(NL)

      ! clear left hand side vector
      do i=1,n
        r(i)=0.0e0
      enddo

      if(mod(ne,NL)/=0) then
        write(*,'(a23,i2)') "ne cannot be divided by",NL
        stop
      endif

      do ieo=1,ne,NL
        ! load right hand side vector and compute components of B and D
        do ie=1,NL
          cny1=cny(1,ieo+ie-1)
          cny2=cny(2,ieo+ie-1)
          ue1(ie)=u(cny1)
          ue2(ie)=u(cny2)
          ds(ie)=coor(cny2)-coor(cny1)
          in=num(ieo+ie-1)
          young(ie)=younglst(in)
        enddo

        ! compute determinant and BDBu
        do ie=1,NL
          detjt=1.0e0/ds(ie)
          bdbu1(ie)=young(ie)*(ue1(ie)-ue2(ie))*detjt
          bdbu2(ie)=young(ie)*(ue2(ie)-ue1(ie))*detjt
        enddo

        ! add BDBu into left hand side vector
        do ie=1,NL
          cny1=cny(1,ieo+ie-1)
          cny2=cny(2,ieo+ie-1)
          r(cny1)=r(cny1)+bdbu1(ie)
          r(cny2)=r(cny2)+bdbu2(ie)
        enddo
      enddo
        
      return
      end
