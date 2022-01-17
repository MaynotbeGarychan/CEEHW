      program main
      implicit none
      integer n,ne,kd
      n=8192+1 ! number of nodes
      ne=n-1   ! number of elements
      kd=1     ! number of materials
      call main_allocate(n,ne,kd)
      end

      subroutine main_allocate(n,ne,kd)
      use omp_lib
      implicit none
      integer n,ne,kd
      integer cny(2,ne),num(ne)
      real*4 coor(n),younglst(kd)
      real*4 u(n),f(n)
      real*4 ds,norm1,norm2
      real*8 t_start,t_end
      integer i,n_iter

      write(*,'(a5,i4,a1,i4)') "n,ne ",n,",",ne

      ! set node coordinates
      ds=1.0e0/ne
      do i=1,n
        coor(i)=(i-1)*ds
      enddo

      ! set element connectivity and material number
      do i=1,ne
        cny(1,i)=i
        cny(2,i)=i+1
        num(i)=1
      enddo

      ! set material property
      younglst(1)=1.0e0

      ! set initial guess
      do i=1,n
        u(i)=coor(i)
      enddo

      ! benchmark
      n_iter=100000
      t_start=omp_get_wtime()
      do i=1,n_iter
      call cal_amat(n,ne,kd,cny,num,coor,younglst,u,f)
      enddo
      t_end=omp_get_wtime()
      write(*,'(a7,f10.8,a2,f10.8,a9,f11.8)') "0 took ",t_end-t_start,&
      " (",(n_iter*9.0*ne)/(t_end-t_start)/1.0e9," GFLOPS) ",f(1)

      ! check norm
      norm1=0.0e0
      norm2=0.0e0
      do i=1,n
        norm1=norm1+u(i)*u(i)
        norm2=norm2+f(i)*f(i)
      enddo
      write(*,'(a5,2E14.7)') "norm ",norm1,norm2

      return
      end
