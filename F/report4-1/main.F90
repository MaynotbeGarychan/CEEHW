      program main
      implicit none
      integer n,ne
      n=4096 ! number of nodes
      ne=8192   ! number of elements
      call main_allocate(n,ne)
      end

      subroutine main_allocate(n,ne)
      use omp_lib
      implicit none
      integer n,ne
      integer cny(ne)
      real*4 u(n),f(n),mat(ne)
      real*4 norm1,norm2
      real*8 t_start,t_end
      integer i,ie,n_iter

      write(*,'(a5,i4,a1,i4)') "n,ne ",n,",",ne

      do ie=1,ne
        cny(ie)=mod(ie-1,n)+1
        mat(ie)=1.0e0
      enddo
      do i=1,n
        u(i)=i
      enddo

      ! benchmark
      n_iter=100000
      t_start=omp_get_wtime()
      do i=1,n_iter
      call func(n,ne,cny,mat,u,f)
      enddo
      t_end=omp_get_wtime()
      write(*,'(a7,f10.5,a1,f10.5)') "0 took ",t_end-t_start," ",f(n)

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
