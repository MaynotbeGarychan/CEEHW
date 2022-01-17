      implicit none
      integer nblock,n,ne,kd
! make nblock x nblock x nblock structured mesh
      nblock=20
      n=(nblock+1)**3
      ne=6*(nblock**3)
      kd=1
      call main_allocate(nblock,n,ne,kd)
      end
      
      subroutine main_allocate(nblock,n,ne,kd)
      use omp_lib
      implicit none
      integer nblock,n,ne,kd
      integer i,j,k,itmp,ie
      integer cny(4,ne),num(ne)
      real coor(3,n),younglst(kd,2)
      real u(3,n),f(3,n)
      real ds,norm1,norm2
      real*8 t1,t2
      integer cnylocal(4,6),localn(8)
      integer niter
       
       write(*,*) 'nblock,n,ne',nblock,n,ne
       
! set node coordinates
      ds=1.0
      itmp=0
      do i=1,nblock+1
      do j=1,nblock+1
      do k=1,nblock+1
      itmp=itmp+1
      coor(1,itmp)=(i-1)*ds
      coor(2,itmp)=(j-1)*ds
      coor(3,itmp)=(k-1)*ds
      enddo
      enddo
      enddo

! set element connectivity and material number
      cnylocal(1,1)=2
      cnylocal(2,1)=7
      cnylocal(3,1)=3
      cnylocal(4,1)=5
      cnylocal(1,2)=4
      cnylocal(2,2)=6
      cnylocal(3,2)=8
      cnylocal(4,2)=7
      cnylocal(1,3)=2
      cnylocal(2,3)=7
      cnylocal(3,3)=4
      cnylocal(4,3)=3
      cnylocal(1,4)=3
      cnylocal(2,4)=2
      cnylocal(3,4)=5
      cnylocal(4,4)=1
      cnylocal(1,5)=4
      cnylocal(2,5)=6
      cnylocal(3,5)=7
      cnylocal(4,5)=2
      cnylocal(1,6)=6
      cnylocal(2,6)=7
      cnylocal(3,6)=2
      cnylocal(4,6)=5

      itmp=0
      do i=1,nblock
      do j=1,nblock
      do k=1,nblock     
      localn(1)=k+    (nblock+1)*(j-1)+(nblock+1)*(nblock+1)*(i-1)
      localn(2)=k+1+(nblock+1)*(j-1)+(nblock+1)*(nblock+1)*(i-1)
      localn(3)=k+    (nblock+1)*(j   )+(nblock+1)*(nblock+1)*(i-1)
      localn(4)=k+1+(nblock+1)*(j   )+(nblock+1)*(nblock+1)*(i-1)
      localn(5)=k+    (nblock+1)*(j-1)+(nblock+1)*(nblock+1)*i
      localn(6)=k+1+(nblock+1)*(j-1)+(nblock+1)*(nblock+1)*i
      localn(7)=k+    (nblock+1)*(j   )+(nblock+1)*(nblock+1)*i
      localn(8)=k+1+(nblock+1)*(j   )+(nblock+1)*(nblock+1)*i
      do ie=1,6
      itmp=itmp+1
      num(itmp)=1      
      cny(1,itmp)=localn(cnylocal(1,ie))
      cny(2,itmp)=localn(cnylocal(2,ie))
      cny(3,itmp)=localn(cnylocal(3,ie))
      cny(4,itmp)=localn(cnylocal(4,ie))
      enddo
      enddo
      enddo
      enddo
      if(itmp.ne.ne)then
      write(*,*) 'error in ne'
      endif
      
! set material property
      younglst(1,1)=1.0
      younglst(1,2)=0.1
 
! set initial guess
      do i=1,n
      u(1,i)=coor(1,i)
      u(2,i)=coor(2,i)
      u(3,i)=coor(3,i)
      enddo
 
! benchmark
      niter=500
      t1=omp_get_wtime()
      do i=1,niter
      call cal_amat_tet4_s(n,ne,kd,cny,num,coor,younglst,u,f)
      enddo
      t2=omp_get_wtime()
      write(*,'(a7,f10.8,a2,f10.6,a9,f11.8)') "0 took ",t2-t1, &
      " (",(niter*246.0*ne)/(t2-t1)/1.0e9," GFLOPS) ",f(1,1)
      
      norm1=0.0
      norm2=0.0
      do i=1,n
      norm1=norm1+f(1,i)**2+f(2,i)**2+f(3,i)**2;
      norm2=norm2+u(1,i)**2+u(2,i)**2+u(3,i)**2;
      enddo
      write(*,*) 'norm',norm1,norm2
            
      end
