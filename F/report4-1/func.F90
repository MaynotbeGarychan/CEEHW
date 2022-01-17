      subroutine func(n,ne,cny,mat,u,r)
      implicit none
      ! input
      integer n,ne
      integer cny(ne)
      real*4 u(n),mat(ne)
      ! output
      real*4 r(n)
      ! temporal
      integer i,ie,cnytmp
      real*4 utmp,mattmp,rtmp

      do i=1,n
        r(i)=0.0
      enddo

      do ie=1,ne
        cnytmp=cny(ie)
        utmp=u(cnytmp)
        mattmp=mat(ie)
        rtmp=((1.0/(utmp+1.0))**3)
        rtmp=rtmp+((2.0/(utmp+2.0))**3)
        rtmp=rtmp+((3.0/(utmp+3.0))**3)
        rtmp=rtmp+((4.0/(utmp+4.0))**3)
        rtmp=rtmp/mattmp
        r(cnytmp)=r(cnytmp)+rtmp
      enddo
        
      return
      end

