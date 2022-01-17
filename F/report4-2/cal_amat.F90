!=======================================================================
      subroutine cal_amat_tet4_s(n,ne,kd,cny,num,coor,younglst,u,r)
!=======================================================================
!:: calculate r = K u by element-by-element method
!-----------------------------------------------------------------------
      implicit none
      ! in
      integer n,ne,kd
      integer cny(4,ne),num(ne)
      real coor(3,n),younglst(kd,2)
      real u(3,n)
      ! out
      real r(3,n)
      ! local      
      integer i,ie,in,cny1,cny2,cny3,cny4
      real xl21,xl22,xl23,xl31,xl32,xl33
      real xl41,xl42,xl43
      real ue11,ue21,ue31,ue12,ue22,ue32
      real ue13,ue23,ue33,ue14,ue24,ue34
      real young,rnyu,detj,detjt,d1,d2,d3
      real ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      real bmc1,bmc2,bmc3,bu1,bu2,bu3
      real dbu1,dbu2,dbu3,dbu4,dbu5,dbu6
      real bdbu11,bdbu21,bdbu31
      real bdbu12,bdbu22,bdbu32
      real bdbu13,bdbu23,bdbu33
      real bdbu14,bdbu24,bdbu34
!-----------------------------------------------------------------------
      do i=1,n
      r(1,i)=0.0
      r(2,i)=0.0
      r(3,i)=0.0
      enddo
      do ie=1,ne
      cny1=cny(1,ie)
      cny2=cny(2,ie)
      cny3=cny(3,ie)
      cny4=cny(4,ie)
      ue11=u(1,cny1)
      ue21=u(2,cny1)
      ue31=u(3,cny1)
      ue12=u(1,cny2)
      ue22=u(2,cny2)
      ue32=u(3,cny2)
      ue13=u(1,cny3)
      ue23=u(2,cny3)
      ue33=u(3,cny3)
      ue14=u(1,cny4)
      ue24=u(2,cny4)
      ue34=u(3,cny4)
      xl21=coor(1,cny2)-coor(1,cny1)
      xl22=coor(2,cny2)-coor(2,cny1)
      xl23=coor(3,cny2)-coor(3,cny1)
      xl31=coor(1,cny3)-coor(1,cny1)
      xl32=coor(2,cny3)-coor(2,cny1)
      xl33=coor(3,cny3)-coor(3,cny1)
      xl41=coor(1,cny4)-coor(1,cny1)
      xl42=coor(2,cny4)-coor(2,cny1)
      xl43=coor(3,cny4)-coor(3,cny1)
      in=num(ie)
      young=younglst(in,1)
      rnyu=younglst(in,2)
      d1=((1.0-rnyu)*young)/((1.0-2.0*rnyu)*(1.0+rnyu))
      d2=d1*rnyu/(1.0-rnyu)
      d3=young/(2.0*(1.0+rnyu))
      detj=-(xl23*xl32*xl41)+xl22*xl33*xl41+ &
       xl23*xl31*xl42-xl21*xl33*xl42- &
       xl22*xl31*xl43+xl21*xl32*xl43
      detjt=1./detj*0.16666666666666667  
      ij11=-(xl33*xl42)+xl32*xl43
      ij12=xl23*xl42-xl22*xl43
      ij13=-(xl23*xl32)+xl22*xl33
      ij21=xl33*xl41-xl31*xl43
      ij22=-(xl23*xl41)+xl21*xl43
      ij23=xl23*xl31-xl21*xl33
      ij31=-(xl32*xl41)+xl31*xl42
      ij32=xl22*xl41-xl21*xl42
      ij33=-(xl22*xl31)+xl21*xl32
      bmc1=-(ij11+ij12+ij13)
      bmc2=-(ij21+ij22+ij23)
      bmc3=-(ij31+ij32+ij33)
      bu1=bmc1*ue11 + ij11*ue12 + ij12*ue13 + ij13*ue14
      bu2=bmc2*ue21 + ij21*ue22 + ij22*ue23 + ij23*ue24
      bu3=bmc3*ue31 + ij31*ue32 + ij32*ue33 + ij33*ue34
      dbu4=(bmc2*ue11 + bmc1*ue21 + &
       ij21*ue12 + ij11*ue22 + &
       ij22*ue13 + ij12*ue23 + &
       ij23*ue14 + ij13*ue24)*d3*detjt
      dbu5=(bmc3*ue21 + bmc2*ue31 + &
       ij31*ue22 + ij21*ue32 + &
       ij32*ue23 + ij22*ue33 + &
       ij33*ue24 + ij23*ue34)*d3*detjt
      dbu6=(bmc3*ue11 + bmc1*ue31 + &
       ij31*ue12 + ij11*ue32 + &
       ij32*ue13 + ij12*ue33 + &
       ij33*ue14 + ij13*ue34)*d3*detjt
      dbu1=(d1*bu1 + d2*(bu2 + bu3))*detjt
      dbu2=(d1*bu2 + d2*(bu1 + bu3))*detjt
      dbu3=(d2*(bu1 + bu2) + d1*bu3)*detjt
      bdbu11=bmc1*dbu1 + bmc2*dbu4 + bmc3*dbu6
      bdbu21=bmc2*dbu2 + bmc1*dbu4 + bmc3*dbu5
      bdbu31=bmc3*dbu3 + bmc2*dbu5 + bmc1*dbu6
      bdbu12=ij11*dbu1 + ij21*dbu4 + ij31*dbu6
      bdbu22=ij21*dbu2 + ij11*dbu4 + ij31*dbu5
      bdbu32=ij31*dbu3 + ij21*dbu5 + ij11*dbu6
      bdbu13=ij12*dbu1 + ij22*dbu4 + ij32*dbu6
      bdbu23=ij22*dbu2 + ij12*dbu4 + ij32*dbu5
      bdbu33=ij32*dbu3 + ij22*dbu5 + ij12*dbu6
      bdbu14=ij13*dbu1 + ij23*dbu4 + ij33*dbu6
      bdbu24=ij23*dbu2 + ij13*dbu4 + ij33*dbu5
      bdbu34=ij33*dbu3 + ij23*dbu5 + ij13*dbu6
      r(1,cny1)=bdbu11+r(1,cny1)
      r(2,cny1)=bdbu21+r(2,cny1)
      r(3,cny1)=bdbu31+r(3,cny1)
      r(1,cny2)=bdbu12+r(1,cny2)
      r(2,cny2)=bdbu22+r(2,cny2)
      r(3,cny2)=bdbu32+r(3,cny2)
      r(1,cny3)=bdbu13+r(1,cny3)
      r(2,cny3)=bdbu23+r(2,cny3)
      r(3,cny3)=bdbu33+r(3,cny3)
      r(1,cny4)=bdbu14+r(1,cny4)
      r(2,cny4)=bdbu24+r(2,cny4)
      r(3,cny4)=bdbu34+r(3,cny4)
      enddo
      end
