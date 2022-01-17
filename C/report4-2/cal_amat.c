void cal_amat(int n,int ne,int kd,const int* restrict cny,const int* restrict num,const float* restrict coor,const float* restrict younglst,const float* const restrict u,float* restrict r){
/* calculate r = K u by element-by-element method */
  int i,ie,in,cny1,cny2,cny3,cny4;
  float xl21,xl22,xl23,xl31,xl32,xl33;
  float xl41,xl42,xl43;
  float ue11,ue21,ue31,ue12,ue22,ue32;
  float ue13,ue23,ue33,ue14,ue24,ue34;
  float young,rnyu,detj,detjt,d1,d2,d3;
  float ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33;
  float bmc1,bmc2,bmc3,bu1,bu2,bu3,bu4,bu5,bu6;
  float dbu1,dbu2,dbu3,dbu4,dbu5,dbu6;
  float bdbu11,bdbu21,bdbu31;
  float bdbu12,bdbu22,bdbu32;
  float bdbu13,bdbu23,bdbu33;
  float bdbu14,bdbu24,bdbu34;

/* clear left hand side vector */
  for( i=0; i < 3*n; ++i )
    r[i]=0.0F;

  for( ie = 0; ie < ne; ++ie){
/* load right hand side vector */
    cny1=cny[4*ie];
    cny2=cny[4*ie+1];
    cny3=cny[4*ie+2];
    cny4=cny[4*ie+3];
    ue11=u[3*cny1];
    ue21=u[3*cny1+1];
    ue31=u[3*cny1+2];
    ue12=u[3*cny2];
    ue22=u[3*cny2+1];
    ue32=u[3*cny2+2];
    ue13=u[3*cny3];
    ue23=u[3*cny3+1];
    ue33=u[3*cny3+2];
    ue14=u[3*cny4];
    ue24=u[3*cny4+1];
    ue34=u[3*cny4+2];

/* compute components of B */
    xl21=coor[3*cny2]-coor[3*cny1];
    xl22=coor[3*cny2+1]-coor[3*cny1+1];
    xl23=coor[3*cny2+2]-coor[3*cny1+2];
    xl31=coor[3*cny3]-coor[3*cny1];
    xl32=coor[3*cny3+1]-coor[3*cny1+1];
    xl33=coor[3*cny3+2]-coor[3*cny1+2];
    xl41=coor[3*cny4]-coor[3*cny1];
    xl42=coor[3*cny4+1]-coor[3*cny1+1];
    xl43=coor[3*cny4+2]-coor[3*cny1+2];
    ij11=-(xl33*xl42)+xl32*xl43;
    ij12=xl23*xl42-xl22*xl43;
    ij13=-(xl23*xl32)+xl22*xl33;
    ij21=xl33*xl41-xl31*xl43;
    ij22=-(xl23*xl41)+xl21*xl43;
    ij23=xl23*xl31-xl21*xl33;
    ij31=-(xl32*xl41)+xl31*xl42;
    ij32=xl22*xl41-xl21*xl42;
    ij33=-(xl22*xl31)+xl21*xl32;
    bmc1=-(ij11+ij12+ij13);
    bmc2=-(ij21+ij22+ij23);
    bmc3=-(ij31+ij32+ij33);

/* compute components of D */
    in=num[ie];
    young=younglst[2*in];
    rnyu=younglst[2*in+1];
    d1=((1.0F-rnyu)*young)/((1.0F-2.0F*rnyu)*(1.0F+rnyu));
    d2=d1*rnyu/(1.0F-rnyu);
    d3=young/(2.0F*(1.0F+rnyu));

/* compute determinant */
    detj=-(xl23*xl32*xl41)+xl22*xl33*xl41+xl23*xl31*xl42-xl21*xl33*xl42-xl22*xl31*xl43+xl21*xl32*xl43;
    detjt=1.0F/detj*0.16666666666666667F;

/* compute Bu */
    bu1=bmc1*ue11 + ij11*ue12 + ij12*ue13 + ij13*ue14;
    bu2=bmc2*ue21 + ij21*ue22 + ij22*ue23 + ij23*ue24;
    bu3=bmc3*ue31 + ij31*ue32 + ij32*ue33 + ij33*ue34;
    bu4=bmc2*ue11 + bmc1*ue21 + ij21*ue12 + ij11*ue22 + 
      ij22*ue13 + ij12*ue23 + ij23*ue14 + ij13*ue24;
    bu5=bmc3*ue21 + bmc2*ue31 + ij31*ue22 + ij21*ue32 + 
      ij32*ue23 + ij22*ue33 + ij33*ue24 + ij23*ue34;
    bu6=bmc3*ue11 + bmc1*ue31 + ij31*ue12 + ij11*ue32 + 
      ij32*ue13 + ij12*ue33 + ij33*ue14 + ij13*ue34;

/* compute DBu */
    dbu1=(d1*bu1 + d2*(bu2 + bu3))*detjt;
    dbu2=(d1*bu2 + d2*(bu1 + bu3))*detjt;
    dbu3=(d2*(bu1 + bu2) + d1*bu3)*detjt;
    dbu4=bu4*d3*detjt;
    dbu5=bu5*d3*detjt;
    dbu6=bu6*d3*detjt;

/* compute BDBu */
    bdbu11=bmc1*dbu1 + bmc2*dbu4 + bmc3*dbu6;
    bdbu21=bmc2*dbu2 + bmc1*dbu4 + bmc3*dbu5;
    bdbu31=bmc3*dbu3 + bmc2*dbu5 + bmc1*dbu6;
    bdbu12=ij11*dbu1 + ij21*dbu4 + ij31*dbu6;
    bdbu22=ij21*dbu2 + ij11*dbu4 + ij31*dbu5;
    bdbu32=ij31*dbu3 + ij21*dbu5 + ij11*dbu6;
    bdbu13=ij12*dbu1 + ij22*dbu4 + ij32*dbu6;
    bdbu23=ij22*dbu2 + ij12*dbu4 + ij32*dbu5;
    bdbu33=ij32*dbu3 + ij22*dbu5 + ij12*dbu6;
    bdbu14=ij13*dbu1 + ij23*dbu4 + ij33*dbu6;
    bdbu24=ij23*dbu2 + ij13*dbu4 + ij33*dbu5;
    bdbu34=ij33*dbu3 + ij23*dbu5 + ij13*dbu6;

/* add BDBu into left side vector */
    r[3*cny1]+=bdbu11;
    r[3*cny1+1]+=bdbu21;
    r[3*cny1+2]+=bdbu31;
    r[3*cny2]+=bdbu12;
    r[3*cny2+1]+=bdbu22;
    r[3*cny2+2]+=bdbu32;
    r[3*cny3]+=bdbu13;
    r[3*cny3+1]+=bdbu23;
    r[3*cny3+2]+=bdbu33;
    r[3*cny4]+=bdbu14;
    r[3*cny4+1]+=bdbu24;
    r[3*cny4+2]+=bdbu34;
  }
}