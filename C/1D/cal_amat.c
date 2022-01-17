void cal_amat(int n,int ne,int kd,const int* cny,const int* num,const float* coor,const float* younglst,const float* const u,float* r){
/* calculate r = K u by element-by-element method */
  int i,ie,in,cny1,cny2;
  float ue1,ue2;
  float young,detjt;
  float bdbu1,bdbu2;
  float ds;

asm volatile("# clear left hand side vector");
  for( i=0; i < n; ++i )
    r[i]=0.0F;

asm volatile("# start ne loop");
  for( ie = 0; ie < ne; ++ie){
/* load right hand side vector */
    cny1=cny[2*ie];
    cny2=cny[2*ie+1];
    ue1=u[cny1];
    ue2=u[cny2];

/* compute components of B */
    ds=coor[cny2]-coor[cny1];

/* compute components of D */
    in=num[ie];
    young=younglst[in];

/* compute determinant */
    detjt=1.0F/ds;

/* compute BDBu */
    bdbu1=young*(ue1-ue2)*detjt;
    bdbu2=young*(ue2-ue1)*detjt;

/* add BDBu into left side vector */
    r[cny1]+=bdbu1;
    r[cny2]+=bdbu2;
  }
asm volatile("# end ne loop");
}
