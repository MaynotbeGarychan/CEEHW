#include <stdio.h>
#include <stdlib.h>
#define NL 16

void cal_amat(int n,int ne,int kd,const int* cny,const int* num,const float* coor,const float* younglst,const float* const u,float* r){
/* calculate r = K u by element-by-element method */
  int i,ie,in,cny1,cny2,ieo;
  float ue1[NL],ue2[NL];
  float young[NL],detjt;
  float bdbu1[NL],bdbu2[NL];
  float ds[NL];

asm volatile("# clear left hand side vector");
  for( i=0; i < n; ++i )
    r[i]=0.0F;

  if(ne%NL!=0){
    printf("ne cannot be divided by %d\n",ne);
    exit(0);
  }

  for( ieo = 0; ieo < ne; ieo += NL){
asm volatile("# load right hand side vector ");
    for( ie = 0; ie < NL; ++ie){
      cny1=cny[2*(ieo+ie)];
      cny2=cny[2*(ieo+ie)+1];
      ue1[ie]=u[cny1];
      ue2[ie]=u[cny2];
      ds[ie]=coor[cny2]-coor[cny1];
      in=num[ieo+ie];
      young[ie]=younglst[in];
    }
asm volatile("# compute determinant and BDBu");
    for( ie = 0; ie < NL; ++ie){
      detjt=1.0F/ds[ie];
      bdbu1[ie]=young[ie]*(ue1[ie]-ue2[ie])*detjt;
      bdbu2[ie]=young[ie]*(ue2[ie]-ue1[ie])*detjt;
    }
asm volatile("# add BDBu into left side vector");
    for( ie = 0; ie < NL; ++ie){
      cny1=cny[2*(ieo+ie)];
      cny2=cny[2*(ieo+ie)+1];
      r[cny1]+=bdbu1[ie];
      r[cny2]+=bdbu2[ie];
    }
  }
}
