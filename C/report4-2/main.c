#include <stdio.h>
#include <omp.h>
#include "fem.h"

int main(){
  int nblock,n,ne,kd;
/* make nblock * nblock * nblock structured mesh */
  nblock=20;
  n=(nblock+1)*(nblock+1)*(nblock+1); /* number of nodes */
  ne=6*nblock*nblock*nblock; /* number of elements */
  kd=1; /* number of materials */
  main_allocate(nblock,n,ne,kd);
  return 0;
}     

void main_allocate(int nblock, int n,int ne,int kd){
  int i,j,k,itmp,ie,niter;
  int cny[4*ne],num[ne];
  float coor[3*n],younglst[2*kd];
  float u[3*n],f[3*n];
  float ds,norm1,norm2;
  double t1,t2;
  int cnylocal[4][6],localn[8];

  printf("nblock,n,ne %d %d %d\n",nblock,n,ne);

/* set node coordinates */
  ds=1.0;
  itmp=0;
  for( i=0; i < nblock+1; ++i ){
    for( j=0; j < nblock+1; ++j ){
      for( k=0; k < nblock+1; ++k ){
        coor[3*itmp]=i*ds; /* x coordinate */
        coor[3*itmp+1]=j*ds; /* y coordinate */
        coor[3*itmp+2]=k*ds; /* z coordinate */
        itmp=itmp+1;
      }
    }
  }

/* set element connectivity and material number */
  cnylocal[0][0]=1;
  cnylocal[1][0]=6;
  cnylocal[2][0]=2;
  cnylocal[3][0]=4;
  cnylocal[0][1]=3;
  cnylocal[1][1]=5;
  cnylocal[2][1]=7;
  cnylocal[3][1]=6;
  cnylocal[0][2]=1;
  cnylocal[1][2]=6;
  cnylocal[2][2]=3;
  cnylocal[3][2]=2;
  cnylocal[0][3]=2;
  cnylocal[1][3]=1;
  cnylocal[2][3]=4;
  cnylocal[3][3]=0;
  cnylocal[0][4]=3;
  cnylocal[1][4]=5;
  cnylocal[2][4]=6;
  cnylocal[3][4]=1;
  cnylocal[0][5]=5;
  cnylocal[1][5]=6;
  cnylocal[2][5]=1;
  cnylocal[3][5]=4;

  itmp=0;
  for( i=0; i < nblock; ++i ){
    for( j=0; j < nblock; ++j ){
      for( k=0; k < nblock; ++k ){
        localn[0]=k+  (nblock+1)*j    +(nblock+1)*(nblock+1)*i;
        localn[1]=k+1+(nblock+1)*j    +(nblock+1)*(nblock+1)*i;
        localn[2]=k+  (nblock+1)*(j+1)+(nblock+1)*(nblock+1)*i;
        localn[3]=k+1+(nblock+1)*(j+1)+(nblock+1)*(nblock+1)*i;
        localn[4]=k+  (nblock+1)*j    +(nblock+1)*(nblock+1)*(i+1);
        localn[5]=k+1+(nblock+1)*j    +(nblock+1)*(nblock+1)*(i+1);
        localn[6]=k+  (nblock+1)*(j+1)+(nblock+1)*(nblock+1)*(i+1);
        localn[7]=k+1+(nblock+1)*(j+1)+(nblock+1)*(nblock+1)*(i+1);
        for( ie=0; ie < 6; ++ie ){
          num[itmp]=0;
          cny[4*itmp+0]=localn[cnylocal[0][ie]];
          cny[4*itmp+1]=localn[cnylocal[1][ie]];
          cny[4*itmp+2]=localn[cnylocal[2][ie]];
          cny[4*itmp+3]=localn[cnylocal[3][ie]];
          itmp=itmp+1;
        }
      }
    }
  }
  if(itmp!=ne) printf("error in ne\n");

/* set material property */
  younglst[0]=1.0;
  younglst[1]=0.1;

/* set initial guess */
  for( i=0; i < 3*n; ++i)
    u[i] = coor[i];

/* benchmark */
  niter=500;
  t1=omp_get_wtime();
  for( i=0; i< niter; ++i )
    cal_amat(n,ne,kd,cny,num,coor,younglst,u,f);
  t2=omp_get_wtime();
  printf("0 took %f (%f GFLOPS) %f\n",t2-t1,(niter*246.0*ne)/(t2-t1)/1.0E9,f[0]);

  norm1=0.0;
  norm2=0.0;
  for( i=0; i < 3*n; ++i )
  {
    norm1+=f[i]*f[i];
    norm2+=u[i]*u[i];
  }
  printf("norm %e %e\n",norm1,norm2);
  
}
