#include <stdio.h>
#include <omp.h>
#include "fem.h"

int main(){
  int n,ne,kd;
  n=8192+1; /* number of nodes */
  ne=n-1; /* number of elements */
  kd=1; /* number of materials */
  main_allocate(n,ne,kd);
  return 0;
}

void main_allocate(int n,int ne,int kd){
  int i,niter;
  int cny[2*ne],num[ne];
  float coor[n],younglst[kd];
  float u[n],f[n];
  float ds,norm1,norm2;
  double t1,t2;

  printf("n,ne %d %d\n",n,ne);

/* set node coordinates */
  ds=1.0/ne;
  for( i=0; i < n; ++i )
    coor[i]=i*ds; /* coordinate */

/* set element connectivity and material number */
  for( i=0; i < ne; ++i ){
    cny[2*i]=i;
    cny[2*i+1]=i+1;
    num[i]=0;
  }

/* set material property */
  younglst[0]=1.0;

/* set initial guess */
  for( i=0; i < n; ++i)
    u[i] = coor[i];

/* benchmark */
  niter=100000;
  t1=omp_get_wtime();
  for( i=0; i< niter; ++i )
    cal_amat(n,ne,kd,cny,num,coor,younglst,u,f);
  t2=omp_get_wtime();
  printf("0 took %f (%f GFLOPS) %f\n",t2-t1,(niter*9.0*ne)/(t2-t1)/1.0E9,f[0]);

/* check norm */
  norm1=0.0;
  norm2=0.0;
  for( i=0; i < n; ++i )
  {
    norm1+=u[i]*u[i];
    norm2+=f[i]*f[i];
  }
  printf("norm %e %e\n",norm1,norm2);
  
}
