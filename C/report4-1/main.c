#include <stdio.h>
#include <omp.h>
#include "header.h"

int main(){
  int n,ne;
  n=4096; /* number of nodes */
  ne=8192; /* number of elements */
  main_allocate(n,ne);
  return 0;
}

void main_allocate(int n,int ne){
  int cny[ne];
  float u[n],f[n],mat[ne];
  float norm1,norm2;
  double t1,t2;
  int i,ie,niter;

  printf("n,ne %d %d\n",n,ne);

  for( ie=0; ie < ne; ++ie ){
    cny[ie]=ie%n;
    mat[ie]=1.0;
  }
  for( i=0; i < n; ++i)
    u[i] = i+1;

/* benchmark */
  niter=100000;
  t1=omp_get_wtime();
  for( i=0; i< niter; ++i )
    func(n,ne,cny,mat,u,f);
  t2=omp_get_wtime();
  printf("0 took %f %f\n",t2-t1,f[n-1]);

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
