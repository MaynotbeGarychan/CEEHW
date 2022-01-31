#include <stdio.h>
#include <stdlib.h>
#define NL 16

void funcSIMD(int n,int ne,const int* cny,const float* mat,const float* u,float* r){
  float utmp[NL],mattmp[NL],rtmp[NL],ttmp[NL];
  int i,ie,cnytmp,ieo;

asm volatile("# clear left hand side vector");
  for( i=0; i < n; ++i )
  {
    r[i]=0.0F;
  }
  
  if (ne%NL!=0)
  {
    printf("ne cannot be divided by %d\n",ne);
    exit(0);
  }
asm volatile("# start ne loop");
  for ( ieo = 0; ieo < ne; ieo += NL)
  {
    for ( ie = 0; ie < NL; ie++)
    {
      cnytmp=cny[ieo+ie];
      utmp[ie]=u[cnytmp];
    }

    for ( ie = 0; ie < NL; ie++)
    {
      ttmp[ie]=1.0F/(utmp[ie]+1.0F);
      rtmp[ie]=ttmp[ie]*ttmp[ie]*ttmp[ie];
      ttmp[ie]=2.0F/(utmp[ie]+2.0F);
      rtmp[ie]+=ttmp[ie]*ttmp[ie]*ttmp[ie];
      ttmp[ie]=3.0F/(utmp[ie]+3.0F);
      rtmp[ie]+=ttmp[ie]*ttmp[ie]*ttmp[ie];
      ttmp[ie]=4.0F/(utmp[ie]+4.0F);
      rtmp[ie]+=ttmp[ie]*ttmp[ie]*ttmp[ie];

      mattmp[ie] = mat[ieo+ie];
      rtmp[ie]/=mattmp[ie];
    }
    for ( ie = 0; ie < NL; ie++)
    {
      cnytmp=cny[ieo+ie];
      r[cnytmp]+=rtmp[ie];
    }
asm volatile("# end ne loop");
  }
}
