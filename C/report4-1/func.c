void func(int n,int ne,const int* cny,const float* mat,const float* u,float* r){
  float utmp,mattmp,rtmp,ttmp;
  int i,ie,cnytmp;

asm volatile("# clear left hand side vector");
  for( i=0; i < n; ++i )
    r[i]=0.0F;

asm volatile("# start ne loop");
  for( ie = 0; ie < ne; ++ie){
    cnytmp=cny[ie];
    utmp=u[cnytmp];
    mattmp=mat[ie];
    ttmp=1.0F/(utmp+1.0F);
    rtmp=ttmp*ttmp*ttmp;
    ttmp=2.0F/(utmp+2.0F);
    rtmp+=ttmp*ttmp*ttmp;
    ttmp=3.0F/(utmp+3.0F);
    rtmp+=ttmp*ttmp*ttmp;
    ttmp=4.0F/(utmp+4.0F);
    rtmp+=ttmp*ttmp*ttmp;
    rtmp/=mattmp;
    r[cnytmp]+=rtmp;
  }
asm volatile("# end ne loop");
}
