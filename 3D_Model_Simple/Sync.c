// -*- c -*-
#include "public.h"

void Sync
(
 int nol, int non, int noe, 
 struct vert_geom *veG_h,
 struct cell *ceR_h
 )
{
  //printf("Sync1\n");
  static int l,n,i;

  for(n=0;n<noe;n++){
    ceR_h[n].h=ceR_h[n].hn;
  }

  //printf("Sync2\n");
  for (n = 0; n<non; n++){
	  //printf("%d\n", n);
	  veG_h[n].h = ((DATA_TYPE)0.0);
    for(i=0;i<veG_h[n].ceN;i++){
		//printf("%d %d\n", i, veG_h[n].ceN);
		veG_h[n].h += ceR_h[veG_h[n].ce[i]].h;
    }
    veG_h[n].h/=veG_h[n].ceN;
  }
  //printf("Sync3\n");
  for (n = 0; n<non; n++){
    veG_h[n].z[nol+1]=veG_h[n].z[nol]+veG_h[n].h;
  }

  
  //printf("Sync4\n");
  for (l = 0; l<nol; l++){
    for(n=0;n<noe;n++){
      ceR_h[n].u[l]=ceR_h[n].un[l];
      ceR_h[n].v[l]=ceR_h[n].vn[l];
      ceR_h[n].w[l]=ceR_h[n].wn[l];
      ceR_h[n].s[l]=ceR_h[n].sn[l];
      ceR_h[n].r[l]=ceR_h[n].rn[l];
      ceR_h[n].p[l]=ceR_h[n].pn[l];
      ceR_h[n].k[l]=ceR_h[n].kn[l];
      ceR_h[n].e[l]=ceR_h[n].en[l];
    }
  }
  
  return;    
}
