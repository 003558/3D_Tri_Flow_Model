// -*- c -*-
#include "public.h"
void w_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 )
{
  static int n,l;
  static DATA_TYPE Ev1,Ev2,dq1,dq2,dz1,dz2;

  static int e;
  static DATA_TYPE f1,f2;
  static DATA_TYPE qxx,qyy;
  static DATA_TYPE dq,dx,dy;
  static DATA_TYPE dqx[3],dqy[3];

  ///////////////////////
  //Horizontal Diffusion
  ///////////////////////

  for(l=0;l<nol;l++){
    for(n=0;n<nos;n++){
      if(suG_h[n].bd[l]==0){
	dq=-((DATA_TYPE)2.0)*ceR_h[suC_h[n].ce[0]].w[l];
	}
	else{
	dq=ceR_h[suC_h[n].ce[1]].w[l]-ceR_h[suC_h[n].ce[1]].w[l];
      }
      
      dx=ceG_h[suC_h[n].ce[1]].xc-ceG_h[suC_h[n].ce[0]].xc;
      dy=ceG_h[suC_h[n].ce[1]].yc-ceG_h[suC_h[n].ce[0]].yc;
      
      if(dx==((DATA_TYPE)0.0)){
	suR_h[n].dqx=((DATA_TYPE)0.0);
      }
      else{
	suR_h[n].dqx=dq/dx;
      }
      if(dy==((DATA_TYPE)0.0)){
	suR_h[n].dqy=((DATA_TYPE)0.0);
      }
      else{
	suR_h[n].dqy=dq/dy;
      }
    }
    
    for(n=0;n<noe;n++){

      for(e=0;e<3;e++){
	dqx[e]=((DATA_TYPE)0.5)*(ceR_h[ceC_h[n].ce[e]].Eh[l]+ceR_h[n].Eh[l])*suR_h[ceC_h[n].su[e]].dqx;
	dqy[e]=((DATA_TYPE)0.5)*(ceR_h[ceC_h[n].ce[e]].Eh[l]+ceR_h[n].Eh[l])*suR_h[ceC_h[n].su[e]].dqy;
      }

      f1=ceG_h[n].dx[0]*dqx[0]+ceG_h[n].dx[1]*dqx[1]+ceG_h[n].dx[2]*dqx[2];
      f2=ceG_h[n].dy[0]*dqx[0]+ceG_h[n].dy[1]*dqx[1]+ceG_h[n].dy[2]*dqx[2];

      qxx=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];

      f1=ceG_h[n].dx[0]*dqy[0]+ceG_h[n].dx[1]*dqy[1]+ceG_h[n].dx[2]*dqy[2];
      f2=ceG_h[n].dy[0]*dqy[0]+ceG_h[n].dy[1]*dqy[1]+ceG_h[n].dy[2]*dqy[2];
      
      qyy=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
      
      ceR_h[n].w_Dif[l]=(qxx+qyy);
    }
  }

  /////////////////////
  //Vertical Diffusion
  /////////////////////


  for(n=0;n<noe;n++){
    if(ceG_h[n].ns+1>=nol-1){
      ceR_h[n].w_Dif[nol-1]=((DATA_TYPE)0.0);
    }
    else{
      for(l=ceG_h[n].ns+1;l<nol-1;l++){
	Ev1=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l+1]+ceR_h[n].Ev[l]);
	Ev2=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l]+ceR_h[n].Ev[l-1]);
	dq1=ceR_h[n].w[l+1]-ceR_h[n].w[l];
	dq2=ceR_h[n].w[l]-ceR_h[n].w[l-1];
	dz1=((DATA_TYPE)0.5)*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l]);
	dz2=((DATA_TYPE)0.5)*(ceG_h[n].dz[l]+ceG_h[n].dz[l-1]);
	
	ceR_h[n].w_Dif[l]+=(Ev1*dq1/dz1-Ev2*dq2/dz2)/ceG_h[n].dz[l];
      }
      
      //bottom
      l=ceG_h[n].ns;
      Ev1=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l+1]+ceR_h[n].Ev[l]);
      dq1=ceR_h[n].w[l+1]-ceR_h[n].w[l];
      dz1=((DATA_TYPE)0.5)*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l]);
      
      ceR_h[n].w_Dif[l]+=(Ev1*dq1/dz1-((DATA_TYPE)0.0))/ceG_h[n].dz[l];
      
      //surface
      l=nol-1;
      Ev2=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l]+ceR_h[n].Ev[l-1]);
      dq2=ceR_h[n].w[l]-ceR_h[n].w[l-1];
      dz2=((DATA_TYPE)0.5)*(ceG_h[n].dz[l]+ceG_h[n].dz[l-1]);
      
      ceR_h[n].w_Dif[l]+=(((DATA_TYPE)0.0)-Ev2*dq2/dz2)/ceG_h[n].dz[l];
    }
  }
  
  return;    
}
