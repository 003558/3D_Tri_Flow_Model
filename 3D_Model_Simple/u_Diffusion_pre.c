// -*- c -*-
#include "public.h"
void u_Diffusion
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
  static DATA_TYPE dq[3],dqx[3],dqy[3];

  ///////////////////////
  //Horizontal Diffusion
  ///////////////////////

  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      
      for(e=0;e<3;e++){
	if(suG_h[ceC_h[n].su[e]].bd[l]==0){
	  dq[e]=-((DATA_TYPE)2.0)*ceR_h[n].u[l];
	}
	else{
	  dq[e]=ceR_h[ceC_h[n].ce[e]].u[l]-ceR_h[n].u[l];
	}
      }
      
      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].qx[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      ceR_h[n].qy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
      
    }
    
    for(n=0;n<noe;n++){

      for(e=0;e<3;e++){
	if(suG_h[ceC_h[n].su[e]].bd[l]==0){
	  dqx[e]=((DATA_TYPE)0.0);
	  dqy[e]=((DATA_TYPE)0.0);
	}
	else{
	  dqx[e]=((DATA_TYPE)0.5)*(ceR_h[ceC_h[n].ce[e]].Eh[l]+ceR_h[n].Eh[l])*(ceR_h[ceC_h[n].ce[e]].qx[l]-ceR_h[n].qx[l]);
	  dqy[e]=((DATA_TYPE)0.5)*(ceR_h[ceC_h[n].ce[e]].Eh[l]+ceR_h[n].Eh[l])*(ceR_h[ceC_h[n].ce[e]].qy[l]-ceR_h[n].qy[l]);
	}
      }

      f1=ceG_h[n].dx[0]*dqx[0]+ceG_h[n].dx[1]*dqx[1]+ceG_h[n].dx[2]*dqx[2];
      f2=ceG_h[n].dy[0]*dqx[0]+ceG_h[n].dy[1]*dqx[1]+ceG_h[n].dy[2]*dqx[2];

      qxx=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];

      f1=ceG_h[n].dx[0]*dqy[0]+ceG_h[n].dx[1]*dqy[1]+ceG_h[n].dx[2]*dqy[2];
      f2=ceG_h[n].dy[0]*dqy[0]+ceG_h[n].dy[1]*dqy[1]+ceG_h[n].dy[2]*dqy[2];
      
      qyy=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
      
      ceR_h[n].u_Dif[l]=(qxx+qyy);
    }
  }


  /////////////////////
  //Vertical Diffusion
  /////////////////////
 
  for(n=0;n<noe;n++){
    if(ceG_h[n].ns+1>=nol-1){
      ceR_h[n].u_Dif[nol-1]=((DATA_TYPE)0.0);
    }
    else{
      for(l=ceG_h[n].ns+1;l<nol-1;l++){
	Ev1=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l+1]+ceR_h[n].Ev[l]);
	Ev2=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l]+ceR_h[n].Ev[l-1]);
	dq1=ceR_h[n].u[l+1]-ceR_h[n].u[l];
	dq2=ceR_h[n].u[l]-ceR_h[n].u[l-1];
	dz1=((DATA_TYPE)0.5)*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l]);
	dz2=((DATA_TYPE)0.5)*(ceG_h[n].dz[l]+ceG_h[n].dz[l-1]);
	
	ceR_h[n].u_Dif[l]+=(Ev1*dq1/dz1-Ev2*dq2/dz2)/ceG_h[n].dz[l];
      }
      
      //bottom
      l=ceG_h[n].ns;
      Ev1=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l+1]+ceR_h[n].Ev[l]);
      dq1=ceR_h[n].u[l+1]-ceR_h[n].u[l];
      dz1=((DATA_TYPE)0.5)*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l]);
      
      ceR_h[n].u_Dif[l]+=(Ev1*dq1/dz1-((DATA_TYPE)0.0))/ceG_h[n].dz[l];
      
      //surface
      l=nol-1;
      Ev2=((DATA_TYPE)0.5)*(ceR_h[n].Ev[l]+ceR_h[n].Ev[l-1]);
      dq2=ceR_h[n].u[l]-ceR_h[n].u[l-1];
      dz2=((DATA_TYPE)0.5)*(ceG_h[n].dz[l]+ceG_h[n].dz[l-1]);
      
      ceR_h[n].u_Dif[l]+=(((DATA_TYPE)0.0)-Ev2*dq2/dz2)/ceG_h[n].dz[l];
    }
  }
  
  return;    
}
