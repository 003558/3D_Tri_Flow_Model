// -*- c -*-
#include "public.h"
////////////////
//Check
////////////////
void Cal_Info
(
 int step, int step_o, DATA_TYPE step_t,
 int nol, int nos,int noe, 
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 )
{
  static int l,n;
  static DATA_TYPE q_total;
  static DATA_TYPE dl,Vel;
  static DATA_TYPE CFL_H,CFL_H_max;
  static DATA_TYPE CFL_V,CFL_V_max;

  //Mass
  q_total=((DATA_TYPE)0.0);
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      q_total+=ceR_h[n].r[l]*ceG_h[n].ds*ceG_h[n].dz[l];
    }
  }
  
  //CFL number of Horizontal
  CFL_H=((DATA_TYPE)0.0);
  CFL_H_max=((DATA_TYPE)0.0);

  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      dl=sqrt(((DATA_TYPE)2.0)*ceG_h[n].ds);
      Vel=sqrt(ceR_h[n].u[l]*ceR_h[n].u[l]+ceR_h[n].v[l]*ceR_h[n].v[l]);
      
      CFL_H=Vel*DT/dl;
      
      if(CFL_H>CFL_H_max){
	CFL_H_max=CFL_H;
      }
    }
  }

  //CFL number of Vertical
  CFL_V=((DATA_TYPE)0.0);
  CFL_V_max=((DATA_TYPE)0.0);

  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      dl=ceG_h[n].dz[l];
      Vel=FABS(ceR_h[n].w[l]);
      
      CFL_V=Vel*DT/dl;
      
      if(CFL_V>CFL_V_max){
	CFL_V_max=CFL_V;
      }
    }
  }

  printf("#######################################\n");
  printf(" Step:%d Time=%6.2e(s)\n",step,step_t);

  printf("  Cal_Check: Mass=%6.4e, CFL_H=%6.4e, CFL_V=%6.4e\n",q_total,CFL_H_max,CFL_V_max);
  printf("  Output File: t_%04d.plt\n",step_o);
  printf("#######################################\n");
  return;
}
