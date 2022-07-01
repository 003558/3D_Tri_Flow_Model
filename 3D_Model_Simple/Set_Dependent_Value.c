// -*- c -*-
#include "public.h"

void Set_Dependent_Value
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound Bc
 )
{
  static int n,l,e;
  static DATA_TYPE dq[3],f1,f2;


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

      ceR_h[n].dux[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      ceR_h[n].duy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];


      for(e=0;e<3;e++){
	if(suG_h[ceC_h[n].su[e]].bd[l]==0){
	  dq[e]=-((DATA_TYPE)2.0)*ceR_h[n].v[l];
	}
	else{
	  dq[e]=ceR_h[ceC_h[n].ce[e]].v[l]-ceR_h[n].v[l];
	}
      }

      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].dvx[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      ceR_h[n].dvy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];


      for(e=0;e<3;e++){
	if(suG_h[ceC_h[n].su[e]].bd[l]==0){
	  dq[e]=-((DATA_TYPE)2.0)*ceR_h[n].w[l];
	}
	else{
	  dq[e]=ceR_h[ceC_h[n].ce[e]].w[l]-ceR_h[n].w[l];
	}
      }

      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].dwx[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      ceR_h[n].dwy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
    }
  }

 
  //duz,dvz,dwz,drz
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol-1;l++){
      ceR_h[n].duz[l]=((DATA_TYPE)0.5)*((ceR_h[n].u[l+1]-ceR_h[n].u[l])/ceG_h[n].dz[l]+(ceR_h[n].u[l]-ceR_h[n].u[l-1])/ceG_h[n].dz[l]);
      ceR_h[n].dvz[l]=((DATA_TYPE)0.5)*((ceR_h[n].v[l+1]-ceR_h[n].v[l])/ceG_h[n].dz[l]+(ceR_h[n].v[l]-ceR_h[n].v[l-1])/ceG_h[n].dz[l]);
      ceR_h[n].dwz[l]=((DATA_TYPE)0.5)*((ceR_h[n].w[l+1]-ceR_h[n].w[l])/ceG_h[n].dz[l]+(ceR_h[n].w[l]-ceR_h[n].w[l-1])/ceG_h[n].dz[l]);
      ceR_h[n].drz[l]=((DATA_TYPE)0.5)*((ceR_h[n].r[l+1]-ceR_h[n].r[l])/ceG_h[n].dz[l]+(ceR_h[n].r[l]-ceR_h[n].r[l-1])/ceG_h[n].dz[l]);
    }

    ceR_h[n].duz[ceG_h[n].ns]=(ceR_h[n].u[ceG_h[n].ns+1]-ceR_h[n].u[ceG_h[n].ns])/ceG_h[n].dz[ceG_h[n].ns];
    ceR_h[n].duz[nol-1]=(ceR_h[n].u[nol-1]-ceR_h[n].u[nol-2])/ceG_h[n].dz[nol-1];
    ceR_h[n].dvz[ceG_h[n].ns]=(ceR_h[n].v[ceG_h[n].ns+1]-ceR_h[n].v[ceG_h[n].ns])/ceG_h[n].dz[ceG_h[n].ns];
    ceR_h[n].dvz[nol-1]=(ceR_h[n].v[nol-1]-ceR_h[n].v[nol-2])/ceG_h[n].dz[nol-1];
    ceR_h[n].dwz[ceG_h[n].ns]=(ceR_h[n].w[ceG_h[n].ns+1]-ceR_h[n].w[ceG_h[n].ns])/ceG_h[n].dz[ceG_h[n].ns];
    ceR_h[n].dwz[nol-1]=(ceR_h[n].w[nol-1]-ceR_h[n].w[nol-2])/ceG_h[n].dz[nol-1];
    ceR_h[n].drz[ceG_h[n].ns]=(ceR_h[n].r[ceG_h[n].ns+1]-ceR_h[n].r[ceG_h[n].ns])/ceG_h[n].dz[ceG_h[n].ns];
    ceR_h[n].drz[nol-1]=(ceR_h[n].r[nol-1]-ceR_h[n].r[nol-2])/ceG_h[n].dz[nol-1];
  }


  //Ev
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].Ev[l]=N_mol+C_mu*ceR_h[n].k[l]*ceR_h[n].k[l]/ceR_h[n].e[l];
      ceR_h[n].Ev[l]=MIN(((DATA_TYPE)0.0005),ceR_h[n].Ev[l]);
      
      /* for(e=0;e<3;e++){ */
      /* 	if(suG_h[ceC_h[n].su[e]].bd[l]==0){ */
      /* 	  ceR_h[n].Ev[l]=N_mol; */
      /* 	} */
      /* } */
    }
  }


  //Pr
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      //ceR_h[n].Pr[l]=(ceR_h[n].Ev[l]-N_mol)*(((DATA_TYPE)2.0)*ceR_h[n].dwz[l]*ceR_h[n].dwz[l]+ceR_h[n].duz[l]*ceR_h[n].duz[l]+ceR_h[n].dvz[l]*ceR_h[n].dvz[l]);
      ceR_h[n].Pr[l]=(ceR_h[n].Ev[l]-N_mol)*(((DATA_TYPE)2.0)*(ceR_h[n].dux[l]*ceR_h[n].dux[l]+ceR_h[n].dvy[l]*ceR_h[n].dvy[l]+ceR_h[n].dwz[l]*ceR_h[n].dwz[l])+(ceR_h[n].dvx[l]+ceR_h[n].duy[l])*(ceR_h[n].dvx[l]*ceR_h[n].duy[l])+(ceR_h[n].dwy[l]+ceR_h[n].dvz[l])*(ceR_h[n].dwy[l]+ceR_h[n].dvz[l])+(ceR_h[n].duz[l]+ceR_h[n].dwx[l])*(ceR_h[n].duz[l]+ceR_h[n].dwx[l]));
    }
  }

  //Gk
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].Gk[l]=MIN(((DATA_TYPE)0.0),ceR_h[n].Ev[l]/S_r*GV/ceR_h[n].r[l]*ceR_h[n].drz[l]);
    }
  }


  //Adt
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].u_Adt[l]=((DATA_TYPE)0.0);
      ceR_h[n].v_Adt[l]=((DATA_TYPE)0.0);
      ceR_h[n].k_Adt[l]=ceR_h[n].Pr[l]-ceR_h[n].e[l]+ceR_h[n].Gk[l];
      if(ceR_h[n].k[l]<=((DATA_TYPE)7.6e-26)){
	ceR_h[n].e_Adt[l]=((DATA_TYPE)0.0);
      }
      else{
	ceR_h[n].e_Adt[l]=(C_1*ceR_h[n].Pr[l]-C_2*ceR_h[n].e[l])*ceR_h[n].e[l]/ceR_h[n].k[l]+C_1*(((DATA_TYPE)1.0)-C_3)*ceR_h[n].e[l]/ceR_h[n].k[l]*ceR_h[n].Gk[l];
      }

    }
    ceR_h[n].u_Adt[nol-1]=C_d*R_a/R_w*Bc.Ux*(Bc.Ux*Bc.Ux+Bc.Uy*Bc.Uy);
    ceR_h[n].v_Adt[nol-1]=C_d*R_a/R_w*Bc.Uy*(Bc.Ux*Bc.Ux+Bc.Uy*Bc.Uy);
  }

  return;
}
