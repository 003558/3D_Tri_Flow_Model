// -*- c -*-
#include "public.h"
void w_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 )
{
  static int l,n,e;
  static DATA_TYPE uM,vM,unM,wM;
  static int ce[2];
  static DATA_TYPE qc[2],uc[2],vc[2],un[2],wc[2];
  static DATA_TYPE fq[2];


  // Flux @ Cell Surface
  
  ///////////////////
  //Horizontal Flux
  ///////////////////
  for(n=0;n<nos;n++){
    for(l=suG_h[n].ns;l<nol;l++){
      
      //boundary(Wall)
      if(suG_h[n].bd[l]==0){
	suR_h[n].fl[l]=((DATA_TYPE)0.0);
      }
      else{
	
	//Value_LR
	for(e=0;e<2;e++){
	  ce[e]=suC_h[n].ce[e];
	  qc[e]=ceR_h[ce[e]].w[l];
	  uc[e]=ceR_h[ce[e]].u[l];
	  vc[e]=ceR_h[ce[e]].v[l];
	  un[e]=uc[e]*suG_h[n].nnx+vc[e]*suG_h[n].nny;
	}
	
	//Flux_LR
	for(e=0;e<2;e++){
	  fq[e]=qc[e]*un[e];
	}
	
	//Velocity
	uM=((DATA_TYPE)0.5)*(uc[0]+uc[1]);
	vM=((DATA_TYPE)0.5)*(vc[0]+vc[1]);
	unM=uM*suG_h[n].nnx+vM*suG_h[n].nny;
	
	suR_h[n].fl[l]=((DATA_TYPE)0.5)*suG_h[n].nnl*(fq[0]+fq[1]-FABS(unM)*(qc[1]-qc[0]));
      }
    }
  }
    
  //Advection term
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].w_Adv[l]=((DATA_TYPE)0.0);
      for(e=0;e<3;e++){
	ceR_h[n].w_Adv[l]+=ceG_h[n].sgn[e]*suR_h[ceC_h[n].su[e]].fl[l]*ceG_h[n].dz[l];
      }
      ceR_h[n].w_Adv[l]/=(ceG_h[n].ds*ceG_h[n].dz[l]);
    }
  }

  ///////////////////
  //Vertical Flux
  ///////////////////
  
  for(n=0;n<noe;n++){

    for(l=ceG_h[n].ns;l<=nol;l++){
     
      //Value_LR
      if(l==ceG_h[n].ns){
	qc[0]=ceR_h[n].w[l];
	wc[0]=-ceR_h[n].w[l];
      }
      else{
	qc[0]=ceR_h[n].w[l-1];
	wc[0]=ceR_h[n].w[l-1];
      }
      if(l==nol){
	qc[1]=ceR_h[n].w[l-1];
	wc[1]=ceR_h[n].w[l-1];
      }
      else{
	qc[1]=ceR_h[n].w[l];
	wc[1]=ceR_h[n].w[l];
      }
      
      //Flux_LR
      for(e=0;e<2;e++){
	fq[e]=qc[e]*wc[e];
      }
	
      //Velocity
      //wM=((DATA_TYPE)0.5)*(wc[0]+wc[1]);
      wM=MAX(FABS(wc[0]),FABS(wc[1]));
	
      ceR_h[n].fl[l]=((DATA_TYPE)0.5)*ceG_h[n].ds*(fq[0]+fq[1]-FABS(wM)*(qc[1]-qc[0]));
    }

    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].w_Adv[l]+=-(ceR_h[n].fl[l+1]-ceR_h[n].fl[l])/(ceG_h[n].dz[l]*ceG_h[n].ds);

    }
    
  }
    
  return;    
}
