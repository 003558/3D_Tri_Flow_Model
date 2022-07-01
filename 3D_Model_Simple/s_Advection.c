// -*- c -*-
#include "public.h"
void s_Advection
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
  static int lf[2];

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
	  qc[e]=ceR_h[ce[e]].s[l];
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

	/* if(suR_h[n].fl[l]*DT>=((DATA_TYPE)0.0)){ */
	/*   if(suR_h[n].fl[l]*DT>qc[0]*ceG_h[ce[0]].dz[l]*ceG_h[ce[0]].ds){ */
	/*     suR_h[n].fl[l]=qc[0]*ceG_h[ce[0]].dz[l]*ceG_h[ce[0]].ds/DT; */
	/*   } */
	/* } */
	/* if(suR_h[n].fl[l]*DT<=((DATA_TYPE)0.0)){ */
	/*   if(suR_h[n].fl[l]*DT>qc[1]*ceG_h[ce[1]].dz[l]*ceG_h[ce[1]].ds){ */
	/*     suR_h[n].fl[l]=qc[1]*ceG_h[ce[1]].dz[l]*ceG_h[ce[1]].ds/DT; */
	/*   } */
	/* } */

      }
    }
  }
    
  //Advection term
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].s_Adv[l]=((DATA_TYPE)0.0);
      for(e=0;e<3;e++){
	ceR_h[n].s_Adv[l]+=ceG_h[n].sgn[e]*suR_h[ceC_h[n].su[e]].fl[l]*ceG_h[n].dz[l];
      }
      ceR_h[n].s_Adv[l]/=(ceG_h[n].ds*ceG_h[n].dz[l]);
    }
  }

  ///////////////////
  //Vertical Flux
  ///////////////////
  
  for(n=0;n<noe;n++){
   for(l=ceG_h[n].ns;l<=nol;l++){

      lf[0]=MAX(ceG_h[n].ns,l-1);
      lf[1]=MIN(nol-1,l);

      for(e=0;e<2;e++){
	qc[e]=ceR_h[n].s[lf[e]];
	wc[e]=ceR_h[n].w[lf[e]];
      }
      
      //Flux_LR
      for(e=0;e<2;e++){
  	fq[e]=qc[e]*wc[e];
      }
	
      //Velocity
      wM=MAX(FABS(wc[0]),FABS(wc[1]));
      
      ceR_h[n].fl[l]=((DATA_TYPE)0.5)*ceG_h[n].ds*(fq[0]+fq[1]-FABS(wM)*(qc[1]-qc[0]));

      /* if(ceR_h[n].fl[l]>=((DATA_TYPE)0.0)){ */
      /* 	if(ceR_h[n].fl[l]*DT>=qc[0]*ceG_h[n].dz[lf[0]]*ceG_h[n].ds){ */
      /* 	  ceR_h[n].fl[l]=qc[0]*ceG_h[n].dz[lf[0]]*ceG_h[n].ds/DT; */
      /* 	} */
      /* } */
      /* if(ceR_h[n].fl[l]<=((DATA_TYPE)0.0)){ */
      /* 	if(ceR_h[n].fl[l]*DT>=qc[1]*ceG_h[n].dz[lf[1]]*ceG_h[n].ds){ */
      /* 	  ceR_h[n].fl[l]=qc[1]*ceG_h[n].dz[lf[1]]*ceG_h[n].ds/DT; */
      /* 	} */
      /* } */

   }
   ceR_h[n].fl[ceG_h[n].ns]=((DATA_TYPE)0.0);
   ceR_h[n].fl[nol]=((DATA_TYPE)0.0);
   
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].s_Adv[l]+=-(ceR_h[n].fl[l+1]-ceR_h[n].fl[l])/(ceG_h[n].dz[l]*ceG_h[n].ds);
    }
    
  }
    
  return;    
}
