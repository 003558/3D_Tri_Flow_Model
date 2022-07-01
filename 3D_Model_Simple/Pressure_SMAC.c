// -*- c -*-
#include "public.h"

void Pressure
(
 int nol, int nos,int noe, 
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound Bc
 )
{
  static int n,l,e;
  static int it;
  static int su[3],ce[3],ve[3];

  static DATA_TYPE x[3],y[3],u[3],v[3];

  static DATA_TYPE r0,r1;
  static DATA_TYPE c1,c2;

  static DATA_TYPE *a,*b,*c,*d;

  static DATA_TYPE KE;
  static DATA_TYPE Sp_Scale_Turb;
  static DATA_TYPE epse;

  a = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  b = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  c = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  d = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);

  //Previous Value
  for(n=0;n<noe;n++){
    ceR_h[n].hn=ceR_h[n].h+DT*(ceR_h[n].h_Adv+ceR_h[n].w[nol-1]);
    ceR_h[n].hm=ceR_h[n].hn;
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].un[l]=ceR_h[n].u[l]+DT*(ceR_h[n].u_Adv[l]+ceR_h[n].u_Adt[l]+ceR_h[n].u_Dif[l]);
      ceR_h[n].vn[l]=ceR_h[n].v[l]+DT*(ceR_h[n].v_Adv[l]+ceR_h[n].v_Adt[l]+ceR_h[n].v_Dif[l]);
      ceR_h[n].wn[l]=ceR_h[n].w[l]+DT*(ceR_h[n].w_Adv[l]-GV+ceR_h[n].w_Dif[l]);
      ceR_h[n].sn[l]=ceR_h[n].s[l]+DT*(ceR_h[n].s_Adv[l]+ceR_h[n].s_Dif[l]);
      //ceR_h[n].sn[l]=ceR_h[n].s[l]+DT*(ceR_h[n].s_Adv[l]);
      ceR_h[n].kn[l]=ceR_h[n].k[l]+DT*(ceR_h[n].k_Adv[l]+ceR_h[n].k_Dif[l]+ceR_h[n].k_Adt[l]);
      ceR_h[n].en[l]=ceR_h[n].e[l]+DT*(ceR_h[n].e_Adv[l]+ceR_h[n].e_Dif[l]+ceR_h[n].e_Adt[l]);
    }
    ceR_h[n].kn[nol-1]=(R_a*C_d*(Bc.Ux*Bc.Ux+Bc.Uy*Bc.Uy))/(R_w*sqrt(C_mu));
    ceR_h[n].en[nol-1]=POW(sqrt(C_mu)*ceR_h[n].kn[nol-1],((DATA_TYPE)1.5))/(((DATA_TYPE)0.5)*Kappa*(ceR_h[n].hn+((DATA_TYPE)0.5)*ceG_h[n].dz[nol-1]));
  }


  for(n=0;n<noe;n++){
    Sp_Scale_Turb=ceR_h[n].h+((DATA_TYPE)0.5)*ceG_h[n].dz[nol-1]+ceG_h[n].zc[nol-1]-ceG_h[n].zc[ceG_h[n].ns];
    for(l=ceG_h[n].ns;l<nol;l++){
      KE=((DATA_TYPE)7.6e-26);
      if(ceR_h[n].kn[l]<=KE){
  	ceR_h[n].kn[l]=KE;
  	ceR_h[n].en[l]=POW(ceR_h[n].kn[l]*sqrt(C_mu),((DATA_TYPE)1.5))/Sp_Scale_Turb;
      }
      epse=POW(sqrt(C_mu)*ceR_h[n].kn[l],((DATA_TYPE)1.5))/Sp_Scale_Turb;
      if(ceR_h[n].en[l]<epse){
  	ceR_h[n].en[l]=epse;
      }
    }
  }

  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns+1;l<nol;l++){
      ceR_h[n].sn[l]=MAX(Smin,MIN(Smax,ceR_h[n].sn[l]));
      ceR_h[n].rn[l]=Density_Value(ceR_h[n].sn[l],((DATA_TYPE)20.0),ceR_h[n].p[l]);
    }
    l=ceG_h[n].ns;
    ceR_h[n].sn[l]=MAX(Smin,MIN(Smax,ceR_h[n].sn[l]));
    ceR_h[n].rn[l]=ceR_h[n].rn[l+1];
  }

  //Boundary
  for(n=0;n<noe;n++){
    for(l=0;l<ceG_h[n].ns;l++){
      ceR_h[n].sn[l]=ceR_h[n].sn[ceG_h[n].ns];
      ceR_h[n].en[l]=ceR_h[n].en[ceG_h[n].ns];
      ceR_h[n].kn[l]=ceR_h[n].kn[ceG_h[n].ns];
      ceR_h[n].rn[l]=ceR_h[n].rn[ceG_h[n].ns];
    }
  }

  //u,v,r at cell boundary
  for(n=0;n<nos;n++){
    ce[0]=suC_h[n].ce[0];
    ce[1]=suC_h[n].ce[1];
    for(l=suG_h[n].ns;l<nol;l++){
      if(suG_h[n].bd[l]==0){
  	suR_h[n].ub[l]=((DATA_TYPE)0.0);
  	suR_h[n].vb[l]=((DATA_TYPE)0.0);
  	if(ceG_h[ce[0]].bb[l]==0){
  	  suR_h[n].rb[l]=ceR_h[ce[0]].rn[l];
  	}
  	else{
  	  suR_h[n].rb[l]=ceR_h[ce[1]].rn[l];
  	}
      }
      else{
  	r0=ceG_h[ce[0]].ds/(ceG_h[ce[0]].ds+ceG_h[ce[1]].ds);;
  	r1=ceG_h[ce[1]].ds/(ceG_h[ce[0]].ds+ceG_h[ce[1]].ds);;
  	suR_h[n].ub[l]=r0*ceR_h[ce[0]].un[l]+r1*ceR_h[ce[1]].un[l];
  	suR_h[n].vb[l]=r0*ceR_h[ce[0]].vn[l]+r1*ceR_h[ce[1]].vn[l];
  	suR_h[n].rb[l]=r0*ceR_h[ce[0]].rn[l]+r1*ceR_h[ce[1]].rn[l];

	if(suR_h[n].rb[l]>((DATA_TYPE)0.0)){
	  suR_h[n].ub[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].p[l]-ceR_h[ce[0]].p[l])*suG_h[n].nnx*suG_h[n].nnl;
	  suR_h[n].vb[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].p[l]-ceR_h[ce[0]].p[l])*suG_h[n].nny*suG_h[n].nnl;
	}
      }

    }

  }

  //w at cell boundary
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns+1;l<nol;l++){
      ceR_h[n].wb[l]=((DATA_TYPE)0.5)*(ceR_h[n].wn[l-1]+ceR_h[n].wn[l]);
      ceR_h[n].rb[l]=((DATA_TYPE)0.5)*(ceR_h[n].rn[l-1]+ceR_h[n].rn[l]);

      if(ceR_h[n].rb[l]>((DATA_TYPE)0.0)){
	ceR_h[n].wb[l]-=DT/(ceR_h[n].rb[l]*ceG_h[n].dz[l])*(ceR_h[n].p[l]-ceR_h[n].p[l-1]);
      }
    }
    
    //ceR_h[n].wb[ceG_h[n].ns]=((DATA_TYPE)0.0);
    ceR_h[n].wb[ceG_h[n].ns]=ceG_h[n].bzdx*ceR_h[n].un[ceG_h[n].ns]+ceG_h[n].bzdy*ceR_h[n].vn[ceG_h[n].ns];
    ceR_h[n].wb[nol]=(ceR_h[n].hn-ceR_h[n].h)/DT;
    ceR_h[n].rb[ceG_h[n].ns]=ceR_h[n].rn[ceG_h[n].ns];
    ceR_h[n].rb[nol]=ceR_h[n].rn[nol-1];
  }

  //It
  for(it=0;it<10;it++){

    for(n=0;n<noe;n++){
      for(l=ceG_h[n].ns;l<nol;l++){
  	for(e=0;e<3;e++){
  	  su[e]=ceC_h[n].su[e];
  	  ve[e]=ceC_h[n].ve[e];

  	  x[e]=veG_h[ve[e]].x;
  	  y[e]=veG_h[ve[e]].y;

  	  u[e]=suR_h[su[e]].ub[l];
  	  v[e]=suR_h[su[e]].vb[l];
  	}
	
  	ceR_h[n].D[l]=(u[0]*(y[1]-y[0])+u[1]*(y[2]-y[1])+u[2]*(y[0]-y[2])+v[0]*(x[0]-x[1])+v[1]*(x[1]-x[2])+v[2]*(x[2]-x[0]))/ceG_h[n].ds;

  	ceR_h[n].D[l]+=(ceR_h[n].wb[l+1]-ceR_h[n].wb[l])/ceG_h[n].dz[l];
      }
    }


    // Pressure Implicit method
    for(n=0;n<noe;n++){
      
      //Initial
      l=ceG_h[n].ns;
      c2=((DATA_TYPE)2.0)/(ceR_h[n].rb[l+1]*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l])*ceG_h[n].dz[l]);
      c[l]=((DATA_TYPE)0.0);
      a[l]=-c2;
      b[l]=c2;
      d[l]=ceR_h[n].D[l]/DT;
      
      for(l=ceG_h[n].ns+1;l<nol-1;l++){
      	c1=((DATA_TYPE)2.0)/(ceR_h[n].rb[l]*(ceG_h[n].dz[l]+ceG_h[n].dz[l-1])*ceG_h[n].dz[l]);
      	c2=((DATA_TYPE)2.0)/(ceR_h[n].rb[l+1]*(ceG_h[n].dz[l+1]+ceG_h[n].dz[l])*ceG_h[n].dz[l]);

      	c[l]=c1;
      	a[l]=-(c1+c2);
      	b[l]=c2;
      	d[l]=ceR_h[n].D[l]/DT;
      }
      
      l=nol-1;
      c[l]=((DATA_TYPE)0.0);
      a[l]=((DATA_TYPE)1.0);
      b[l]=((DATA_TYPE)0.0);
      d[nol-1]=ceR_h[n].rn[nol-1]*GV*(ceR_h[n].hn-ceR_h[n].h);

      //LU
      b[ceG_h[n].ns]/=a[ceG_h[n].ns];
      d[ceG_h[n].ns]/=a[ceG_h[n].ns];

      for(l=ceG_h[n].ns+1;l<nol;l++){
      	a[l]-=c[l]*b[l-1];
      	d[l]-=c[l]*d[l-1];
      	d[l]/=a[l];
      	b[l]/=a[l];
      }
      for(l=nol-2;l>=ceG_h[n].ns;l--){
      	d[l]-=b[l]*d[l+1];
      }


      //Update P
      for(l=ceG_h[n].ns;l<nol;l++){
  	ceR_h[n].Psi[l]=d[l];
      }
      
      for(l=ceG_h[n].ns;l<nol;l++){
  	ceR_h[n].pn[l]=ceR_h[n].p[l]+ceR_h[n].Psi[l];
      }
    }

    //Update ub,vb
    for(n=0;n<nos;n++){
      ce[0]=suC_h[n].ce[0];
      ce[1]=suC_h[n].ce[1];

      for(l=ceG_h[ce[0]].ns;l<nol;l++){
  	if(suG_h[n].bd[l]!=0){
  	  suR_h[n].ub[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].Psi[l]-ceR_h[ce[0]].Psi[l])*suG_h[n].nnx*suG_h[n].nnl;
  	  suR_h[n].vb[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].Psi[l]-ceR_h[ce[0]].Psi[l])*suG_h[n].nny*suG_h[n].nnl;
  	}
  	else{
  	  suR_h[n].ub[l]=((DATA_TYPE)0.0);
  	  suR_h[n].vb[l]=((DATA_TYPE)0.0);
  	}

      }
    }

    //Update wb
    for(n=0;n<noe;n++){
      for(l=ceG_h[n].ns+1;l<nol;l++){
  	ceR_h[n].wb[l]-=DT/(ceR_h[n].rb[l]*ceG_h[n].dz[l])*(ceR_h[n].Psi[l]-ceR_h[n].Psi[l-1]);
      }
      l=nol;
    }


    for(n=0;n<noe;n++){
      for(l=ceG_h[n].ns;l<nol-1;l++){
  	ceR_h[n].p[l]=ceR_h[n].pn[l];
      }
    }
  }


  //Update u,v
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns;l<nol;l++){
      ceR_h[n].un[l]=(suR_h[ceC_h[n].su[0]].ub[l]+suR_h[ceC_h[n].su[1]].ub[l]+suR_h[ceC_h[n].su[2]].ub[l])/((DATA_TYPE)3.0);
      ceR_h[n].vn[l]=(suR_h[ceC_h[n].su[0]].vb[l]+suR_h[ceC_h[n].su[1]].vb[l]+suR_h[ceC_h[n].su[2]].vb[l])/((DATA_TYPE)3.0);
    }
  }

  //Update w
  for(n=0;n<noe;n++){
    for(l=ceG_h[n].ns+1;l<nol-1;l++){
      ceR_h[n].wn[l]=((DATA_TYPE)0.5)*(ceR_h[n].wb[l]+ceR_h[n].wb[l+1]);
    }
    ceR_h[n].wn[ceG_h[n].ns]=((DATA_TYPE)0.0);
    ceR_h[n].wn[nol-1]=ceR_h[n].wb[nol-1];
  }
  
  return;    
}
