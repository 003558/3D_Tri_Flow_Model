// -*- c -*-
#include "public.h"


void Pressure
(
 int nol, int nos,int noe, 
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 )
{
  static int n,m,l,e;
  static int it,ik;
  static int ch;
  static int su[3],ce[3];

  static DATA_TYPE dq[3];
  static DATA_TYPE dqx,dqy,dqz;
  static DATA_TYPE f1,f2;
  static DATA_TYPE r0,r1;

  static DATA_TYPE *a,*b,*c,*d;

  a = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  b = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  c = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);
  d = (DATA_TYPE *)malloc(sizeof(DATA_TYPE)*nol);

  //Previous Value
  for(n=0;n<noe;n++){
    ceR_h[n].hn=ceR_h[n].h+DT*(ceR_h[n].h_Adv+ceR_h[n].w[nol-1]);
    for(l=0;l<nol;l++){
      ceR_h[n].un[l]=ceR_h[n].u[l]+DT*(ceR_h[n].u_Adv[l]+ceR_h[n].u_Dif[l]+ceR_h[n].u_Adt[l]);
      ceR_h[n].vn[l]=ceR_h[n].v[l]+DT*(ceR_h[n].v_Adv[l]+ceR_h[n].v_Dif[l]+ceR_h[n].v_Adt[l]);
      ceR_h[n].wn[l]=ceR_h[n].w[l]+DT*(ceR_h[n].w_Adv[l]+ceR_h[n].w_Dif[l]-GV);
      //ceR_h[n].rn[l]=ceR_h[n].r[l]+DT*(ceR_h[n].r_Adv[l]+ceR_h[n].r_Dif[l]);
      ceR_h[n].kn[l]=ceR_h[n].k[l]+DT*(ceR_h[n].k_Adv[l]+ceR_h[n].k_Dif[l]+ceR_h[n].k_Adt[l]);
      ceR_h[n].en[l]=ceR_h[n].e[l]+DT*(ceR_h[n].e_Adv[l]+ceR_h[n].e_Dif[l]+ceR_h[n].e_Adt[l]);

    }
  }

  //u,v,r at cell boundary
  for(n=0;n<nos;n++){
    for(l=0;l<nol;l++){
      ce[0]=suC_h[n].ce[0];
      ce[1]=suC_h[n].ce[1];
      if(suG_h[n].bd==0){
	suR_h[n].ub[l]=((DATA_TYPE)0.0);
	suR_h[n].vb[l]=((DATA_TYPE)0.0);
	suR_h[n].rb[l]=((DATA_TYPE)0.5)*(ceR_h[ce[0]].rn[l]+ceR_h[ce[1]].rn[l]);
      }
      else{
	r0=ceG_h[ce[0]].ds/(ceG_h[ce[0]].ds+ceG_h[ce[1]].ds);;
	r1=ceG_h[ce[1]].ds/(ceG_h[ce[0]].ds+ceG_h[ce[1]].ds);;
	suR_h[n].ub[l]=r0*ceR_h[ce[0]].un[l]+r1*ceR_h[ce[1]].un[l];
	suR_h[n].vb[l]=r0*ceR_h[ce[0]].vn[l]+r1*ceR_h[ce[1]].vn[l];
	suR_h[n].rb[l]=r0*ceR_h[ce[0]].rn[l]+r1*ceR_h[ce[1]].rn[l];

	suR_h[n].ub[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].p[l]-ceR_h[ce[0]].p[l])*suG_h[n].nnx*suG_h[n].nnl;
	suR_h[n].vb[l]-=DT/(suR_h[n].rb[l]*suG_h[n].ds)*(ceR_h[ce[1]].p[l]-ceR_h[ce[0]].p[l])*suG_h[n].nny*suG_h[n].nnl;
      }
    }
    suR_h[n].ub[0]=((DATA_TYPE)0.0);
    suR_h[n].vb[0]=((DATA_TYPE)0.0);
  }

  
  //w at cell boundary
  for(n=0;n<noe;n++){
    for(l=1;l<nol;l++){
      ceR_h[n].wb[l]=((DATA_TYPE)0.5)*(ceR_h[n].wn[l-1]+ceR_h[n].wn[l]);
      ceR_h[n].rb[l]=((DATA_TYPE)0.5)*(ceR_h[n].rn[l-1]+ceR_h[n].rn[l]);
      
      ceR_h[n].wb[l]-=DT/(ceR_h[n].rb[l]*ceG_h[n].dz[l])*(ceR_h[n].p[l]-ceR_h[n].p[l-1]);
    }

    ceR_h[n].wb[0]=((DATA_TYPE)0.0);
    //ceR_h[n].wb[nol]=ceR_h[n].wn[nol-1]+DT*GV;
    ceR_h[n].wb[nol]=(ceR_h[n].hn-ceR_h[n].h)/DT;
  }

  //Divergence V at cell center
  for(n=0;n<noe;n++){
    for(l=0;l<nol;l++){
      //dudx
      f1=ceG_h[n].dx[0]*suR_h[ceC_h[n].su[0]].ub[l]+ceG_h[n].dx[1]*suR_h[ceC_h[n].su[1]].ub[l]+ceG_h[n].dx[2]*suR_h[ceC_h[n].su[2]].ub[l];
      f2=ceG_h[n].dy[0]*suR_h[ceC_h[n].su[0]].ub[l]+ceG_h[n].dy[1]*suR_h[ceC_h[n].su[1]].ub[l]+ceG_h[n].dy[2]*suR_h[ceC_h[n].su[2]].ub[l];
      
      dqx=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      
      //dvdy
      f1=ceG_h[n].dx[0]*suR_h[ceC_h[n].su[0]].vb[l]+ceG_h[n].dx[1]*suR_h[ceC_h[n].su[1]].vb[l]+ceG_h[n].dx[2]*suR_h[ceC_h[n].su[2]].vb[l];
      f2=ceG_h[n].dy[0]*suR_h[ceC_h[n].su[0]].vb[l]+ceG_h[n].dy[1]*suR_h[ceC_h[n].su[1]].vb[l]+ceG_h[n].dy[2]*suR_h[ceC_h[n].su[2]].vb[l];
      
      dqy=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
      
      //dwdz
      dqz=(ceR_h[n].wb[l+1]-ceR_h[n].wb[l])/ceG_h[n].dz[l];
      
      ceR_h[n].D[l]=dqx+dqy+dqz;
    }
  }

  //Horizontal grad of Pressure
  for(l=0;l<nol;l++){
    //dpx,dpy
    for(n=0;n<noe;n++){
      for(e=0;e<3;e++){
  	if(suG_h[ceC_h[n].su[e]].bd==0){
  	  dq[e]=((DATA_TYPE)0.0);
  	}
  	else{
  	  dq[e]=ceR_h[ceC_h[n].ce[e]].p[l]-ceR_h[n].p[l];
  	}
      }
      
      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].dpx[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
      ceR_h[n].dpy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
    }

    //ddpx
    for(n=0;n<noe;n++){
      for(e=0;e<3;e++){
  	if(suG_h[ceC_h[n].su[e]].bd==0){
  	  dq[e]=((DATA_TYPE)0.0);
  	}
  	else{
  	  dq[e]=ceR_h[ceC_h[n].ce[e]].dpx[l]-ceR_h[n].dpx[l];
  	}
      }
      
      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].ddpx[l]=(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3];
    }

    //ddpy
    for(n=0;n<noe;n++){
      for(e=0;e<3;e++){
  	if(suG_h[ceC_h[n].su[e]].bd==0){
  	  dq[e]=((DATA_TYPE)0.0);
  	}
  	else{
  	  dq[e]=ceR_h[ceC_h[n].ce[e]].dpy[l]-ceR_h[n].dpy[l];
  	}
      }
      
      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];

      ceR_h[n].ddpy[l]=(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3];
    }

    //Mod Divergence
    /* for(n=0;n<noe;n++){ */
    /*   ceR_h[n].D[l]-=(ceR_h[n].ddpx[l]+ceR_h[n].ddpy[l])*ceG_h[n].dz[l]*ceG_h[n].dz[l]*DT; */
    /* } */
  }

  
  for(n=0;n<noe;n++){
    
    //Initial
    c[0]=((DATA_TYPE)0.0);
    a[0]=-((DATA_TYPE)1.0);
    b[0]=((DATA_TYPE)1.0);
    d[0]=ceR_h[n].rn[0]*ceR_h[n].D[0]/DT*ceG_h[n].dz[0]*ceG_h[n].dz[0];

    for(l=ceG_h[n].ns+1;l<nol-1;l++){
      c[l]=((DATA_TYPE)1.0);
      a[l]=-((DATA_TYPE)2.0);
      b[l]=((DATA_TYPE)1.0);
      d[l]=ceR_h[n].rn[l]*ceR_h[n].D[l]/DT*ceG_h[n].dz[l]*ceG_h[n].dz[l];
    }

    c[nol-1]=((DATA_TYPE)0.0);
    a[nol-1]=((DATA_TYPE)1.0);
    b[nol-1]=((DATA_TYPE)0.0);
    d[nol-1]=ceR_h[n].rn[nol-1]*GV*(ceR_h[n].hn-ceR_h[n].h)/DT;

    b[0]/=a[0];
    d[0]/=a[0];
    for(l=1;l<nol;l++){
      a[l]-=c[l]*b[l-1];
      d[l]-=c[l]*d[l-1];
      d[l]/=a[l];
      b[l]/=a[l];
    }
    for(l=nol-2;l>=0;l--){
      d[l]-=b[l]*d[l+1];
    }

    //Update P    
    for(l=0;l<nol;l++){
      ceR_h[n].pn[l]=ceR_h[n].p[l]+d[l];
    }
    ceR_h[n].pn[nol-1]=ceR_h[n].rn[nol-1]*GV*(ceG_h[n].dz[nol-1]+ceR_h[n].hn);
    
  }

  //u,v
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      for(e=0;e<3;e++){
	if(suG_h[ceC_h[n].su[e]].bd==1){
	  dq[e]=ceR_h[ceC_h[n].ce[e]].pn[l]-ceR_h[n].pn[l];
	}
	else{
	  dq[e]=ceR_h[ceC_h[n].ce[e]].pn[l]-ceR_h[n].pn[l];
	}
	
      }
      
      f1=ceG_h[n].dx[0]*dq[0]+ceG_h[n].dx[1]*dq[1]+ceG_h[n].dx[2]*dq[2];
      f2=ceG_h[n].dy[0]*dq[0]+ceG_h[n].dy[1]*dq[1]+ceG_h[n].dy[2]*dq[2];
      
      ceR_h[n].un[l]-=DT/ceR_h[n].rn[l]*(ceG_h[n].dl[1]*f1-ceG_h[n].dl[2]*f2)*ceG_h[n].dl[3]/ceG_h[n].ds;
      ceR_h[n].vn[l]-=DT/ceR_h[n].rn[l]*(ceG_h[n].dl[0]*f2-ceG_h[n].dl[2]*f1)*ceG_h[n].dl[3]/ceG_h[n].ds;
    }
  }

  //w
  for(n=0;n<noe;n++){
    for(l=1;l<nol-1;l++){
      ceR_h[n].wn[l]-=((DATA_TYPE)0.5)*DT/ceR_h[n].rn[l]*(ceR_h[n].pn[l+1]-ceR_h[n].pn[l-1])/ceG_h[n].dz[l];
    }
    ceR_h[n].wn[0]=((DATA_TYPE)0.0);
    ceR_h[n].wn[nol-1]=ceR_h[n].wn[nol-2];
  }

  return;    
}
