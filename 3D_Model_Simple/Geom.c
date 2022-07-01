#pragma warning(disable: 4996)
// -*- c -*-
/**************************

 Geometryデータの設定

 Last Update: 2017.05.16

***************************/

#include "public.h"

void Geom 
(
 int nol, int non, int nos, int noe, 
 struct vert_geom *veG_h, 
 struct surf_conn *suC_h, struct surf_geom *suG_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h
 )
{
  static int l,n,e;
  static DATA_TYPE xs,ys;
  static DATA_TYPE nnx,nny;
  static DATA_TYPE xp[3],yp[3];
  static DATA_TYPE a,b;
  FILE *fin;
  static int ve[3];
  static DATA_TYPE ax,ay,az,bx,by,bz,A,B,C;

  printf("# Geom Start ... ");

  //Vertex Geometry
  for(l=0;l<=nol;l++){
    for(n=0;n<non;n++){
      veG_h[n].z[l]=(Zmax-Zmin)/((DATA_TYPE)nol)*l+Zmin+ZZ;
    }
  }

  fin=fopen("Data/Vert_Bottom.dat","r");
  //fin=fopen("Data/Vert_Bottom_case14.dat","r");
  for(n=0;n<non;n++){
    fscanf(fin,"%lf",&veG_h[n].bz);
    if(veG_h[n].bz>=-((DATA_TYPE)1.0)){
      veG_h[n].bz=((DATA_TYPE)-1.0);
    }
    veG_h[n].bz+=ZZ;
  }
  fclose(fin);

  //Surface Geometry
  for(n=0;n<nos;n++){
    xs=((DATA_TYPE)0.5)*(veG_h[suC_h[n].ve[0]].x+veG_h[suC_h[n].ve[1]].x);
    ys=((DATA_TYPE)0.5)*(veG_h[suC_h[n].ve[0]].y+veG_h[suC_h[n].ve[1]].y);
    nnx=veG_h[suC_h[n].ve[1]].y-veG_h[suC_h[n].ve[0]].y;
    nny=veG_h[suC_h[n].ve[0]].x-veG_h[suC_h[n].ve[1]].x;

    suG_h[n].xs=xs;
    suG_h[n].ys=ys;
    suG_h[n].nnl=sqrt(nnx*nnx+nny*nny);
    suG_h[n].nnx=nnx/suG_h[n].nnl;
    suG_h[n].nny=nny/suG_h[n].nnl;
  }

  for(l=0;l<=nol;l++){
    for(n=0;n<non;n++){
      suG_h[n].zs[l]=(Zmax-Zmin)/((DATA_TYPE)nol)*l+Zmin+ZZ;
    }
  }
  
  
  //Cell Geometry
  for(n=0;n<noe;n++){
    ceC_h[n].ve[3]=ceC_h[n].ve[0];

    for(e=0;e<3;e++){
      xp[e]=veG_h[ceC_h[n].ve[e]].x;
      yp[e]=veG_h[ceC_h[n].ve[e]].y;
    }

    ceG_h[n].xc=(xp[0]+xp[1]+xp[2])/((DATA_TYPE)3.0);
    ceG_h[n].yc=(yp[0]+yp[1]+yp[2])/((DATA_TYPE)3.0);
    ceG_h[n].ds=((DATA_TYPE)0.5)*((xp[0]-xp[1])*(yp[0]-yp[2])-(xp[0]-xp[2])*(yp[0]-yp[1]));
  }

  

  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      ceG_h[n].zc[l]=(Zmax-Zmin)/((DATA_TYPE)nol)*l+Zmin+ZZ;
    }
  }

  for(n=0;n<noe;n++){
    ceG_h[n].bz=(veG_h[ceC_h[n].ve[0]].bz+veG_h[ceC_h[n].ve[1]].bz+veG_h[ceC_h[n].ve[2]].bz)/((DATA_TYPE)3.0);
  }


  for(n=0;n<noe;n++){
    for(l=1;l<nol-1;l++){
      ceG_h[n].dz[l]=ceG_h[n].zc[l+1]-ceG_h[n].zc[l];
    }
    ceG_h[n].dz[0]=ceG_h[n].dz[1];
    ceG_h[n].dz[nol-1]=ceG_h[n].dz[nol-2];
  }

  //bottom
  for(n=0;n<non;n++){
    for(l=0;l<nol;l++){
      veG_h[n].bb[l]=0;
    }
    for(l=0;l<nol;l++){
      if(veG_h[n].z[l]<=veG_h[n].bz&&veG_h[n].z[l+1]>veG_h[n].bz){
	veG_h[n].bb[l]=1;
      }
    }
  }

  for(n=0;n<noe;n++){
    for(l=0;l<nol;l++){
      if(ceG_h[n].zc[l]<=ceG_h[n].bz){
	ceG_h[n].bb[l]=1;
	for(e=0;e<3;e++){
	  suG_h[ceC_h[n].su[e]].bd[l]=0;
	}
      }
      else{
	ceG_h[n].bb[l]=0;
      }
    }
  }

  for(n=0;n<noe;n++){
    for(l=0;l<nol;l++){
      if(ceG_h[n].bb[l]==0){
	ceG_h[n].ns=l;
	break;
      }
    }
  }

  for(n=0;n<nos;n++){
    suG_h[n].ns=MIN(ceG_h[suC_h[n].ce[0]].ns,ceG_h[suC_h[n].ce[1]].ns);
  }


  for(n=0;n<noe;n++){
    for(e=0;e<3;e++){
      ve[e]=ceC_h[n].ve[e];
    }
    ax=(veG_h[ve[1]].x-veG_h[ve[0]].x);
    ay=(veG_h[ve[1]].y-veG_h[ve[0]].y);
    az=(veG_h[ve[1]].bz-veG_h[ve[0]].bz);
    bx=(veG_h[ve[2]].x-veG_h[ve[0]].x);
    by=(veG_h[ve[2]].y-veG_h[ve[0]].y);
    bz=(veG_h[ve[2]].bz-veG_h[ve[0]].bz);

    A=ay*bz-by*az;
    B=az*bx-bz*ax;
    C=ax*by-bx*ay;

    ceG_h[n].bzdx=-A/C;
    ceG_h[n].bzdy=-B/C;
  }

  for(n=0;n<noe;n++){
    for(e=0;e<3;e++){
      if(ceC_h[n].ce[e]==n){
	if(veG_h[suC_h[ceC_h[n].su[e]].ve[0]].x==veG_h[suC_h[ceC_h[n].su[e]].ve[1]].x){
	  ceG_h[n].dx[e]=((DATA_TYPE)2.0)*(suG_h[ceC_h[n].su[e]].xs-ceG_h[n].xc);
	  ceG_h[n].dy[e]=((DATA_TYPE)0.0);
	}
	else if(veG_h[suC_h[ceC_h[n].su[e]].ve[0]].y==veG_h[suC_h[ceC_h[n].su[e]].ve[1]].y){
	  ceG_h[n].dx[e]=((DATA_TYPE)0.0);
	  ceG_h[n].dy[e]=((DATA_TYPE)2.0)*(suG_h[ceC_h[n].su[e]].ys-ceG_h[n].yc);
	}
	else{
	  a=(veG_h[suC_h[ceC_h[n].su[e]].ve[1]].y-veG_h[suC_h[ceC_h[n].su[e]].ve[0]].y)/(veG_h[suC_h[ceC_h[n].su[e]].ve[1]].x-veG_h[suC_h[ceC_h[n].su[e]].ve[0]].x);
	  b=-a*veG_h[suC_h[ceC_h[n].su[e]].ve[0]].x+veG_h[suC_h[ceC_h[n].su[e]].ve[0]].y;
	  
	  ceG_h[n].dx[e]=(a*(ceG_h[n].yc-b)-a*a*ceG_h[n].xc)/(a*a+((DATA_TYPE)1.0));
	  ceG_h[n].dy[e]=(a*ceG_h[n].xc-ceG_h[n].yc-a*b)/(a*a+((DATA_TYPE)1.0))+b;
	}
      }
      else{
	ceG_h[n].dx[e]=ceG_h[ceC_h[n].ce[e]].xc-ceG_h[n].xc;
	ceG_h[n].dy[e]=ceG_h[ceC_h[n].ce[e]].yc-ceG_h[n].yc;
      }
    }

    ceG_h[n].dl[0]=ceG_h[n].dx[0]*ceG_h[n].dx[0]+ceG_h[n].dx[1]*ceG_h[n].dx[1]+ceG_h[n].dx[2]*ceG_h[n].dx[2];
    ceG_h[n].dl[1]=ceG_h[n].dy[0]*ceG_h[n].dy[0]+ceG_h[n].dy[1]*ceG_h[n].dy[1]+ceG_h[n].dy[2]*ceG_h[n].dy[2];
    ceG_h[n].dl[2]=ceG_h[n].dx[0]*ceG_h[n].dy[0]+ceG_h[n].dx[1]*ceG_h[n].dy[1]+ceG_h[n].dx[2]*ceG_h[n].dy[2];
    ceG_h[n].dl[3]=FABS(((DATA_TYPE)1.0)/(ceG_h[n].dl[0]*ceG_h[n].dl[1]-ceG_h[n].dl[2]*ceG_h[n].dl[2]));
  }


  for(n=0;n<nos;n++){
    suG_h[n].ds=ceG_h[suC_h[n].ce[0]].ds+ceG_h[suC_h[n].ce[1]].ds;
  }

  printf("End\n");

  return;
}
