#pragma warning(disable: 4996)
// -*- c -*-
#include "public.h"
////////////////
//Data
////////////////
void Dat_Tecp
(
 int step_o,
 int nol, int non, int nos, int noe,
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h 
 )
{
  static int l,n;
  static char f_name[50];
  static FILE *fout;
  static DATA_TYPE SS,CC;

  sprintf(f_name,"Output_1/t_%04d.plt",step_o);
  fout = fopen(f_name,"w");
  
  fprintf(fout,"TITLE = \"FE_Cell\"\n");
  fprintf(fout,"VARIABLES =\"X\", \"Y\",\"Z\", \"U\", \"V\",\"W\",\"P\",\"D\",\"R\",\"s\",\"k\",\"e\",\"Ev\", \"Pr\",\"Gk\"\n");
  fprintf(fout,"ZONE T=\"t_%04d\", N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",step_o,non*(nol+2),noe*(nol+1));
  fprintf(fout,"VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED,[7]=CELLCENTERED,[8]=CELLCENTERED,[9]=CELLCENTERED,[10]=CELLCENTERED,[11]=CELLCENTERED,[12]=CELLCENTERED,[13]=CELLCENTERED,[14]=CELLCENTERED,[15]=CELLCENTERED)\n");
  
  //1:x
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)veG_h[n].x);
    }
  }


  //2:y
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)veG_h[n].y);
    }
  }

  //3:z
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)veG_h[n].z[l]-ZZ);
    }
  }

  //4:u
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].u[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].u[nol-1]);
  }

  //5:v
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].v[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].v[nol-1]);
  }


  //6:w
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].w[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].w[nol-1]);
  }


  //7:p
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].p[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"0.0\n");
  }
  
  
  //8:D
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].D[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"0.0\n");
  }


 //9:Rho
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].r[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].r[nol-1]);
  }

 //10:s
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].s[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].s[nol-1]);
  }

 //11:k
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].k[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].k[nol-1]);
  }


 //12:e
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].e[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].e[nol-1]);
  }

 //13:Ev
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Ev[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Ev[nol-1]);
  }

 //14:Pr
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Pr[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Pr[nol-1]);
  }

 //15:Gk
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Gk[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Gk[nol-1]);
  }

  /* for(l=0;l<nol;l++){ */
  /*   for(n=0;n<noe;n++){ */
  /*     fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdx); */
  /*   } */
  /* } */
  /* for(n=0;n<noe;n++){ */
  /*   fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdx); */
  /* } */
  /* for(l=0;l<nol;l++){ */
  /*   for(n=0;n<noe;n++){ */
  /*     fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdy); */
  /*   } */
  /* } */
  /* for(n=0;n<noe;n++){ */
  /*   fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdy); */
  /* } */


  //Cell Connecitity
  for(l=0;l<=nol;l++){
    for(n=0;n<noe;n++){
      fprintf(fout,"%d %d %d %d ",non*l+ceC_h[n].ve[0]+1,non*l+ceC_h[n].ve[1]+1,non*l+ceC_h[n].ve[2]+1,non*l+ceC_h[n].ve[2]+1);
      fprintf(fout,"%d %d %d %d\n",non*(l+1)+ceC_h[n].ve[0]+1,non*(l+1)+ceC_h[n].ve[1]+1,non*(l+1)+ceC_h[n].ve[2]+1,non*(l+1)+ceC_h[n].ve[2]+1);

    }
  }
    
  fclose(fout);


  sprintf(f_name,"Output_2/t_%04d.plt",step_o);
  fout = fopen(f_name,"w");
  
  fprintf(fout,"TITLE = \"FE_Cell\"\n");
  fprintf(fout,"VARIABLES =\"X\", \"Y\",\"Z\", \"U\", \"V\",\"W\",\"P\",\"D\",\"R\",\"s\",\"k\",\"e\",\"Ev\", \"Pr\",\"Gk\",\"s_Adv\",\"s_Dif\"\n");
  fprintf(fout,"ZONE T=\"t_%04d\", N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n",step_o,non*(nol+2),noe*(nol+1));
  fprintf(fout,"VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED,[7]=CELLCENTERED,[8]=CELLCENTERED,[9]=CELLCENTERED,[10]=CELLCENTERED,[11]=CELLCENTERED,[12]=CELLCENTERED,[13]=CELLCENTERED,[14]=CELLCENTERED,[15]=CELLCENTERED,[16]=CELLCENTERED,[17]=CELLCENTERED)\n");

  CC=cos(-22.0/180.0*M_PI);
  SS=sin(-22.0/180.0*M_PI);

  //1:x1
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)(veG_h[n].x*CC+veG_h[n].y*SS));
    }
  }


  //2:y1
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)(-veG_h[n].x*SS+veG_h[n].y*CC));
    }
  }

  //3:z
  for(l=0;l<=(nol+1);l++){
    for(n=0;n<non;n++){
      fprintf(fout,"%lf\n",(double)veG_h[n].z[l]-ZZ);
    }
  }

  //4:u
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)(ceR_h[n].u[l]*CC+ceR_h[n].v[l]*SS));
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)(ceR_h[n].u[nol-1]*CC+ceR_h[n].v[nol-1]*SS));
  }

  //5:v
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)(-ceR_h[n].u[l]*SS+ceR_h[n].v[l]*CC));
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)(-ceR_h[n].u[nol-1]*SS+ceR_h[n].v[nol-1]*CC));
  }


  //6:w
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].w[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].w[nol-1]);
  }


  //7:p
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].p[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"0.0\n");
  }
  
  
  //8:D
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].D[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"0.0\n");
  }


 //9:Rho
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].r[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].r[nol-1]);
  }

 //10:s
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].s[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].s[nol-1]);
  }

 //11:k
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].k[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].k[nol-1]);
  }


 //12:e
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].e[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].e[nol-1]);
  }

 //13:Ev
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Ev[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Ev[nol-1]);
  }

 //14:Pr
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Pr[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Pr[nol-1]);
  }

 //15:Gk
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].Gk[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].Gk[nol-1]);
  }

 //15:s_Adv
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].s_Adv[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].s_Adv[nol-1]);
  }

 //15:s_Dif
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      if(l<ceG_h[n].ns){
	fprintf(fout,"0.0\n");
      }
      else{
	fprintf(fout,"%10.4e\n",(double)ceR_h[n].s_Dif[l]);
      }
    }
  }
  for(n=0;n<noe;n++){
    fprintf(fout,"%10.4e\n",(double)ceR_h[n].s_Dif[nol-1]);
  }

  /* for(l=0;l<nol;l++){ */
  /*   for(n=0;n<noe;n++){ */
  /*     fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdx); */
  /*   } */
  /* } */
  /* for(n=0;n<noe;n++){ */
  /*   fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdx); */
  /* } */
  /* for(l=0;l<nol;l++){ */
  /*   for(n=0;n<noe;n++){ */
  /*     fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdy); */
  /*   } */
  /* } */
  /* for(n=0;n<noe;n++){ */
  /*   fprintf(fout,"%10.4e\n",(double)ceG_h[n].bzdy); */
  /* } */


  //Cell Connecitity
  for(l=0;l<=nol;l++){
    for(n=0;n<noe;n++){
      fprintf(fout,"%d %d %d %d ",non*l+ceC_h[n].ve[0]+1,non*l+ceC_h[n].ve[1]+1,non*l+ceC_h[n].ve[2]+1,non*l+ceC_h[n].ve[2]+1);
      fprintf(fout,"%d %d %d %d\n",non*(l+1)+ceC_h[n].ve[0]+1,non*(l+1)+ceC_h[n].ve[1]+1,non*(l+1)+ceC_h[n].ve[2]+1,non*(l+1)+ceC_h[n].ve[2]+1);

    }
  }
    
  fclose(fout);

  return;
}
