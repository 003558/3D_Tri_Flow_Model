// -*- c -*-
/**************************

 Main Program

 Last Update: 2017.05.16

***************************/

#include "public.h"

int main(int argc, char *argv[])
{
  int nol,non,nos,noe;

  int step;
  int step_o;
  DATA_TYPE step_t;

  //Vertex
  struct vert_geom *veG_h;

  //Surface
  struct surf_conn *suC_h;
  struct surf_geom *suG_h;
  struct surf *suR_h;

  //Cell
  struct cell_conn *ceC_h;
  struct cell_geom *ceG_h;
  struct cell *ceR_h;

  //Boundary condition
  int nt;
  struct bound Bc;

  /* FILE *fout; */
  /* int nn,ln; */

  /* FILE *fout1; */
  /* int nn1,ln1; */

  /* FILE *fout2; */
  /* int nn2,ln2; */


  ////////////////////////
  // Pre Calculation
  ////////////////////////

  printf("###############################################\n");

  Read
    (
     &nol,&non,&nos,&noe,
     &veG_h,
     &suC_h,&suG_h,&suR_h,
     &ceC_h,&ceG_h,&ceR_h
     );
    
  Geom
    (
     nol,non,nos,noe,
     veG_h,
     suC_h,suG_h,
     ceC_h,ceG_h 
     );
  
  Init
    (
     nol,non,nos,noe,
     veG_h,
     suC_h,suG_h,suR_h,
     ceC_h,ceG_h,ceR_h,
     &Bc
     );

  printf("###############################################\n");


  
  ////////////////////////
  // Main Calculation
  ////////////////////////

  step=0;
  step_o=0;
  step_t=0.0;

  Cal_Info
    (
     step,step_o,step_t,
     nol,nos,noe,
     veG_h,
     suC_h,suG_h,suR_h,
     ceC_h,ceG_h,ceR_h
     );
  
  Dat_Tecp
    (
     step_o, 
     nol, non, nos, noe,
     veG_h,
     suC_h,suG_h,suR_h,
     ceC_h,ceG_h,ceR_h
     );
  
  /* fout=fopen("Data/Check_Value_20.csv","w"); */
  /* fprintf(fout,"step,u,v,w,p,D,R,S,k,e,Ev,Pr,Gk,s_Adv,s_Dif\n"); */
  /* nn=1478; */
  /* ln=20; */

  /* fout1=fopen("Data/Check_Value_21.csv","w"); */
  /* fprintf(fout1,"step,u,v,w,p,D,R,S,k,e,Ev,Pr,Gk,s_Adv,s_Dif\n"); */
  /* nn1=1478; */
  /* ln1=21; */

  /* fout2=fopen("Data/Check_Value_22.csv","w"); */
  /* fprintf(fout2,"step,u,v,w,p,D,R,S,k,e,Ev,Pr,Gk,s_Adv,s_Dif\n"); */
  /* nn2=1478; */
  /* ln2=22; */
  
  nt=0;

  for(step=1;step<=STEP_MAX;step++){
    step_t=step_t+DT;

    if(step_t>DDT*nt){
      nt++;
      printf("nt=%d\n",nt);
    }

    Bc.Ux=(Bc.Wx[nt]-Bc.Wx[nt-1])/DDT*(((DATA_TYPE)step_t-DDT*(nt-1)))+Bc.Wx[nt-1];
    Bc.Uy=(Bc.Wy[nt]-Bc.Wy[nt-1])/DDT*(((DATA_TYPE)step_t-DDT*(nt-1)))+Bc.Wy[nt-1];
    Bc.H=(Bc.WL[nt]-Bc.WL[nt-1])/DDT*(((DATA_TYPE)step_t-DDT*(nt-1)))+Bc.WL[nt-1];

	//printf("Set_Dependent_Value Start\n");
    Set_Dependent_Value
      (
       nol, nos,noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h,
       Bc
       );
    
	//printf("Advection Start\n");
	//////////////////
    // Advection
    //////////////////
    h_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    u_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    v_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    w_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    s_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    k_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    e_Advection
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );



	//printf("Diffusion Start\n");
	//////////////////
    // Diffusion
    //////////////////
    u_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    v_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    w_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    s_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    k_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );

    e_Diffusion
      (
       nol, nos, noe,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h
       );


	//printf("Pressure Start\n");
	Pressure
      (
       nol,nos,noe,
       veG_h,
       suC_h,suG_h,suR_h,
       ceC_h,ceG_h,ceR_h,
       Bc
       );

	//printf("Sync Start\n");
	Sync
      (
       nol,non,noe,
       veG_h,
       ceR_h
       );
	//printf("Sync End\n");


    if(!(step%DAT_OUT)){

      step_o=step_o+1;

	  //printf("Cal_Info Start\n");
	  Cal_Info
  	(
  	 step,step_o,step_t,
	 nol,nos,noe,
	 veG_h,
	 suC_h,suG_h,suR_h,
	 ceC_h,ceG_h,ceR_h
  	 );

	  //printf("Dat_Tecp Start\n");
	  Dat_Tecp
  	(
  	 step_o,
  	 nol, non, nos, noe,
  	 veG_h,
  	 suC_h,suG_h,suR_h,
  	 ceC_h,ceG_h,ceR_h
  	 );
    }

    ///
    /* fprintf(fout,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",step,ceR_h[nn].u[ln],ceR_h[nn].v[ln],ceR_h[nn].w[ln],ceR_h[nn].p[ln],ceR_h[nn].D[ln],ceR_h[nn].r[ln],ceR_h[nn].s[ln],ceR_h[nn].k[ln],ceR_h[nn].e[ln],ceR_h[nn].Ev[ln],ceR_h[nn].Pr[ln],ceR_h[nn].Gk[ln],ceR_h[nn].s_Adv[ln],ceR_h[nn].s_Dif[ln]); */

    /* fprintf(fout1,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",step,ceR_h[nn1].u[ln1],ceR_h[nn1].v[ln1],ceR_h[nn1].w[ln1],ceR_h[nn1].p[ln1],ceR_h[nn1].D[ln1],ceR_h[nn1].r[ln1],ceR_h[nn1].s[ln1],ceR_h[nn1].k[ln1],ceR_h[nn1].e[ln1],ceR_h[nn1].Ev[ln1],ceR_h[nn1].Pr[ln1],ceR_h[nn1].Gk[ln1],ceR_h[nn1].s_Adv[ln1],ceR_h[nn1].s_Dif[ln1]); */

    /* fprintf(fout2,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",step,ceR_h[nn2].u[ln2],ceR_h[nn2].v[ln2],ceR_h[nn2].w[ln2],ceR_h[nn2].p[ln2],ceR_h[nn2].D[ln2],ceR_h[nn2].r[ln2],ceR_h[nn2].s[ln2],ceR_h[nn2].k[ln2],ceR_h[nn2].e[ln2],ceR_h[nn2].Ev[ln2],ceR_h[nn2].Pr[ln2],ceR_h[nn2].Gk[ln2],ceR_h[nn2].s_Adv[ln2],ceR_h[nn2].s_Dif[ln2]); */

  }
  /* fclose(fout); */
  /* fclose(fout1); */
  /* fclose(fout2); */

  return 0;
}
