#pragma warning(disable: 4996)
// -*- c -*-
/**************************

 三角形格子データファイル"Unst.dat"の読込み

 Last Update: 2017.05.16

***************************/
#include "public.h"
void Read
(
 int *nol, int *non, int *nos, int *noe,
 struct vert_geom **veG_h, 
 struct surf_conn **suC_h, struct surf_geom **suG_h, struct surf **suR_h,
 struct cell_conn **ceC_h, struct cell_geom **ceG_h, struct cell **ceR_h
 )
{
  int n,e,i,l;
  int bd;
  FILE *fin;

  printf("# Read Start ... ");

  (*nol)=NOZ;
  
  fin=fopen("Data/Unst.dat","r");
  fscanf(fin,"%d %d %d",non,nos,noe);

  //*veG_h=(struct vert_geom *)malloc(sizeof(struct vert_geom)*(*non));
  //*suC_h=(struct surf_conn *)malloc(sizeof(struct surf_conn)*(*nos));
  //*suG_h=(struct surf_geom *)malloc(sizeof(struct surf_geom)*(*nos));
  //*suR_h=(struct surf *)malloc(sizeof(struct surf)*(*nos));
  //*ceC_h=(struct cell_conn *)malloc(sizeof(struct cell_conn)*(*noe));
  //*ceG_h=(struct cell_geom *)malloc(sizeof(struct cell_geom)*(*noe));
  //*ceR_h=(struct cell *)malloc(sizeof(struct cell)*(*noe));
  *veG_h = (struct vert_geom *)calloc((*non), sizeof(struct vert_geom));
  *suC_h = (struct surf_conn *)calloc((*nos), sizeof(struct surf_conn));
  *suG_h = (struct surf_geom *)calloc((*nos), sizeof(struct surf_geom));
  *suR_h = (struct surf *)calloc((*nos), sizeof(struct surf));
  *ceC_h = (struct cell_conn *)calloc((*noe), sizeof(struct cell_conn));
  *ceG_h = (struct cell_geom *)calloc((*noe), sizeof(struct cell_geom));
  *ceR_h = (struct cell *)calloc((*noe), sizeof(struct cell));

  for(n=0;n<(*non);n++){
    fscanf(fin,"%lf %lf",&(*veG_h)[n].x,&(*veG_h)[n].y);
  }

  for(n=0;n<(*nos);n++){
    fscanf(fin,"%d %d %d %d %d",
	   &bd,
	   &(*suC_h)[n].ve[0],&(*suC_h)[n].ve[1],
	   &(*suC_h)[n].ce[0],&(*suC_h)[n].ce[1]
	   );

    for(l=0;l<(*nol);l++){
      (*suG_h)[n].bd[l]=bd;
    }

  }

  for(n=0;n<(*noe);n++){
    fscanf(fin,"%d %d %d %d %d %d %d %d %d",
	   &(*ceC_h)[n].ve[0],&(*ceC_h)[n].ve[1],&(*ceC_h)[n].ve[2],
	   &(*ceC_h)[n].su[0],&(*ceC_h)[n].su[1],&(*ceC_h)[n].su[2],
	   &(*ceC_h)[n].ce[0],&(*ceC_h)[n].ce[1],&(*ceC_h)[n].ce[2]
	   );
  }
  fclose(fin);

  for(n=0;n<(*noe);n++){
    for(e=0;e<3;e++){
      if((*suC_h)[(*ceC_h)[n].su[e]].ce[0]==n){
	(*ceG_h)[n].sgn[e]=((DATA_TYPE)-1.0);
      }
      else{
	(*ceG_h)[n].sgn[e]=((DATA_TYPE)1.0);
      }
    }
  }
  
  for(n=0;n<(*noe);n++){
    for(e=0;e<3;e++){
      if((*ceC_h)[n].ce[e]==-1){
	(*ceC_h)[n].ce[e]=n;
      }
      if((*ceC_h)[n].su[e]==-1){
	printf("\n");
	printf("ceC_su %d %d ERROR\n",n,e);
      }
    }
  }

  for(n=0;n<(*nos);n++){
    if(((*suC_h)[n].ce[0]==-1)&&((*suC_h)[n].ce[1]!=-1)){
      (*suC_h)[n].ce[0]=(*suC_h)[n].ce[1];
    }
    if(((*suC_h)[n].ce[0]!=-1)&&((*suC_h)[n].ce[1]==-1)){
      (*suC_h)[n].ce[1]=(*suC_h)[n].ce[0];
    }
    if(((*suC_h)[n].ce[0]==-1)&&((*suC_h)[n].ce[1]==-1)){
      printf("\n");
      printf("suC_ce %d ERROR\n",n);
    }
  }



  for(n=0;n<(*noe);n++){
    for(i=0;i<3;i++){
      (*veG_h)[(*ceC_h)[n].ve[i]].ce[(*veG_h)[(*ceC_h)[n].ve[i]].ceN]=n;
      (*veG_h)[(*ceC_h)[n].ve[i]].ceN++;
    }
  }



  printf("End\n");

  return;

}

