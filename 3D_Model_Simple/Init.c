#pragma warning(disable: 4996)
// // -*- c -*-
/**************************

 初期条件の設定

 Last Update: 2017.05.16

***************************/
////////////////
//Init
////////////////
#include "public.h"
void Init
(
 int nol, int non, int nos, int noe,
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound *Bc
 )
{
  static int l,n;
  static DATA_TYPE Sp_Scale_Turb;
  static DATA_TYPE epse;
  FILE *fin;
  static DATA_TYPE Rad,Wind,WL;
	static DATA_TYPE x_naka, y_naka, x_yona, y_yona;
	static DATA_TYPE a_0, a_1, b_naka, b_yona;
	static DATA_TYPE x_t_naka, y_t_naka, x_t_yona, y_t_yona;
	static DATA_TYPE z1_naka, z2_naka, z3_naka, z4_naka, z5_naka;
	static DATA_TYPE z1_yona, z2_yona, z3_yona, z4_yona, z5_yona;
	static DATA_TYPE s1_naka, s2_naka, s3_naka, s4_naka, s5_naka;
	static DATA_TYPE s1_yona, s2_yona, s3_yona, s4_yona, s5_yona;
	static DATA_TYPE d_naka, d_yona;
	static DATA_TYPE w_naka, w_yona;
	static DATA_TYPE z1_c, z2_c, z3_c, z4_c, z5_c;
	static DATA_TYPE s1_c, s2_c, s3_c, s4_c, s5_c;

  printf("# Init Start ... ");

  //r
	
	//観測所平面座標
  x_naka = 93003.47;
	y_naka = -58717.81;
	x_yona = 102930.45;
	y_yona = -62431.23;

	//観測高(H.P.m)
	z1_naka = 0.025;
	z2_naka = -0.475;
	z3_naka = -2.875;
	z4_naka = -4.975;
	z5_naka = -6.000;

	z1_yona = 0.018;
	z2_yona = -0.484;
	z3_yona = -1.584;
	z4_yona = -2.684;
	z5_yona = -3.700;

	//塩分(psu)
	s1_naka = 7.2;
	s2_naka = 7.3;
	s3_naka = 20.0;
	s4_naka = 29.9;
	s5_naka = 30.2;

	s1_yona = 7.8;
	s2_yona = 9.0;
	s3_yona = 12.4;
	s4_yona = 15.4;
	s5_yona = 18.6;

	//観測所平面位置を通る直線の式
	a_0 = (y_yona - y_naka) / (x_yona - x_naka); //湖心-米子湾の傾き
	a_1 = -1.0 / a_0; //a_0の垂線の傾き
	b_naka = y_naka - a_1 * x_naka; //湖心を通る直線の切片
	b_yona = y_yona - a_1 * x_yona; //米子湾を通る直線の切片

	//初期塩分分布設定
  for(n=0;n<noe;n++){

		x_t_naka = (b_naka - ceG_h[n].yc) / a_1;
		y_t_naka = a_1 * ceG_h[n].xc + b_naka;

		x_t_yona = (b_yona - ceG_h[n].yc) / a_1;
		y_t_yona = a_1 * ceG_h[n].xc + b_yona;

		//湖心より大橋川側は湖心と同じ塩分分布
		if (ceG_h[n].xc <= x_t_naka &&  ceG_h[n].yc >= y_t_naka){
			for (l = 0; l < nol; l++){
				if (ceG_h[n].zc[l] <= z5_naka + ZZ){
					ceR_h[n].s[l] = s5_naka;
				}
				else if (ceG_h[n].zc[l] <= z4_naka + ZZ) {
					ceR_h[n].s[l] = (s5_naka - s4_naka) / (z5_naka - z4_naka) * (ceG_h[n].zc[l] - ZZ - z4_naka) + s4_naka;
				}
				else if (ceG_h[n].zc[l] <= z3_naka + ZZ) {
					ceR_h[n].s[l] = (s4_naka - s3_naka) / (z4_naka - z3_naka) * (ceG_h[n].zc[l] - ZZ - z3_naka) + s3_naka;
				}
				else if (ceG_h[n].zc[l] <= z2_naka + ZZ) {
					ceR_h[n].s[l] = (s3_naka - s2_naka) / (z3_naka - z2_naka) * (ceG_h[n].zc[l] - ZZ - z2_naka) + s2_naka;
				}
				else{
					ceR_h[n].s[l] = (s2_naka - s1_naka) / (z2_naka - z1_naka) * (ceG_h[n].zc[l] - ZZ - z1_naka) + s1_naka;
				}
			}
		}

		//米子湾より湾奥は米子湾と同じ塩分分布
		else if (ceG_h[n].xc >= x_t_yona &&  ceG_h[n].yc <= y_t_yona){
			for (l = 0; l < nol; l++){
				if (ceG_h[n].zc[l] <= z5_naka + ZZ){
					ceR_h[n].s[l] = s5_naka;
				}
				else if (ceG_h[n].zc[l] <= z5_yona + ZZ){
					ceR_h[n].s[l] = (s5_naka - s5_yona) / (z5_naka - z5_yona) * (ceG_h[n].zc[l] - ZZ - z5_yona) + s5_yona;
				}
				else if (ceG_h[n].zc[l] <= z4_yona + ZZ) {
					ceR_h[n].s[l] = (s5_yona - s4_yona) / (z5_yona - z4_yona) * (ceG_h[n].zc[l] - ZZ - z4_yona) + s4_yona;
				}
				else if (ceG_h[n].zc[l] <= z3_yona + ZZ) {
					ceR_h[n].s[l] = (s4_yona - s3_yona) / (z4_yona - z3_yona) * (ceG_h[n].zc[l] - ZZ - z3_yona) + s3_yona;
				}
				else if (ceG_h[n].zc[l] <= z2_yona + ZZ) {
					ceR_h[n].s[l] = (s3_yona - s2_yona) / (z3_yona - z2_yona) * (ceG_h[n].zc[l] - ZZ - z2_yona) + s2_yona;
				}
				else{
					ceR_h[n].s[l] = (s2_yona - s1_yona) / (z2_yona - z1_yona) * (ceG_h[n].zc[l] - ZZ - z1_yona) + s1_yona;
				}
			}
		}

		//湖心～米子湾の間は線形補間
		else{
			d_naka = sqrt(POW(ceG_h[n].xc - x_naka, 2) + POW(ceG_h[n].yc - y_naka, 2));
			d_yona = sqrt(POW(ceG_h[n].xc - x_yona, 2) + POW(ceG_h[n].yc - y_yona, 2));
			
			w_naka = 1.0 / d_naka;
			w_yona = 1.0 / d_yona;
			
			//観測高の重み付き平均
			z1_c = (w_naka * z1_naka + w_yona * z1_yona) / (w_naka + w_yona);
			z2_c = (w_naka * z2_naka + w_yona * z2_yona) / (w_naka + w_yona);
			z3_c = (w_naka * z3_naka + w_yona * z3_yona) / (w_naka + w_yona);
			z4_c = (w_naka * z4_naka + w_yona * z4_yona) / (w_naka + w_yona);
			z5_c = (w_naka * z5_naka + w_yona * z5_yona) / (w_naka + w_yona);

			//塩分の重み付き平均
			s1_c = (w_naka * s1_naka + w_yona * s1_yona) / (w_naka + w_yona);
			s2_c = (w_naka * s2_naka + w_yona * s2_yona) / (w_naka + w_yona);
			s3_c = (w_naka * s3_naka + w_yona * s3_yona) / (w_naka + w_yona);
			s4_c = (w_naka * s4_naka + w_yona * s4_yona) / (w_naka + w_yona);
			s5_c = (w_naka * s5_naka + w_yona * s5_yona) / (w_naka + w_yona);

			for (l = 0; l < nol; l++){
				if (ceG_h[n].zc[l] <= z5_naka + ZZ){
					ceR_h[n].s[l] = s5_naka;
				}
				else if (ceG_h[n].zc[l] <= z5_c + ZZ){
					ceR_h[n].s[l] = (s5_naka - s5_c) / (z5_naka - z5_c) * (ceG_h[n].zc[l] - ZZ - z5_c) + s5_c;
				}
				else if (ceG_h[n].zc[l] <= z4_c + ZZ) {
					ceR_h[n].s[l] = (s5_c - s4_c) / (z5_c - z4_c) * (ceG_h[n].zc[l] - ZZ - z4_c) + s4_c;
				}
				else if (ceG_h[n].zc[l] <= z3_c + ZZ) {
					ceR_h[n].s[l] = (s4_c - s3_c) / (z4_c - z3_c) * (ceG_h[n].zc[l] - ZZ - z3_c) + s3_c;
				}
				else if (ceG_h[n].zc[l] <= z2_yona + ZZ) {
					ceR_h[n].s[l] = (s3_c - s2_c) / (z3_c - z2_c) * (ceG_h[n].zc[l] - ZZ - z2_c) + s2_c;
				}
				else{
					ceR_h[n].s[l] = (s2_c - s1_c) / (z2_c - z1_c) * (ceG_h[n].zc[l] - ZZ - z1_c) + s1_c;
				}
			}

		}
  }

  //h
  for(n=0;n<noe;n++){
    ceR_h[n].h=((DATA_TYPE)0.2);
  }

  for(n=0;n<non;n++){
    veG_h[n].h=((DATA_TYPE)0.2);
    //veG_h[n].h/=veG_h[n].ceN;
  }
  for(n=0;n<non;n++){
    veG_h[n].z[nol+1]=veG_h[n].z[nol]+veG_h[n].h;
  }

  //u,v,Eh
  for(l=0;l<nol;l++){
    for(n=0;n<noe;n++){
      ceR_h[n].u[l]=((DATA_TYPE)0.0);
      ceR_h[n].v[l]=((DATA_TYPE)0.0);
      ceR_h[n].w[l]=((DATA_TYPE)0.0);
      ceR_h[n].k[l]=((DATA_TYPE)0.0);
      ceR_h[n].e[l]=((DATA_TYPE)0.0);
      ceR_h[n].Eh[l]=((DATA_TYPE)0.01)*POW(D_h,((DATA_TYPE)(4.0/3.0)));
    }
  }

  for(n=0;n<noe;n++){
    Sp_Scale_Turb=ceR_h[n].h+((DATA_TYPE)0.5)*ceG_h[n].dz[nol-1]+ceG_h[n].zc[nol-1]-ceG_h[n].zc[ceG_h[n].ns];
    for(l=ceG_h[n].ns;l<nol;l++){
      if(ceR_h[n].k[l]<=((DATA_TYPE)7.6e-26)){
	ceR_h[n].k[l]=((DATA_TYPE)7.6e-26);
	ceR_h[n].e[l]=POW(ceR_h[n].k[l]*sqrt(C_mu),((DATA_TYPE)1.5))/Sp_Scale_Turb;
      }
      epse=POW(sqrt(C_mu)*ceR_h[n].k[l],((DATA_TYPE)1.5))/Sp_Scale_Turb;
      if(ceR_h[n].e[l]<epse){
	ceR_h[n].e[l]=epse;
      }
    }
  }

  for(n=0;n<noe;n++){
    ceR_h[n].r[nol-1]=Density_Value(ceR_h[n].s[nol-1],((DATA_TYPE)20.0),((DATA_TYPE)0.0));
    ceR_h[n].p[nol-1]=ceR_h[n].r[nol-1]*GV*(ceR_h[n].h+((DATA_TYPE)0.5)*ceG_h[n].dz[nol-1]);
    for(l=nol-1;l>=0;l--){
      ceR_h[n].r[l]=Density_Value(ceR_h[n].s[l],((DATA_TYPE)20.0),ceR_h[n].p[l+1]);
      ceR_h[n].p[l]+=ceR_h[n].p[l+1]+GV*((DATA_TYPE)0.5)*(ceR_h[n].r[l]*ceG_h[n].dz[l]+ceR_h[n].r[l+1]*ceG_h[n].dz[l+1]);
    }

  }
  printf("End\n");

  fin=fopen("Data/Wind.dat","r");
  for(n=0;n<NOB;n++){
    fscanf(fin,"%lf %lf",&Rad,&Wind);
    Rad=((DATA_TYPE)180.0)-Rad;
    if(Rad<((DATA_TYPE)0.0)){
      Rad=((DATA_TYPE)360.0)+Rad;
    }
    (*Bc).Wx[n]=Wind*cos(Rad/((DATA_TYPE)180.0)*M_PI);
    (*Bc).Wy[n]=Wind*sin(Rad/((DATA_TYPE)180.0)*M_PI);
  }
  fclose(fin);

  fin=fopen("Data/WL.dat","r");
  for(n=0;n<NOB;n++){
    fscanf(fin,"%lf",&WL);
    (*Bc).WL[n]=WL;
  }
  fclose(fin);

  return;
}
