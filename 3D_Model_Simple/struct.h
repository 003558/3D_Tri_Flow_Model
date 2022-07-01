// -*-c-*-
/**************************

 構造体の定義

 Last Update: 2017.05.16

***************************/

//////////////////
//Connectivity用の構造体
//////////////////

struct surf_conn //Surface
{
  int ve[2]; //両端のVertex番号
  int ce[2]; //Surfaceを含むCell番号
};

struct cell_conn //Cell
{
  int ve[4]; //Cellに含まれるVertex番号
  int su[3]; //Cellに含まれるSurface番号
  int ce[3]; //Cellに隣接するCell番号
};

//////////////////
//Geometry用の構造体
//////////////////
struct vert_geom //Vertex
{
  int ceN;
  int ce[20];

  DATA_TYPE x; //x座標値
  DATA_TYPE y; //y座標値
  DATA_TYPE z[NOZ+2];//z座標
  DATA_TYPE h;

  int bb[NOZ+2];
  DATA_TYPE bz;
};

struct surf_geom //Surface
{
  int ns;
  int bd[NOZ]; //数値流束の計算方法選択用Flag[-1:通常,0:壁境界,1:堤防境界,2:建物境界]
  DATA_TYPE xs; //Surface中心のx座標
  DATA_TYPE ys; //Surface中心のy座標
  DATA_TYPE zs[NOZ+1]; //Surface中心のz座標
  DATA_TYPE nnx; //垂直右向き単位ベクトルのx方向成分(m)
  DATA_TYPE nny; //垂直右向き単位ベクトルのy方向成分(m)
  DATA_TYPE nnl; //Surface長さ(m)

  DATA_TYPE ds;//隣接するセルの面積和
  DATA_TYPE dz[NOZ];
};

struct cell_geom //Cell
{
  int bd[NOZ];
  int ns;
  int bb[NOZ];
  DATA_TYPE sgn[3]; //各境界におけるFluxの流入方向の符号
  DATA_TYPE xc; //Cell中心のx座標値
  DATA_TYPE yc; //Cell中心のy座標値
  DATA_TYPE zc[NOZ]; //Cell中心のz座標
  DATA_TYPE ds; //Cell表面積(m^2)(三角形の面積)
  DATA_TYPE dx[3];
  DATA_TYPE dy[3];
  DATA_TYPE dl[4];
  DATA_TYPE dz[NOZ];
  DATA_TYPE bz;
  DATA_TYPE bzdx,bzdy;
};

//////////////////
//Valiables用の構造体
//////////////////
struct surf
{
  DATA_TYPE ub[NOZ];
  DATA_TYPE vb[NOZ];
  DATA_TYPE rb[NOZ];
  DATA_TYPE pb[NOZ];
  DATA_TYPE fl[NOZ];//Advection Flux

  DATA_TYPE dqx;
  DATA_TYPE dqy;
};

struct cell //Cell
{
  DATA_TYPE h;//水面の高さ[m]
  DATA_TYPE u[NOZ];//流速のx方向成分(m/s)
  DATA_TYPE v[NOZ];//流速のy方向成分(m/s)
  DATA_TYPE w[NOZ];//流速のz方向成分(m/s)
  DATA_TYPE r[NOZ];//密度
  DATA_TYPE p[NOZ];//圧力
  DATA_TYPE k[NOZ];//乱流エネルギー
  DATA_TYPE e[NOZ];//乱流散逸率
  DATA_TYPE s[NOZ];

  DATA_TYPE D[NOZ];//Divergence
  DATA_TYPE Psi[NOZ];//圧力補正

  DATA_TYPE hm;
  DATA_TYPE hn;
  DATA_TYPE un[NOZ];
  DATA_TYPE vn[NOZ];
  DATA_TYPE wn[NOZ];
  DATA_TYPE rn[NOZ];
  DATA_TYPE pn[NOZ];
  DATA_TYPE kn[NOZ];
  DATA_TYPE en[NOZ];
  DATA_TYPE sn[NOZ];

  DATA_TYPE Eh[NOZ];//Horizontal Diffusion Coefficient
  DATA_TYPE Ev[NOZ];//Vertical Diffusion Coefficient

  DATA_TYPE Pr[NOZ];
  DATA_TYPE Gk[NOZ];

  DATA_TYPE qx[NOZ];//dq/dx
  DATA_TYPE qy[NOZ];//dq/dy
  DATA_TYPE qz[NOZ];//dq/dz

  DATA_TYPE dux[NOZ];
  DATA_TYPE duy[NOZ];
  DATA_TYPE duz[NOZ];
  DATA_TYPE dvx[NOZ];
  DATA_TYPE dvy[NOZ];
  DATA_TYPE dvz[NOZ];
  DATA_TYPE dwx[NOZ];
  DATA_TYPE dwy[NOZ];
  DATA_TYPE dwz[NOZ];
  DATA_TYPE drx[NOZ];
  DATA_TYPE dry[NOZ];
  DATA_TYPE drz[NOZ];

  DATA_TYPE dpx[NOZ];
  DATA_TYPE dpy[NOZ];
  DATA_TYPE ddpx[NOZ];
  DATA_TYPE ddpy[NOZ];

  DATA_TYPE fl[NOZ+1];

  DATA_TYPE h_Adv;//Advection term
  DATA_TYPE u_Adv[NOZ];//Advection term
  DATA_TYPE v_Adv[NOZ];//Advection term
  DATA_TYPE w_Adv[NOZ];//Advection term
  DATA_TYPE s_Adv[NOZ];//Advection term
  DATA_TYPE k_Adv[NOZ];//Advection term
  DATA_TYPE e_Adv[NOZ];//Advection term

  DATA_TYPE u_Dif[NOZ];//Diffusion term
  DATA_TYPE v_Dif[NOZ];//Diffusion term
  DATA_TYPE w_Dif[NOZ];//Diffusion term
  DATA_TYPE s_Dif[NOZ];//Diffusion term
  DATA_TYPE k_Dif[NOZ];//Diffusion term
  DATA_TYPE e_Dif[NOZ];//Diffusion term

  DATA_TYPE u_Adt[NOZ];
  DATA_TYPE v_Adt[NOZ];
  DATA_TYPE k_Adt[NOZ];
  DATA_TYPE e_Adt[NOZ];

  DATA_TYPE wb[NOZ+1];
  DATA_TYPE rb[NOZ+1];
};

//Boundary

struct bound
{
  DATA_TYPE Wx[NOB];
  DATA_TYPE Wy[NOB];

  DATA_TYPE Ux;
  DATA_TYPE Uy;

  DATA_TYPE WL[NOB];
  DATA_TYPE H;
};
