// -*-c-*-
/**************************

 関数のプロトタイプ宣言

 Last Update: 2017.5.16

***************************/

////////////////////////
// Pre Calculation
////////////////////////
void Read
(
 int *nol, int *non, int *nos, int *noe,
 struct vert_geom **veG_h, 
 struct surf_conn **suC_h, struct surf_geom **suG_h, struct surf **suR_h,
 struct cell_conn **ceC_h, struct cell_geom **ceG_h, struct cell **ceR_h 
 );

//Geometry
void Geom 
(
 int nol, int non, int nos, int noe, 
 struct vert_geom *veG_h, 
 struct surf_conn *suC_h, struct surf_geom *suG_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h
 );

//Initial
void Init
(
 int nol, int non, int nos, int noe,
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound *Bc
 );


//Output Data 
void Dat_Tecp
(
 int step_o,
 int nol, int non, int nos, int noe,
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h 
 );

//Calculation Information
void Cal_Info
(
 int step, int step_o, DATA_TYPE step_t,
 int nol, int nos,int noe, 
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );



////////////////////////
// Main Calculation
////////////////////////
void Set_Dependent_Value
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound Bc
 );

void h_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );

void u_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void v_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void w_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void s_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void k_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );

void e_Advection
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );


void u_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void v_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void w_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void s_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );
void k_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );

void e_Diffusion
(
 int nol, int nos, int noe,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h
 );


void Pressure 
(
 int nol, int nos,int noe, 
 struct vert_geom *veG_h,
 struct surf_conn *suC_h, struct surf_geom *suG_h, struct surf *suR_h,
 struct cell_conn *ceC_h, struct cell_geom *ceG_h, struct cell *ceR_h,
 struct bound Bc
 );

void Sync
(
 int nol, int non, int noe,
 struct vert_geom *veG_h,
 struct cell *ceR_h
 );

DATA_TYPE Density_Value
(
 DATA_TYPE S,
 DATA_TYPE T,
 DATA_TYPE P
 );
