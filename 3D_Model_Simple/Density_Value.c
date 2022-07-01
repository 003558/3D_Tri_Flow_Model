// -*- c -*-
#include "public.h"

DATA_TYPE Density_Value
(
 DATA_TYPE S,
 DATA_TYPE T,
 DATA_TYPE P
 )
{
  static DATA_TYPE a[6];
  static DATA_TYPE b[5];
  static DATA_TYPE c[3];
  static DATA_TYPE d[1];
  static DATA_TYPE e[5];
  static DATA_TYPE f[4];
  static DATA_TYPE g[3];
  static DATA_TYPE h[4];
  static DATA_TYPE i[3];
  static DATA_TYPE j[1];
  static DATA_TYPE k[3];
  static DATA_TYPE m[3];

  static DATA_TYPE t,p;
  static DATA_TYPE R0,Rw;
  static DATA_TYPE Kp,K0;
  static DATA_TYPE A,B;
  static DATA_TYPE Kw,Aw,Bw;

  static DATA_TYPE Rho;

  a[0]=((DATA_TYPE)999.842594);
  a[1]=((DATA_TYPE)6.793952*1.0e-2);
  a[2]=((DATA_TYPE)-9.095290*1.0e-3);
  a[3]=((DATA_TYPE)1.001685*1.0e-4);
  a[4]=((DATA_TYPE)-1.120083*1.0e-6);
  a[5]=((DATA_TYPE)6.536332*1.0e-9);

  b[0]=((DATA_TYPE)8.24493*1.0e-1);
  b[1]=((DATA_TYPE)-4.0899*1.0e-3);
  b[2]=((DATA_TYPE)7.6438*1.0e-5);
  b[3]=((DATA_TYPE)-8.2467*1.0e-7);
  b[4]=((DATA_TYPE)5.3875*1.0e-9);

  c[0]=((DATA_TYPE)-5.72466*1.0e-3);
  c[1]=((DATA_TYPE)1.0227*1.0e-4);
  c[2]=((DATA_TYPE)-1.6546*1.0e-6);

  d[0]=((DATA_TYPE)4.8314*1.0e-4);

  e[0]=((DATA_TYPE)19652.21);
  e[1]=((DATA_TYPE)148.4206);
  e[2]=((DATA_TYPE)-2.327105);
  e[3]=((DATA_TYPE)1.360477*1.0e-2);
  e[4]=((DATA_TYPE)-5.155288*1.0e-5);

  f[0]=((DATA_TYPE)54.6746);
  f[1]=((DATA_TYPE)-0.603459);
  f[2]=((DATA_TYPE)1.09987*1.0e-2);
  f[3]=((DATA_TYPE)-6.1670*1.0e-5);

  g[0]=((DATA_TYPE)7.944*1.0e-2);
  g[1]=((DATA_TYPE)1.6483*1.0e-2);
  g[2]=((DATA_TYPE)-5.3009*1.0e-4);

  h[0]=((DATA_TYPE)3.239908);
  h[1]=((DATA_TYPE)1.43713*1.0e-3);
  h[2]=((DATA_TYPE)1.16092*1.0e-4);
  h[3]=((DATA_TYPE)-5.77905*1.0e-7);

  i[0]=((DATA_TYPE)2.2838*1.0e-3);
  i[1]=((DATA_TYPE)-1.0981*1.0e-5);
  i[2]=((DATA_TYPE)-1.6078*1.0e-6);

  j[0]=((DATA_TYPE)1.91075*1.0e-4);

  k[0]=((DATA_TYPE)8.50935*1.0e-5);
  k[1]=((DATA_TYPE)-6.12293*1.0e-6);
  k[2]=((DATA_TYPE)5.2787*1.0e-8);

  m[0]=((DATA_TYPE)-9.9348*1.0e-7);
  m[1]=((DATA_TYPE)2.0816*1.0e-8);
  m[2]=((DATA_TYPE)9.1697*1.0e-10);


  t=T*((DATA_TYPE)1.00024);
  p=P*((DATA_TYPE)1.0e-5);

  Bw=k[0]+k[1]*t+k[2]*t*t;
  Aw=h[0]+h[1]*t+h[2]*t*t+h[3]*t*t*t;
  Kw=e[0]+e[1]*t+e[2]*t*t+e[3]*t*t*t+e[4]*t*t*t*t;

  B=Bw+(m[0]+m[1]*t+m[2]*t*t)*S;
  A=Aw+(i[0]+i[1]*t+i[2]*t*t)*S+j[0]*POW(S,((DATA_TYPE)1.5));

  K0=Kw+(f[0]+f[1]*t+f[2]*t*t+f[3]*t*t*t)*S+(g[0]+g[1]*t+g[2]*t*t)*POW(S,((DATA_TYPE)1.5));
  Kp=K0+A*p+B*p*p;

  Rw=a[0]+a[1]*t+a[2]*t*t+a[3]*t*t*t+a[4]*t*t*t*t+a[5]*t*t*t*t*t;
  R0=Rw+(b[0]+b[1]*t+b[2]*t*t+b[3]*t*t*t+b[4]*t*t*t*t)*S+(c[0]+c[1]*t+c[2]*t*t)*POW(S,((DATA_TYPE)1.5))+d[0]*S*S;

  Rho=R0/(((DATA_TYPE)1.0)-p/Kp);
  
  return Rho;    
}
