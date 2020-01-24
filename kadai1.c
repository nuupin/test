
#include <stdio.h>
#include <math.h>
#define M 600

double x[M];
double g[M];
double mass[M];
int nlist;
int list[M][2];
double energy1(int natom, double x[], int nlist, int list[][2], double g[]);
void getfx(double x[], double fx[],int mm);

void initialize(double xo[], int mm)
{
   int i,nn;
   nn = mm/2;
   for (i=0; i<nn; i++) xo[i] = x[i];
   for (i=0; i<nn; i++) xo[nn+i] = 0;
   for (i=0; i<nn; i++) mass[i] = 1;
}

int input_mol(double x[], int list[][2], int *nlist)
{
  int natom,i;
  scanf("%d",&natom);
  scanf("%d",nlist);
  for(i=0;i<*nlist;i++)
    scanf("%d%d",&list[i][0],&list[i][1]);
  for(i=0;i<natom;i++) {
    scanf("%lf%lf%lf",&x[3*i],&x[3*i+1],&x[3*i+2]);
  }
  return natom;
}

void prt_mol1(double x[], int natom)
{
  int i;
  for(i=0; i<natom; i++) {
    printf("%10.5f %10.5f %10.5f\n",x[3*i],x[3*i+1],x[3*i+2]);
  }
}

int main(void)
{
  int i,n,k,mm;
  int natom;
  double xo[M],xn[M],dt;
  double fx0[M],fx1[M],fx2[M],fx3[M],x1[M],x2[M],x3[M];
  natom = input_mol(x,list,&nlist);
  mm = 6*natom;
  n = 10000;
  dt = 0.01;
  initialize(xo,mm);

  for(i=0;i<n;i++){
     prt_mol1(xo,natom);

     getfx(xo,fx0,mm);
     for(k=0;k<mm;k++)  x1[k] = xo[k] + fx0[k] * dt/2;
     getfx(x1,fx1,mm);
     for(k=0;k<mm;k++)  x2[k] = xo[k] + fx1[k] * dt/2;
     getfx(x2,fx2,mm);
     for(k=0;k<mm;k++)  x3[k] = xo[k] + fx2[k] * dt;
     getfx(x3,fx3,mm);
     for(k=0;k<mm;k++)
        xn[k] = xo[k] + (fx0[k] + 2*fx1[k] + 2*fx2[k] + fx3[k])*dt/6;

     for(k=0;k<mm;k++) xo[k] = xn[k];
  }
  return 0;
}

void getfx(double x[], double fx[],int mm)
{
   int n,i;
   double b = 0.001;
   n = mm/2;
   for (i=0; i<n; i++)
      fx[i] = x[n+i];
   for (i=0; i<n; i++) g[i] = 0;
   energy1(n/3,x,nlist,list, g);
   for (i=0; i<n; i++) {
      fx[n+i] = -(g[i]+b*x[n+i])/mass[i];
   }
}

double distance(double x[], int i, int j)
{
  return sqrt((x[3*i  ]-x[3*j  ])*(x[3*i  ]-x[3*j  ])
	   +  (x[3*i+1]-x[3*j+1])*(x[3*i+1]-x[3*j+1])
	   +  (x[3*i+2]-x[3*j+2])*(x[3*i+2]-x[3*j+2]) );
}

double univ(int i, int j, double x[], double eij[3])
{
  double rij;
  rij = distance(x,i,j);
  if (rij < 1e-9) return 0;
  eij[0] = (x[3*i  ] - x[3*j  ])/rij;
  eij[1] = (x[3*i+1] - x[3*j+1])/rij;
  eij[2] = (x[3*i+2] - x[3*j+2])/rij;
  return rij;
}

double energy1(int natom, double x[], int nlist, int list[][2], double g[])
{
  int i,j,m;
  double k,re,rij,e,eij[3],e1;
  k = 1.0;
  re = 1.0;
  e = 0;
  for (m=0; m<nlist; m++) {
    i=list[m][0];
    j=list[m][1];
    rij = distance(x,i,j);
    e += k*(rij-re)*(rij-re);
    e1 = 2* k*(rij-re);
    univ(i,j,x,eij);
    g[3*i]   +=  e1*eij[0];
    g[3*i+1] +=  e1*eij[1];
    g[3*i+2] +=  e1*eij[2];
    g[3*j]   += -e1*eij[0];
    g[3*j+1] += -e1*eij[1];
    g[3*j+2] += -e1*eij[2];

  }
  return e;
}
