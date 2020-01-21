/// \file 
/// \brief simple example of using a function from the gsl library

#include <cstdio>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_matrix.h>
#include <gsl_odeiv2.h>
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>


using namespace std;
const double b=0.25;
const double a=0.5;
const double F_0=2.0;
const double omega=2.4;
const double gam=0.1;
const double m=1.0; 
bool first = true;

int func(double x, const double y[], double f[], void *params){
  f[0] = y[1];
  f[1] = (-gam*y[1] + 2*a*y[0] - 4*b*pow(y[0],3) + F_0*cos(omega*x) ) / m;
  return GSL_SUCCESS;
}


int ode(double x0, double xf,int npoints, double y[], double x, gsl_odeiv2_driver * d, TGraph *g){
  for(int i = 1;i <= npoints;i++){
    double xi = x0 + i*(xf-x0)/(double)npoints;
    gsl_odeiv2_driver_apply(d,&x,xi,y);

    g->SetPoint(i,x,y[0]);
  }

  return 0;
}

int main (int argc, char *argv[])
{
  TApplication theApp("App", &argc, argv);

  int dim = 2;
  gsl_odeiv2_system sys = {func, NULL, dim, NULL};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
 
  double x0 = 0.0;
  double xf = 100.0;
  double x = x0;
  double y[2] = {0.5,0};
  int npoints = 1000;

  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  ode(x0,xf,npoints,y,x,d,g1);

  y[0] = 0.5001;
  ode(x0,xf,npoints,y,x,d,g2);

  g1->Draw("AL");
  g2->Draw("L");
  g2->SetLineColor(kRed);

  gsl_odeiv2_driver_free(d);

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q"); 
                                                                                
  theApp.Run();
  
  return(0);
}

