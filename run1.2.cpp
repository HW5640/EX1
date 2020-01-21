///
/// @file 
/// @brief Test of Runge-Kutta solver for series of ODEs
/// @author Bob Hirosky
/// @date 31 Dec 2019 
/// 
/// Use the Rk4 solver for coupled ODEs to solve for projectile 
/// motion with air resistance
///
/// Definition of our variables
/// x    = time <br>
/// y[0] = position along i axis  ; f_ri = dri/dt => velocity along i axis  <br>
/// y[1] = velocity along i axis  ; f_vi = dvi/dt => acceleration along i axis <br>
/// y[2] = position along j axis  ; f_rj = drj/dt => velocity along j axis <br>
/// y[3] = velocity along j axis  ; f_vj = dvj/dt => acceleration along j axis <br>



#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;


int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************
  ///////////////// Problem 1.2a  ////////////////////
  std::vector<double> x,y;
  double log_rg[10],logN[10];   

  TH2F *h2 = new TH2F("DLA","",150,-5,5,150,-5,5);
  h2->SetBinContent(75,75,1);
  
  TRandom2 tr(1);
  int npoints = 5000;
  int i = 1;
  
  while(i<npoints){
    int bound = tr.Integer(596) + 1;
    int bound_x, bound_y;
    if(bound <= 150) {
      bound_x = bound;
      bound_y = 0;
    }
    if(bound > 150 && bound < 300) {
      bound_x = 150;
      bound_y = bound - 150;
    }
    if(bound >= 300 && bound <= 450) {
      bound_x = bound - 300;
      bound_y = 150;
    }
    if(bound > 450 && bound < 600) {
      bound_x = 0;
      bound_y = bound - 450;
    }

    bool site = false;
    while(!site){
      int num = tr.Integer(4); //0 = left, 1 = right, 2 = down, 3 = up
      if(bound_x == 0 && num == 0) continue;
      if(bound_x == 150 && num == 1) continue;
      if(bound_y == 0 && num == 2) continue;
      if(bound_y == 150 && num == 3) continue;

      if(num == 0) bound_x--;
      if(num == 1) bound_x++;
      if(num == 2) bound_y--;
      if(num == 3) bound_y++;
      
      if(h2->GetBinContent(bound_x,bound_y,1) == 1){
	
	if(num == 0) bound_x++;
	if(num == 1) bound_x--;
	if(num == 2) bound_y++;
	if(num == 3) bound_y--;
	h2->SetBinContent(bound_x,bound_y,1);
	x.push_back(bound_x);
	y.push_back(bound_y);
	site = true;
	i++;
      }
    }
  }

  int palette[1] = {kBlack};
  gStyle->SetPalette(1,palette);
  gStyle->SetOptStat(0);

  h2->Draw("cola");

  //////////  Problem 1.2b  /////////////////////////
  FILE* const fout = fopen("DLA_Rg.dat","w");

  for(int j = 0; j<10;j++){
    double r_0 = 0;
    double sum = 0;
    int N = (j+1)*200;
    for(int i=0;i<N;i++){
      r_0 += 1./N*sqrt(x[i]*x[i] + y[i]*y[i]);
    }

    for(int i=0;i<N;i++){
      sum += x[i]*x[i] + y[i]*y[i];
    }

    double R_g = 1./N*sum - r_0*r_0;

    logN[j] = log(N);
    log_rg[j] = log(sqrt(R_g));

    fprintf(fout,"%lf %lf\n", log_rg[j], logN[j]);
  }
  fclose(fout);

  TCanvas *c2 = new TCanvas();
  TGraph *g = new TGraph(10,log_rg,logN);
  g->Draw("AP");
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.5);
  g->Fit("pol1");
  gStyle->SetOptFit();
  g->SetTitle(";log(R_{g}^{2});log(N)");

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

