#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "RooDataHist.h"
#include "TH1.h"
#include "TH2.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooNLLVar.h"
#include "TFile.h"
#include "RooJohnson.h"
#include "RooHypatia2.h"

using namespace RooFit ;
using namespace std;

void fittingD0()
{
  TFile * f1 = new TFile("output.root");   
  TH1 *h1 = (TH1*)f1-> Get("MD0");		//Combinational Mass with no conditions
  
  RooRealVar mass("mass","", 0, 6);
  RooDataHist data("data","data", RooArgList(mass),h1);

  //signal
  RooRealVar mu_J1("#mu_{J}","mean_johnson", 2.288, 0, 6);                                                   
  RooRealVar sigma_J1 ("#sigma_{J}", "sigma_johnson", 0.03, 0.0000001, 1);                                                              
  RooRealVar gamma_J1("#gamma_{J}", "gamma", 	1.0, -10, 10);
  RooRealVar delta_J1("#delta_{J}", "delta", 12, 0.0000001, 20);
  RooJohnson sig("johnson1", "Johnson PDF", mass, mu_J1, sigma_J1, gamma_J1, delta_J1);

  //background
  RooRealVar a1("a1", "Slope1 of Polynomial", 0.4, -1e6, 1e6);
  RooChebychev bkg("bkg","Chebyshev Polynomial", mass,RooArgList(a1));

  // RooRealVar a1("a1", "Slope1 of Polynomial", 0);
  // RooChebychev bkg("bkg","Chebshev Polynomial", mass,RooArgList(a1));
  
  
    
  //total
  RooRealVar n_sig("N_{sig}", "n_{s}", 2e4, 0, 1e10);
  RooRealVar n_bkg("N_{bkg}", "n_{b}", 2e1, 0, 1e2);
  // RooRealVar n_bkg("N_{bkg}", "n_{b}", 0);

  RooAddPdf pdf("pdf", "two component model",RooArgList(sig, bkg), RooArgList(n_sig, n_bkg));

  //fit
  RooFitResult *fitresult = pdf.fitTo(data, Save(true), Strategy(2), Extended(true));
      
  RooPlot* massplot = mass.frame();
  data.plotOn(massplot);
  pdf.plotOn(massplot);
  pdf.paramOn(massplot);
  cout << "chi^2 = " << massplot->chiSquare() << endl;
  double nbins = 12*massplot->chiSquare(8) / (massplot->chiSquare(8) - massplot->chiSquare());// 10 is the no. of fit parameters
  cout << "chi^2/ndf = " << massplot->chiSquare()*nbins << "/" << nbins-8 << endl;
  cout << "chi^2/ndf = " << massplot->chiSquare(8)<< endl;
  cout << "500bin_histo: chi^2/ndf = " << massplot->chiSquare()*500 << "/" << 500-8 << endl;
  
  RooHist* hpull = massplot->pullHist();
  hpull->SetFillStyle(1001);
  hpull->SetFillColor(1);
  for(int i=0;i<hpull->GetN();++i) hpull->SetPointError(i,0.0,0.0,0.0,0.0);
  RooPlot* pullplot = mass.frame(Title("Pull Plot"));
  pullplot->addPlotable(hpull,"B");
  //pullplot->SetYTitle("Pull");
  pullplot->SetXTitle("#it{m(D0)} [GeV/c^{2}]");
  pullplot->SetMinimum(-150.);
  pullplot->SetMaximum(150.);
  pullplot->GetXaxis()->SetLabelSize(0.1);
  pullplot->GetXaxis()->SetTitleSize(0.1);
  pullplot->GetYaxis()->SetLabelSize(0.07);
  
  TCanvas *canvas = new TCanvas("canvas","canvas", 600, 600);

  Double_t xlow, ylow, xup, yup;
  canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);
  canvas->Divide(1,2);

  TVirtualPad *upPad = canvas->GetPad(1);
  upPad->SetPad(xlow,ylow+0.25*(yup-ylow),xup,yup);
  
  TVirtualPad *dwPad = canvas->GetPad(2);
  dwPad->SetPad(xlow,ylow,xup,ylow+0.25*(yup-ylow));

  canvas->Update();
  canvas->cd(1);
  
  pdf.plotOn(massplot,Components(RooArgSet(sig)),LineColor(kRed),LineStyle(kDashed));
  //pdf.plotOn(massplot,Components(johnson),LineStyle(kDashed),LineColor(kCyan));
  //pdf.plotOn(massplot,Components(gauss),LineStyle(kDashed),LineColor(kBlack));
  //pdf.plotOn(massplot,Components(gauss1),LineStyle(kDashed),LineColor(kCyan));
  //pdf.plotOn(massplot,Components(gauss2),LineStyle(kDashed),LineColor(kBlack));
  pdf.plotOn(massplot,Components(bkg),LineStyle(kDashed),LineColor(kGreen));


  massplot->Draw();
  canvas->cd(2);
  pullplot->Draw();
  
  canvas->Update();
   
}

  



  
