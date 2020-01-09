#include <iostream>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TUnfold.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TMarker.h>

using namespace std;

void exercise11(void) {
  //==============================================================
  // Unfolding exercise 11
  //  unfold data using  MC1 with TUnfold, L-curve method
  //    40 bins on detector level
  //==============================================================

  // the event samples are divided into subsets
  //  for most exercises we divide the event samples into 30 subsets
  //  and use only the first subset
  int nSubset=30;
  int iSubset=0;

  // calculate proper uncertainties when filling weighted events
  //   into histograms
  TH1::SetDefaultSumw2();

  // change to memory
  gDirectory->cd("Rint:/");
  TDirectory *histogramDir=gDirectory;

  // book histograms here
  int nBinGen=10;
  double r0Gen=0.;
  double r1Gen=1.0;
  int nBinRec=40;
  double r0Rec=0.;
  double r1Rec=1.0;
  // data
  TH1D *hist_data_rec=new TH1D("hist_data_rec",";r_{rec}",nBinRec,r0Rec,r1Rec);
  TH1D *hist_data_truth=new TH1D
    ("hist_data_truth",";r_{gen}",nBinGen,r0Gen,r1Gen);
  // MC1
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBinGen,r0Gen,r1Gen);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBinGen,r0Gen,r1Gen,nBinRec,r0Rec,r1Rec);
  // derived histograms
  // unfolding results
  TH1D *hist_data_LCURVE=new TH1D("hist_data_LCURVE",";r_{LCURVE}",
				  nBinGen,r0Gen,r1Gen);
  TH2D *hist_data_LCURVErho=
    new TH2D("hist_data_LCURVErho",";r_{LCURVE};r_{LCURVE}",
	     nBinGen,r0Gen,r1Gen,nBinGen,r0Gen,r1Gen);

  // loop over data events
  {
    // open input file and  get access to the events
    TFile inputFile("data.root");
    TTree *tree;
    inputFile.GetObject("rec",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString recWeight="w*isTrig";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_data_rec","rRec",recWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // loop over data truth events
  {
    // open input file and  get access to the events
    TFile inputFile("datatruth.root");
    TTree *tree;
    inputFile.GetObject("gen",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString genWeight="w*(1-isBgr)";
    TString bgrWeight="w*isBgr";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_data_truth","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("+hist_data_truth","0.95",bgrWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // loop over mc1 events
  {
    // open input file and  get access to the events
    TFile inputFile("mc1.root");
    TTree *tree;
    inputFile.GetObject("recgen",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString genWeight="w*(1-isBgr)";
    TString bgrWeight="w*isBgr";
    TString recWeight="w*isTrig";
    TString recbgrWeight="w*isTrig*isBgr";
    TString recsigWeight="w*isTrig*(1-isBgr)";
    TString norecbgrWeight="w*(1-isTrig)*isBgr";
    TString norecsigWeight="w*(1-isTrig)*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("+hist_mc1_gen","0.95",bgrWeight,"",nEvent,firstEvent);
    // event with generator and reco information
    tree->Project("hist_mc1_recgen",
		  "rRec:rGen",recsigWeight,"",nEvent,firstEvent);
    // background event with reco information
    tree->Project("+hist_mc1_recgen",
		  "rRec:0.95",recbgrWeight,"",nEvent,firstEvent);
    // events with generator but no reco information
    tree->Project("+hist_mc1_recgen",
		  "-1.0:rGen",norecsigWeight,"",nEvent,firstEvent);
    // background event without reco information
    tree->Project("+hist_mc1_recgen",
		  "-1.0:0.95",norecbgrWeight,"",nEvent,firstEvent);
    // clean up
    delete tree;
  }

  // loop over mc2 events
  {
    // open input file and  get access to the events
    TFile inputFile("mc2.root");
    TTree *tree;
    inputFile.GetObject("recgen",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString genWeight="w*(1-isBgr)";
    TString bgrWeight="w*isBgr";
    TString recWeight="w*isTrig";

    // fill histograms

    // clean up
    delete tree;
  }

  // run the unfolding and retreive results
  TUnfold unfold(hist_mc1_recgen,
		 TUnfold::kHistMapOutputHoriz,
		 TUnfold::kRegModeSize,
		 TUnfold::kEConstraintNone);
  unfold.SetInput(hist_data_rec,1.0);
  TGraph *lcurve=0;
  TSpline *logTauX=0,*logTauY=0;
  unfold.ScanLcurve(100,0.,0.,&lcurve,&logTauX,&logTauY);
  unfold.GetOutput(hist_data_LCURVE);
  unfold.GetRhoIJ(hist_data_LCURVErho);
  double tau=unfold.GetTau();
  double logTau=TMath::Log10(tau);
  double lcurveX=logTauX->Eval(logTau);
  double lcurveY=logTauY->Eval(logTau);
  cout<<logTau<<" "<<lcurveX<<" "<<lcurveY<<"\n";

  // make plot
  TCanvas *canvas=
    new TCanvas("TUnfold L-curve scan","",600,600);
  canvas->Divide(2,2);

  double ymax=3500.;

  canvas->cd(1);
  hist_data_LCURVE->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_LCURVE->SetLineColor(kRed);
  hist_data_LCURVE->SetMarkerColor(kRed);
  hist_data_LCURVE->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9,"L-curve scan");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_LCURVE,"data unfolded","l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  TVirtualPad *p3=canvas->cd(3);
  p3->SetRightMargin(0.15);
  hist_data_LCURVErho->Draw("COLZ");
  hist_data_LCURVErho->Draw("BOX SAME");

  canvas->cd(2);
  lcurve->SetLineColor(kCyan+2);
  lcurve->SetLineWidth(2);
  lcurve->SetTitle("L curve;X;Y");
  lcurve->Draw("AL");
  TMarker *bestPoint=new TMarker(lcurveX,lcurveY,20);
  bestPoint->SetMarkerSize(1.3);
  bestPoint->SetMarkerColor(kRed);
  bestPoint->Draw();

  canvas->cd(4);
  TVirtualPad *p4=gPad;
  p4->Divide(1,2);
  p4->cd(1);
  logTauX->SetTitle(";log_{10}(#tau);X");
  logTauX->SetLineColor(kCyan);
  logTauX->SetLineWidth(2);
  logTauX->Draw();
  TMarker *bestPointX=new TMarker(logTau,lcurveX,20);
  bestPointX->SetMarkerSize(1.3);
  bestPointX->SetMarkerColor(kRed);
  bestPointX->Draw();
  p4->cd(2);
  logTauY->SetTitle(";log_{10}(#tau);Y");
  logTauY->SetLineColor(kCyan);
  logTauY->SetLineWidth(2);
  logTauY->Draw();
  TMarker *bestPointY=new TMarker(logTau,lcurveY,20);
  bestPointY->SetMarkerSize(1.3);
  bestPointY->SetMarkerColor(kRed);
  bestPointY->Draw();
}
