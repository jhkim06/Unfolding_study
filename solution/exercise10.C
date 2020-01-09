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

using namespace std;

void exercise10(void) {
  //==============================================================
  // Unfolding exercise 10
  //  unfold data using  MC1 with TUnfold (tau=0) and
  //    10,20,100 bins on detector level
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

#define NUM_UNFOLD 3

  // book histograms here
  int nBinGen=10;
  double r0Gen=0.;
  double r1Gen=1.0;
  int nBinRec[NUM_UNFOLD]={10,20,100};
  double r0Rec=0.;
  double r1Rec=1.0;
  // data
  TH1D *hist_data_rec[NUM_UNFOLD];
  for(int iUnfold=0;iUnfold<NUM_UNFOLD;iUnfold++) {
    hist_data_rec[iUnfold]=
      new TH1D(TString::Format("hist_data_rec%d",iUnfold),
	       ";r_{rec}",nBinRec[iUnfold],r0Rec,r1Rec);
  }
  TH1D *hist_data_truth=new TH1D
    ("hist_data_truth",";r_{gen}",nBinGen,r0Gen,r1Gen);
  // MC1
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBinGen,r0Gen,r1Gen);
  TH2D *hist_mc1_recgen[NUM_UNFOLD];
  for(int iUnfold=0;iUnfold<NUM_UNFOLD;iUnfold++) {
    hist_mc1_recgen[iUnfold]=
      new TH2D(TString::Format("hist_mc1_recgen%d",iUnfold),";r_{gen};r_{rec}",
	       nBinGen,r0Gen,r1Gen,nBinRec[iUnfold],r0Rec,r1Rec);
  }

  // derived histograms
  // unfolding results
  TH1D *hist_data_NOTAU[NUM_UNFOLD];
  TH2D *hist_data_NOTAUrho[NUM_UNFOLD];
  for(int iUnfold=0;iUnfold<NUM_UNFOLD;iUnfold++) {
    hist_data_NOTAU[iUnfold]=
      new TH1D(TString::Format("hist_data_NOTAU%d",iUnfold),";r_{NOTAU}",
	       nBinGen,r0Gen,r1Gen);
    hist_data_NOTAUrho[iUnfold]=
      new TH2D(TString::Format("hist_data_NOTAUrho%d",iUnfold),
	       ";r_{NOTAU};r_{NOTAU}",nBinGen,r0Gen,r1Gen,nBinGen,r0Gen,r1Gen);
  }

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
    for(int iUnfold=0;iUnfold<NUM_UNFOLD;iUnfold++) {
      tree->Project(TString::Format("hist_data_rec%d",iUnfold),
		    "rRec",recWeight,"",nEvent,firstEvent);
    }
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
    for(int iUnfold=0;iUnfold<NUM_UNFOLD;iUnfold++) {
      // event with generator and reco information
      tree->Project(TString::Format("hist_mc1_recgen%d",iUnfold),
		    "rRec:rGen",recsigWeight,"",nEvent,firstEvent);
      // background event with reco information
      tree->Project(TString::Format("+hist_mc1_recgen%d",iUnfold),
		    "rRec:0.95",recbgrWeight,"",nEvent,firstEvent);
      // events with generator but no reco information
      tree->Project(TString::Format("+hist_mc1_recgen%d",iUnfold),
		    "-1.0:rGen",norecsigWeight,"",nEvent,firstEvent);
      // background event without reco information
      tree->Project(TString::Format("+hist_mc1_recgen%d",iUnfold),
		    "-1.0:0.95",norecbgrWeight,"",nEvent,firstEvent);
    }
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
  for(int iUnfold=0;iUnfold<3;iUnfold++) {
    TUnfold unfold(hist_mc1_recgen[iUnfold],
		   TUnfold::kHistMapOutputHoriz,
		   TUnfold::kRegModeSize,
		   TUnfold::kEConstraintNone);
    unfold.SetInput(hist_data_rec[iUnfold]);
    unfold.DoUnfold(0.0);
    unfold.GetOutput(hist_data_NOTAU[iUnfold]);
    unfold.GetRhoIJ(hist_data_NOTAUrho[iUnfold]);
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("TUnfold tau-0 10,15,20 bins","",900,600);
  canvas->Divide(3,2);

  double ymax=3500.;

  canvas->cd(1);
  hist_data_NOTAU[0]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_NOTAU[0]->SetLineColor(kRed);
  hist_data_NOTAU[0]->SetMarkerColor(kRed);
  hist_data_NOTAU[0]->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9,"TUnfold, #tau=0");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_NOTAU[0],
		    TString::Format("data nRec=%d",
				    hist_data_rec[0]->GetNbinsX()),"l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  TVirtualPad *p4=canvas->cd(4);
  p4->SetRightMargin(0.15);
  hist_data_NOTAUrho[0]->Draw("COLZ");
  hist_data_NOTAUrho[0]->Draw("BOX SAME");

  canvas->cd(2);
  hist_data_NOTAU[1]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_NOTAU[1]->SetLineColor(kRed);
  hist_data_NOTAU[1]->SetMarkerColor(kRed);
  hist_data_NOTAU[1]->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9,"TUnfold, #tau=0");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_data_NOTAU[1],
		    TString::Format("data nRec=%d",
				    hist_data_rec[1]->GetNbinsX()),"l");
  legend2->AddEntry(hist_data_truth,"data truth","l");
  legend2->AddEntry(mc1_truth,"MC1 truth","pl");
  legend2->Draw();

  TVirtualPad *p5=canvas->cd(5);
  p5->SetRightMargin(0.15);
  hist_data_NOTAUrho[1]->Draw("COLZ");
  hist_data_NOTAUrho[1]->Draw("BOX SAME");

  canvas->cd(3);
  hist_data_NOTAU[2]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_NOTAU[2]->SetLineColor(kRed);
  hist_data_NOTAU[2]->SetMarkerColor(kRed);
  hist_data_NOTAU[2]->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9,"TUnfold, #tau=0");
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_data_NOTAU[2],
		    TString::Format("data nRec=%d",
				    hist_data_rec[2]->GetNbinsX()),"l");
  legend3->AddEntry(hist_data_truth,"data truth","l");
  legend3->AddEntry(mc1_truth,"MC1 truth","pl");
  legend3->Draw();

  TVirtualPad *p6=canvas->cd(6);
  p6->SetRightMargin(0.15);
  hist_data_NOTAUrho[2]->Draw("COLZ");
  hist_data_NOTAUrho[2]->Draw("BOX SAME");
}
