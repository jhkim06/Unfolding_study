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

using namespace std;

void exercise5(void) {
  //==============================================================
  // Unfolding exercise 5
  //  compare matrix of probabilities for MC1 and MC2
  //  also, compare the efficiencies
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
  int nBin=9;
  double r0=0.;
  double r1=0.9;
  // MC1
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBin,r0,r1);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // MC2
  TH1D *hist_mc2_gen=new TH1D("hist_mc2_gen",";r_{gen}",nBin,r0,r1);
  TH2D *hist_mc2_recgen=new TH2D("hist_mc2_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // derived histograms
  // matrix of probabilities
  TH2D *hist_mc1_prob=new TH2D("hist_mc1_prob","MC1 response matrix;r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  TH2D *hist_mc2_prob=new TH2D("hist_mc2_prob","MC2 response matrix;r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc1_eff=new TH1D("hist_mc1_eff","Efficiency;r_{gen}",nBin,r0,r1);
  TH1D *hist_mc2_eff=new TH1D("hist_mc2_eff","Efficiency;r_{gen}",nBin,r0,r1);
  TH1D *hist_diff_prob=new TH1D("hist_diff_prob","MC1-MC2 differences;%Delta A_{ij}",20,-0.1,0.1);

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

    // fill histograms
    histogramDir->cd();

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
    TString recsigWeight="w*isTrig*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recgen","rRec:rGen",recsigWeight,
		  "",nEvent,firstEvent);

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
    TString recsigWeight="w*isTrig*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_recgen","rRec:rGen",recsigWeight,
		  "",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // fill derived histograms
  for(int iGen=1;iGen<=nBin;iGen++) {
    double nGen=hist_mc1_gen->GetBinContent(iGen);
    if(nGen>0.) {
      double nRec=0.0;
      for(int iRec=1;iRec<=nBin;iRec++) {
	double c=hist_mc1_recgen->GetBinContent(iGen,iRec);
	nRec += c;
	hist_mc1_prob->SetBinContent(iGen,iRec,c/nGen);
      }
      hist_mc1_eff->SetBinContent(iGen,nRec/nGen);
    }
  }
  for(int iGen=1;iGen<=nBin;iGen++) {
    double nGen=hist_mc2_gen->GetBinContent(iGen);
    if(nGen>0.) {
      double nRec=0.0;
      for(int iRec=1;iRec<=nBin;iRec++) {
	double c=hist_mc2_recgen->GetBinContent(iGen,iRec);
	nRec += c;
	hist_mc2_prob->SetBinContent(iGen,iRec,c/nGen);
      }
      hist_mc2_eff->SetBinContent(iGen,nRec/nGen);
    }
  }
  // make histogram of differences
  for(int iGen=1;iGen<=nBin;iGen++) {
    for(int iRec=1;iRec<=nBin;iRec++) {
      double dA=
	hist_mc1_prob->GetBinContent(iGen,iRec)-
	hist_mc2_prob->GetBinContent(iGen,iRec);
      hist_diff_prob->Fill(dA);
    }
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("probabilities and efficiencies","",600,600);
  canvas->Divide(2,2);

  TVirtualPad *p1=canvas->cd(1);
  p1->SetRightMargin(0.15);
  hist_mc1_prob->SetMinimum(-0.001);
  hist_mc1_prob->Draw("COLZ");
  hist_mc1_prob->Draw("BOX SAME");

  TVirtualPad *p2=canvas->cd(2);
  p2->SetRightMargin(0.15);
  hist_mc2_prob->SetMinimum(-0.001);
  hist_mc2_prob->Draw("COLZ");
  hist_mc2_prob->Draw("BOX SAME");

  canvas->cd(3);
  hist_diff_prob->GetYaxis()->SetRangeUser(0.,55.);
  hist_diff_prob->SetMarkerColor(kCyan);
  hist_diff_prob->SetLineColor(kCyan);
  hist_diff_prob->SetFillColor(kCyan-10);
  hist_diff_prob->SetFillStyle(1001);
  hist_diff_prob->Draw("HIST");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9);
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_diff_prob,"A_{ij}^{MC1}-A_{ij}^{MC2}","lf");
  legend3->Draw();

  canvas->cd(4);
  hist_mc1_eff->GetYaxis()->SetRangeUser(0.,1.5);
  hist_mc1_eff->SetLineColor(kBlue);
  hist_mc1_eff->Draw("HIST");
  hist_mc2_eff->SetLineColor(kRed);
  hist_mc2_eff->Draw("HIST SAME");
  TLegend *legend4=new TLegend(0.25,0.7,0.85,0.9);
  legend4->SetBorderSize(0);
  legend4->SetFillStyle(0);
  legend4->AddEntry(hist_mc1_eff,"efficiency MC1","l");
  legend4->AddEntry(hist_mc2_eff,"efficiency MC2","l");
  legend4->Draw();
}
