#include <iostream>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>

using namespace std;

void exercise4(void) {
  //==============================================================
  // Unfolding exercise 4
  //  plot purtity, fraction of events from neightbour bin,
  //  fraction of signal events
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
  TH1D *hist_mc1_rec=new TH1D("hist_mc1_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_recsig=new TH1D("hist_mc1_recsig",";r_{rec}",nBin,r0,r1);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // MC2
  TH1D *hist_mc2_rec=new TH1D("hist_mc2_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc2_recsig=new TH1D("hist_mc2_recsig",";r_{rec}",nBin,r0,r1);
  TH2D *hist_mc2_recgen=new TH2D("hist_mc2_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // derived histograms
  // purity
  TH1D *hist_mc1_purity=new TH1D("hist_mc1_purity",";r_{rec}",nBin,r0,r1);
  // signal fraction
  TH1D *hist_mc1_signal=new TH1D("hist_mc1_signal",";r_{rec}",nBin,r0,r1);
  // migration from neightbour bins
  TH1D *hist_mc1_near=new TH1D("hist_mc1_near",";r_{rec}",nBin,r0,r1);
  // purity
  TH1D *hist_mc2_purity=new TH1D("hist_mc2_purity",";r_{rec}",nBin,r0,r1);
  // signal fraction
  TH1D *hist_mc2_signal=new TH1D("hist_mc2_signal",";r_{rec}",nBin,r0,r1);
  // migration from neightbour bins
  TH1D *hist_mc2_near=new TH1D("hist_mc2_near",";r_{rec}",nBin,r0,r1);

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
    TString recWeight="w*isTrig";
    TString recsigWeight="w*isTrig*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recsig","rRec",recsigWeight,"",nEvent,firstEvent);
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
    TString recWeight="w*isTrig";
    TString recsigWeight="w*isTrig*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_recsig","rRec",recsigWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_recgen","rRec:rGen",recsigWeight,
		  "",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // fill derived histograms
  for(int iRec=1;iRec<=nBin;iRec++) {
    // all reconstructed events (includes background)
    double nRec=hist_mc1_rec->GetBinContent(iRec);
    if(nRec>0.0) {
      // fraction of all signal events
      double nSignal=hist_mc1_recsig->GetBinContent(iRec);
      hist_mc1_signal->SetBinContent(iRec,nSignal/nRec);

      // migration from neighbour bins
      double nNear=0.;
      for(int iGen=iRec-1;iGen<=iRec+1;iGen++) {
	double c=hist_mc1_recgen->GetBinContent(iGen,iRec);
	nNear +=c;
      }
      hist_mc1_near->SetBinContent(iRec,nNear/nRec);

      // purity
      double nGenRec=hist_mc1_recgen->GetBinContent(iRec,iRec);
      hist_mc1_purity->SetBinContent(iRec,nGenRec/nRec);
    }
  }
  for(int iRec=1;iRec<=nBin;iRec++) {
    // all reconstructed events (includes background)
    double nRec=hist_mc2_rec->GetBinContent(iRec);
    if(nRec>0.0) {
      // fraction of all signal events
      double nSignal=hist_mc2_recsig->GetBinContent(iRec);
      hist_mc2_signal->SetBinContent(iRec,nSignal/nRec);

      // migration from neighbour bins
      double nNear=0.;
      for(int iGen=iRec-1;iGen<=iRec+1;iGen++) {
	double c=hist_mc2_recgen->GetBinContent(iGen,iRec);
	nNear += c;
      }
      hist_mc2_near->SetBinContent(iRec,nNear/nRec);

      // purity
      double nGenRec=hist_mc2_recgen->GetBinContent(iRec,iRec);
      hist_mc2_purity->SetBinContent(iRec,nGenRec/nRec);
    }
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("purity and migrations","",300,600);
  canvas->Divide(1,2);

  canvas->cd(1);
  hist_mc1_signal->GetYaxis()->SetRangeUser(0.,1.5);
  hist_mc1_signal->SetLineColor(kCyan);
  hist_mc1_signal->SetFillColor(kCyan-10);
  hist_mc1_signal->SetFillStyle(1001);
  hist_mc1_signal->Draw("HIST");
  hist_mc1_near->SetLineColor(kBlue);
  hist_mc1_near->SetFillColor(kBlue-10);
  hist_mc1_near->SetFillStyle(1001);
  hist_mc1_near->Draw("HIST SAME");
  hist_mc1_purity->SetLineColor(kMagenta);
  hist_mc1_purity->SetFillColor(kMagenta-10);
  hist_mc1_purity->SetFillStyle(1001);
  hist_mc1_purity->Draw("HIST SAME");
  TLine *line=new TLine(r0,1.0,r1,1.);
  line->SetLineStyle(2);
  line->Draw();
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_mc1_signal,"MC1 signal fraction","fl");
  legend1->AddEntry(hist_mc1_near,"MC1 rec=gen#pm1","fl");
  legend1->AddEntry(hist_mc1_purity,"MC1 rec=gen [purity]","fl");
  legend1->Draw();

  canvas->cd(2);
  hist_mc2_signal->GetYaxis()->SetRangeUser(0.,1.5);
  hist_mc2_signal->SetLineColor(kCyan);
  hist_mc2_signal->SetFillColor(kCyan-10);
  hist_mc2_signal->SetFillStyle(1001);
  hist_mc2_signal->Draw("HIST");
  hist_mc2_near->SetLineColor(kBlue);
  hist_mc2_near->SetFillColor(kBlue-10);
  hist_mc2_near->SetFillStyle(1001);
  hist_mc2_near->Draw("HIST SAME");
  hist_mc2_purity->SetLineColor(kMagenta);
  hist_mc2_purity->SetFillColor(kMagenta-10);
  hist_mc2_purity->SetFillStyle(1001);
  hist_mc2_purity->Draw("HIST SAME");
  line->Draw();
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_mc2_signal,"MC2 signal fraction","fl");
  legend2->AddEntry(hist_mc2_near,"MC2 rec=gen#pm1","fl");
  legend2->AddEntry(hist_mc2_purity,"MC2 rec=gen [purity]","fl");
  legend2->Draw();

}
