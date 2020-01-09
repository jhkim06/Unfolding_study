#include <iostream>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

void exercise0(void) {
  //==============================================================
  // Unfolding exercise 0 
  //  plot rRec in data and MC
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
  TH1D *hist_data_total=new TH1D("hist_data_total",";r_{rec}",50,0.,1.);
  TH1D *hist_mc1_total=new TH1D("hist_mc1_total",";r_{rec}",50,0.,1.);
  TH1D *hist_mc2_total=new TH1D("hist_mc2_total",";r_{rec}",50,0.,1.);

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
    tree->Project("hist_data_total","rRec",recWeight,"",nEvent,firstEvent);

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
    TString recWeight="w*isTrig";
    TString bgrWeight="w*isTrig*isBgr";
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
    TString bgrWeight="w*isTrig*isBgr";
    TString signalWeight="w*isTrig*(1-isBgr)";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_total","rRec",recWeight,"",nEvent,firstEvent);

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
    TString bgrWeight="w*isTrig*isBgr";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_total","rRec",recWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("DATA-MC comparison (variable rRec)","rRec",300,600);
  canvas->Divide(1,2);
  canvas->cd(1);
  hist_data_total->SetLineColor(kRed);
  hist_data_total->SetMarkerColor(kRed);
  hist_data_total->Draw();
  hist_mc1_total->SetLineColor(kBlue);
  hist_mc1_total->Draw("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.7,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_total,"data","pl");
  legend1->AddEntry(hist_mc1_total,"MC1","l");
  legend1->Draw();

  canvas->cd(2);
  hist_data_total->SetLineColor(kRed);
  hist_data_total->SetMarkerColor(kRed);
  hist_data_total->Draw();
  hist_mc2_total->SetLineColor(kBlue);
  hist_mc2_total->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.7,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_data_total,"data","pl");
  legend2->AddEntry(hist_mc2_total,"MC2","l");
  legend2->Draw();

}
