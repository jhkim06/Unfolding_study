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

void exercise1(void) {
  //==============================================================
  // Unfolding exercise 1
  //  plot rRec in data and MC
  //  superimpose background
  //  add panels for generated signal
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
  TH1D *hist_data_truth=new TH1D("hist_data_truth",";r_{gen}",50,0.,1.);
  TH1D *hist_mc1_total=new TH1D("hist_mc1_total",";r_{rec}",50,0.,1.);
  TH1D *hist_mc1_bgr=new TH1D("hist_mc1_bgr",";r_{rec}",50,0.,1.);
  TH1D *hist_mc1_signalGen=new TH1D("hist_mc1_signalGen",";r_{gen}",50,0.,1.);
  TH1D *hist_mc2_total=new TH1D("hist_mc2_total",";r_{rec}",50,0.,1.);
  TH1D *hist_mc2_bgr=new TH1D("hist_mc2_bgr",";r_{rec}",50,0.,1.);
  TH1D *hist_mc2_signalGen=new TH1D("hist_mc2_signalGen",";r_{gen}",50,0.,1.);

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
    tree->Project("hist_data_truth","rGen",genWeight,"",nEvent,firstEvent);

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
    tree->Project("hist_mc1_bgr","rRec",bgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_signalGen","rGen",genWeight,"",nEvent,firstEvent);

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
    tree->Project("hist_mc2_bgr","rRec",bgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_signalGen","rGen",genWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("DATA-MC comparisons","rRec and rGen",600,600);
  canvas->Divide(2,2);
  canvas->cd(1);
  hist_data_total->SetLineColor(kRed);
  hist_data_total->SetMarkerColor(kRed);
  hist_data_total->Draw();
  hist_mc1_total->SetLineColor(kBlue);
  hist_mc1_total->Draw("HIST SAME");
  hist_mc1_bgr->SetLineColor(kBlack);
  hist_mc1_bgr->SetFillColor(kGray);
  hist_mc1_bgr->SetFillStyle(1001);
  hist_mc1_bgr->Draw("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.7,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_total,"data","pl");
  legend1->AddEntry(hist_mc1_total,"MC1","l");
  legend1->Draw();

  canvas->cd(2);
  hist_mc1_signalGen->SetLineColor(kBlue);
  hist_mc1_signalGen->Draw("HIST");
  hist_data_truth->Draw("HIST SAME");
  TLegend *legend1gen=new TLegend(0.25,0.7,0.7,0.9);
  legend1gen->SetBorderSize(0);
  legend1gen->SetFillStyle(0);
  legend1gen->AddEntry(hist_mc1_signalGen,"MC1 signal","l");
  legend1gen->AddEntry(hist_data_truth,"data truth","l");
  legend1gen->Draw();

  canvas->cd(3);
  hist_data_total->SetLineColor(kRed);
  hist_data_total->SetMarkerColor(kRed);
  hist_data_total->Draw();
  hist_mc2_total->SetLineColor(kBlue);
  hist_mc2_total->Draw("HIST SAME");
  hist_mc2_bgr->SetLineColor(kBlack);
  hist_mc2_bgr->SetFillColor(kGray);
  hist_mc2_bgr->SetFillStyle(1001);
  hist_mc2_bgr->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.7,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_data_total,"data","pl");
  legend2->AddEntry(hist_mc2_total,"MC2","l");
  legend2->Draw();

  canvas->cd(4);
  hist_mc2_signalGen->SetLineColor(kBlue);
  hist_mc2_signalGen->Draw("HIST");
  hist_data_truth->Draw("HIST SAME");
  TLegend *legend2gen=new TLegend(0.25,0.7,0.7,0.9);
  legend2gen->SetBorderSize(0);
  legend2gen->SetFillStyle(0);
  legend2gen->AddEntry(hist_mc2_signalGen,"MC2 signal","l");
  legend2gen->AddEntry(hist_data_truth,"data truth","l");
  legend2gen->Draw();
}
