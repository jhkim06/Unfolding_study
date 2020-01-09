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

void exercise3a(void) {
  //==============================================================
  // Unfolding exercise 3a
  //  unfold data, MC1,MC2
  //  using the correction factors obtained from MC1
  //==============================================================

  // the event samples are divided into subsets
  //  for most exercises we divide the event samples into 100 subsets
  //  and use only the first subset
  int nSubset=100;
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
  // data
  TH1D *hist_data_rec=new TH1D("hist_data_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_data_truth=new TH1D("hist_data_truth",";r_{gen}",nBin,r0,r1);
  // MC1
  TH1D *hist_mc1_rec=new TH1D("hist_mc1_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_recsig=new TH1D("hist_mc1_recsig",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_recbgr=new TH1D("hist_mc1_recbgr",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBin,r0,r1);
  // MC2
  TH1D *hist_mc2_rec=new TH1D("hist_mc2_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc2_gen=new TH1D("hist_mc2_gen",";r_{gen}",nBin,r0,r1);
  // unfolding result
  //  correction factor, simple Bin-By-Bin (BBB)
  TH1D *hist_mc1_fBBB=new TH1D("hist_mc1_BBB",";r",nBin,r0,r1);
  //  correction factor, Bin-By-Bin with Background Subtraction (BBBBS)
  TH1D *hist_mc1_fBBBBS=new TH1D("hist_mc1_BBBBS",";r",nBin,r0,r1);
  // results for BBB
  TH1D *hist_data_BBB=new TH1D("hist_data_BBB",";r_{BBB}",nBin,r0,r1);
  TH1D *hist_mc1_BBB=new TH1D("hist_mc1_BBB",";r_{BBB}",nBin,r0,r1);
  TH1D *hist_mc2_BBB=new TH1D("hist_mc2_BBB",";r_{BBB}",nBin,r0,r1);
  // results for BBBBS
  TH1D *hist_data_BBBBS=new TH1D("hist_data_BBBBS",";r_{BBBBS}",nBin,r0,r1);
  TH1D *hist_mc1_BBBBS=new TH1D("hist_mc1_BBBBS",";r_{BBBBS}",nBin,r0,r1);
  TH1D *hist_mc2_BBBBS=new TH1D("hist_mc2_BBBBS",";r_{BBBBS}",nBin,r0,r1);

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

    // fill histogram
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

    // fill histogram
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
    TString recbgrWeight="w*isTrig*isBgr";
    TString recsigWeight="w*isTrig*(1-isBgr)";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recbgr","rRec",recbgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recsig","rRec",recsigWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_gen","rGen",genWeight,"",nEvent,firstEvent);

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
    TString recbgrWeight="w*isTrig*isBgr";
    TString recsigWeight="w*isTrig*(1-isBgr)";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_gen","rGen",genWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // calculate correction factors BBB
  hist_mc1_fBBB->Divide(hist_mc1_gen,hist_mc1_rec);
  // calculate correction factors BBBBS
  hist_mc1_fBBBBS->Divide(hist_mc1_gen,hist_mc1_recsig);

  // results for BBB
  hist_data_BBB->Multiply(hist_mc1_fBBB,hist_data_rec);
  hist_mc1_BBB->Multiply(hist_mc1_fBBB,hist_mc1_rec);
  hist_mc2_BBB->Multiply(hist_mc1_fBBB,hist_mc2_rec);

  // results for BBBBS
  hist_data_BBBBS->Add(hist_data_rec,hist_mc1_recbgr,1.,-1.);
  hist_mc1_BBBBS->Add(hist_mc1_rec,hist_mc1_recbgr,1.,-1.);
  hist_mc2_BBBBS->Add(hist_mc2_rec,hist_mc1_recbgr,1.,-1.);
  hist_data_BBBBS->Multiply(hist_mc1_fBBBBS);
  hist_mc1_BBBBS->Multiply(hist_mc1_fBBBBS);
  hist_mc2_BBBBS->Multiply(hist_mc1_fBBBBS);

  // make plot
  TCanvas *canvas=
    new TCanvas("bin-by-bin correction","",900,300);
  canvas->Divide(3,1);

  double ymax=2200.;

  canvas->cd(1);
  hist_data_BBB->SetLineColor(kRed);
  hist_data_BBB->SetMarkerColor(kRed);
  hist_data_BBB->Draw();
  hist_data_truth->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->SetLineWidth(4);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_BBB,"BBB corrected data","pl");
  legend1->AddEntry(hist_data_truth,"data truth","pl");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  canvas->cd(2);
  hist_mc1_BBB->SetLineColor(kRed);
  hist_mc1_BBB->SetMarkerColor(kRed);
  hist_mc1_BBB->Draw();
  hist_mc1_gen->GetYaxis()->SetRangeUser(0.,ymax);
  hist_mc1_gen->SetLineColor(kBlue);
  hist_mc1_gen->SetLineWidth(4);
  hist_mc1_gen->SetLineStyle(1);
  hist_mc1_gen->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_mc1_BBB,"BBB corrected MC1","pl");
  legend2->AddEntry(hist_mc1_gen,"MC1 truth","l");
  legend2->Draw();

  canvas->cd(3);
  hist_mc2_BBB->SetLineColor(kRed);
  hist_mc2_BBB->SetMarkerColor(kRed);
  hist_mc2_BBB->GetYaxis()->SetRangeUser(0.,2200.);
  hist_mc2_BBB->Draw();
  hist_mc2_gen->GetYaxis()->SetRangeUser(0.,ymax);
  hist_mc2_gen->SetLineColor(kBlue);
  hist_mc2_gen->SetLineWidth(4);
  hist_mc2_gen->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9);
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_mc2_BBB,"BBB corrected MC2","pl");
  legend3->AddEntry(hist_mc2_gen,"MC2 truth","l");
  legend3->AddEntry(mc1_truth,"MC1 truth","pl");
  legend3->Draw();


}
