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

using namespace std;

void exercise7(void) {
  //==============================================================
  // Unfolding exercise 7
  //  unfold by simple matrix inversion
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
  int nBin=10;
  double r0=0.;
  double r1=1.0;
  // data
  TH1D *hist_data_rec=new TH1D("hist_data_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_data_truth=new TH1D("hist_data_truth",";r_{gen}",nBin,r0,r1);
  // MC1
  TH1D *hist_mc1_rec=new TH1D("hist_mc1_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBin,r0,r1);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // MC2
  TH1D *hist_mc2_rec=new TH1D("hist_mc2_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc2_gen=new TH1D("hist_mc2_gen",";r_{gen}",nBin,r0,r1);

  // derived histograms
  // matrix of probabilities
  TH2D *hist_mc1_prob=new TH2D("hist_mc1_prob",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // unfolding results
  TH1D *hist_data_INVERT=new TH1D("hist_data_INVERT",";r_{INVERT}",nBin,r0,r1);
  TH2D *hist_data_INVERTrho=new TH2D
    ("hist_data_INVERTrho",";r_{INVERT};r_{INVERT}",nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc1_INVERT=new TH1D("hist_mc1_INVERT",";r_{INVERT}",nBin,r0,r1);
  TH2D *hist_mc1_INVERTrho=new TH2D
    ("hist_mc1_INVERTrho",";r_{INVERT};r_{INVERT}",nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc2_INVERT=new TH1D("hist_mc2_INVERT",";r_{INVERT}",nBin,r0,r1);
  TH2D *hist_mc2_INVERTrho=new TH2D
    ("hist_mc2_INVERTrho",";r_{INVERT};r_{INVERT}",nBin,r0,r1,nBin,r0,r1);

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

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("+hist_mc1_gen","0.95",bgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recgen","rRec:rGen",recsigWeight,
		  "",nEvent,firstEvent);
    tree->Project("+hist_mc1_recgen","rRec:0.95",recbgrWeight,
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
    TString bgrWeight="w*isBgr";
    TString recWeight="w*isTrig";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("+hist_mc2_gen","0.95",bgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc2_rec","rRec",recWeight,"",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // fill derived histograms
  // matrix of probabilities
  TMatrixD A(nBin,nBin);
  for(int iGen=1;iGen<=nBin;iGen++) {
    double nGen=hist_mc1_gen->GetBinContent(iGen);
    if(nGen>0.) {
      for(int iRec=1;iRec<=nBin;iRec++) {
	double c=hist_mc1_recgen->GetBinContent(iGen,iRec);
	hist_mc1_prob->SetBinContent(iGen,iRec,c/nGen);
	A(iRec-1,iGen-1)=c/nGen;
      }
    }
  }
  TMatrixD Ainv(TMatrixD::kInverted,A);

  // fill vectors with rec for data.mc1,mc2
  TVectorD data(nBin),mc1(nBin),mc2(nBin);
  for(int iRec=1;iRec<=nBin;iRec++) {
    data[iRec-1]=hist_data_rec->GetBinContent(iRec);
    mc1[iRec-1]=hist_mc1_rec->GetBinContent(iRec);
    mc2[iRec-1]=hist_mc2_rec->GetBinContent(iRec);
  }
  
  // calculate diagonal covariance matrix of rec-level uncertainties
  TMatrixDSym dataRecEmatrix(nBin),mc1RecEmatrix(nBin),mc2RecEmatrix(nBin);
  for(int iRec=1;iRec<=nBin;iRec++) {
    int i=iRec-1;
    dataRecEmatrix[i][i]=TMath::Power(hist_data_rec->GetBinError(iRec),2.);
    mc1RecEmatrix[i][i]=TMath::Power(hist_mc1_rec->GetBinError(iRec),2.);
    mc2RecEmatrix[i][i]=TMath::Power(hist_mc2_rec->GetBinError(iRec),2.);
  }
  // calculate unfolding results and covariance matrices
  TVectorD dataUNFOLD(Ainv*data);
  TVectorD mc1UNFOLD(Ainv*mc1);
  TVectorD mc2UNFOLD(Ainv*mc2);
  TMatrixD AinvT(TMatrixD::kTransposed,Ainv);
  TMatrixD dataUNFOLDcov(Ainv,TMatrixD::kMult,dataRecEmatrix*AinvT);
  TMatrixD mc1UNFOLDcov(Ainv,TMatrixD::kMult,mc1RecEmatrix*AinvT);
  TMatrixD mc2UNFOLDcov(Ainv,TMatrixD::kMult,mc2RecEmatrix*AinvT);

  // fill results into histograms
  for(int iGen=1;iGen<=nBin;iGen++) {
    int i=iGen-1;
    hist_data_INVERT->SetBinContent(iGen,dataUNFOLD[i]);
    hist_data_INVERT->SetBinError(iGen,TMath::Sqrt(dataUNFOLDcov[i][i]));
    hist_mc1_INVERT->SetBinContent(iGen,mc1UNFOLD[i]);
    hist_mc1_INVERT->SetBinError(iGen,TMath::Sqrt(mc1UNFOLDcov[i][i]));
    hist_mc2_INVERT->SetBinContent(iGen,mc2UNFOLD[i]);
    hist_mc2_INVERT->SetBinError(iGen,TMath::Sqrt(mc2UNFOLDcov[i][i]));
    for(int jGen=1;jGen<=nBin;jGen++) {
      int j=jGen-1;
      hist_data_INVERTrho->SetBinContent
	(iGen,jGen,dataUNFOLDcov[i][j]/
	 TMath::Sqrt(dataUNFOLDcov[i][i]*dataUNFOLDcov[j][j]));
      hist_mc1_INVERTrho->SetBinContent
	(iGen,jGen,mc1UNFOLDcov[i][j]/
	 TMath::Sqrt(mc1UNFOLDcov[i][i]*mc1UNFOLDcov[j][j]));
      hist_mc2_INVERTrho->SetBinContent
	(iGen,jGen,mc2UNFOLDcov[i][j]/
	 TMath::Sqrt(mc2UNFOLDcov[i][i]*mc2UNFOLDcov[j][j]));
    }
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("matrix inversion method","",900,600);
  canvas->Divide(3,2);

  canvas->cd(1);
  //hist_data_INVERT->GetXaxis()->SetRangeUser(0.,0.89);
  hist_data_INVERT->SetLineColor(kRed);
  hist_data_INVERT->SetMarkerColor(kRed);
  hist_data_INVERT->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9,"INVERT method");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_INVERT,"data unfolded","l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->Draw();

  TVirtualPad *p4=canvas->cd(4);
  p4->SetRightMargin(0.15);
  hist_data_INVERTrho->Draw("COLZ");
  hist_data_INVERTrho->Draw("BOX SAME");

  canvas->cd(2);
  //hist_mc1_INVERT->GetXaxis()->SetRangeUser(0.,0.89);
  hist_mc1_INVERT->SetLineColor(kRed);
  hist_mc1_INVERT->SetMarkerColor(kRed);
  hist_mc1_INVERT->Draw();
  hist_mc1_gen->SetLineColor(kBlue);
  hist_mc1_gen->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9,"INVERT method");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_mc1_INVERT,"MC1 unfolded","l");
  legend2->AddEntry(hist_mc1_gen,"MC1 truth","l");
  legend2->Draw();

  TVirtualPad *p5=canvas->cd(5);
  p5->SetRightMargin(0.15);
  hist_mc1_INVERTrho->Draw("COLZ");
  hist_mc1_INVERTrho->Draw("BOX SAME");

  canvas->cd(3);
  //hist_mc2_INVERT->GetXaxis()->SetRangeUser(0.,0.89);
  //hist_mc2_INVERT->GetYaxis()->SetRangeUser(-500.,2000.);
  hist_mc2_INVERT->SetLineColor(kRed);
  hist_mc2_INVERT->SetMarkerColor(kRed);
  hist_mc2_INVERT->Draw();
  hist_mc2_gen->SetLineColor(kBlue);
  hist_mc2_gen->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9,"INVERT method");
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_mc2_INVERT,"MC2 unfolded","l");
  legend3->AddEntry(hist_mc2_gen,"MC2 truth","l");
  legend3->Draw();

  TVirtualPad *p6=canvas->cd(6);
  p6->SetRightMargin(0.15);
  hist_mc2_INVERTrho->Draw("COLZ");
  hist_mc2_INVERTrho->Draw("BOX SAME");

}
