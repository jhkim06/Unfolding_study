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

void exercise8(void) {
  //==============================================================
  // Unfolding exercise 8
  //  unfold using D'Agostini method
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
  TH1D *hist_mc1_eff=new TH1D("hist_mc1_eff",";r_{gen}",nBin,r0,r1);

  // unfolding results
  TH1D *hist_data_AGO=new TH1D("hist_data_AGO",";r_{AGO}",nBin,r0,r1);
  TH2D *hist_data_AGOrho=new TH2D
    ("hist_data_AGOrho",";r_{AGO};r_{AGO}",nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc1_AGO=new TH1D("hist_mc1_AGO",";r_{AGO}",nBin,r0,r1);
  TH2D *hist_mc1_AGOrho=new TH2D
    ("hist_mc1_AGOrho",";r_{AGO};r_{AGO}",nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc2_AGO=new TH1D("hist_mc2_AGO",";r_{AGO}",nBin,r0,r1);
  TH2D *hist_mc2_AGOrho=new TH2D
    ("hist_mc2_AGOrho",";r_{AGO};r_{AGO}",nBin,r0,r1,nBin,r0,r1);

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
  // matrix of probabilities and 
  // histogram of efficiencies
  for(int iGen=1;iGen<=nBin;iGen++) {
    double nGen=hist_mc1_gen->GetBinContent(iGen);
    if(nGen>0.) {
      double recSum=0.0;
      for(int iRec=1;iRec<=nBin;iRec++) {
	double c=hist_mc1_recgen->GetBinContent(iGen,iRec);
	hist_mc1_prob->SetBinContent(iGen,iRec,c/nGen);
	recSum +=c;
      }
      hist_mc1_eff->SetBinContent(iGen,recSum/nGen);
    }
  }
  // unfold data, then MC1, then MC2
  TH1D *input[3];
  TH1D *output[3];
  TH2D *rho[3];
  input[0]=hist_data_rec;
  output[0]=hist_data_AGO;
  rho[0]=hist_data_AGOrho;
  input[1]=hist_mc1_rec;
  output[1]=hist_mc1_AGO;
  rho[1]=hist_mc1_AGOrho;
  input[2]=hist_mc2_rec;
  output[2]=hist_mc2_AGO;
  rho[2]=hist_mc2_AGOrho;

  // repeat unfolding for three inputs
  double *recRatio=new double[nBin+1];
  double *recRatioErrorSquared=new double[nBin+1];

  for(int iDS=0;iDS<3;iDS++) {
    for(int iRec=1;iRec<=nBin;iRec++) {
      double nRecMC=hist_mc1_rec->GetBinContent(iRec);
      recRatio[iRec]=input[iDS]->GetBinContent(iRec)/nRecMC;
      recRatioErrorSquared[iRec]=TMath::Power
	(input[iDS]->GetBinError(iRec)/nRecMC,2.);
    }
    for(int jGen=1;jGen<=nBin;jGen++) {
      double fj=0.0;
      for(int iRec=1;iRec<=nBin;iRec++) {
	fj += hist_mc1_prob->GetBinContent(jGen,iRec)*recRatio[iRec];
      }
      double xj=hist_mc1_gen->GetBinContent(jGen)*fj/
	hist_mc1_eff->GetBinContent(jGen);
      output[iDS]->SetBinContent(jGen,xj);
      // loopp such that the diagonal is calculated first
      for(int kGen=jGen;kGen>=0;kGen--) {
	double v_jk=0.0;
	for(int iRec=1;iRec<=nBin;iRec++) {
	  v_jk += hist_mc1_prob->GetBinContent(jGen,iRec)*
	    hist_mc1_prob->GetBinContent(kGen,iRec)*
	    recRatioErrorSquared[iRec];
	}
	v_jk *=
	  hist_mc1_gen->GetBinContent(jGen)*
	  hist_mc1_gen->GetBinContent(kGen)/
	  (hist_mc1_eff->GetBinContent(jGen)*
	   hist_mc1_eff->GetBinContent(kGen));
	  
	if(jGen==kGen) {
	  output[iDS]->SetBinError(jGen,TMath::Sqrt(v_jk));
	}
	double e_j=output[iDS]->GetBinError(jGen);
	double e_k=output[iDS]->GetBinError(kGen);
	double rho_jk=v_jk/(e_j*e_k);
	rho[iDS]->SetBinContent(jGen,kGen,rho_jk);
	rho[iDS]->SetBinContent(kGen,jGen,rho_jk);
      }
    }
  }
  delete [] recRatio;
  delete [] recRatioErrorSquared;

  // make plot
  TCanvas *canvas=
    new TCanvas("D'Agostini method","",900,600);
  canvas->Divide(3,2);

  double ymax=3500.;
  canvas->cd(1);
  hist_data_AGO->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_AGO->SetLineColor(kRed);
  hist_data_AGO->SetMarkerColor(kRed);
  hist_data_AGO->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_AGO,"data unfolded","l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  TVirtualPad *p4=canvas->cd(4);
  p4->SetRightMargin(0.15);
  hist_data_AGOrho->Draw("COLZ");
  hist_data_AGOrho->Draw("BOX SAME");

  canvas->cd(2);
  hist_mc1_AGO->GetYaxis()->SetRangeUser(0.,ymax);
  hist_mc1_AGO->SetLineColor(kRed);
  hist_mc1_AGO->SetMarkerColor(kRed);
  hist_mc1_AGO->Draw();
  hist_mc1_gen->SetLineColor(kBlue);
  hist_mc1_gen->SetLineStyle(1);
  hist_mc1_gen->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_mc1_AGO,"MC1 unfolded","l");
  legend2->AddEntry(hist_mc1_gen,"MC1 truth","l");
  legend2->Draw();

  TVirtualPad *p5=canvas->cd(5);
  p5->SetRightMargin(0.15);
  hist_mc1_AGOrho->Draw("COLZ");
  hist_mc1_AGOrho->Draw("BOX SAME");

  canvas->cd(3);
  hist_mc2_AGO->GetYaxis()->SetRangeUser(0.,ymax);
  hist_mc2_AGO->SetLineColor(kRed);
  hist_mc2_AGO->SetMarkerColor(kRed);
  hist_mc2_AGO->Draw();
  hist_mc2_gen->SetLineColor(kBlue);
  hist_mc2_gen->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_mc2_AGO,"MC2 unfolded","l");
  legend3->AddEntry(hist_mc2_gen,"MC2 truth","l");
  legend3->AddEntry(mc1_truth,"MC1 truth","pl");
  legend3->Draw();

  TVirtualPad *p6=canvas->cd(6);
  p6->SetRightMargin(0.15);
  hist_mc2_AGOrho->Draw("COLZ");
  hist_mc2_AGOrho->Draw("BOX SAME");
}
