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
#include <TRandom3.h>

using namespace std;

void exercise9(void) {
  //==============================================================
  // Unfolding exercise 9
  //  unfold (data only) using iterative method
  //  use data replicas to estimate covariances
  // show result after 1,10,100 iterations
  //==============================================================

  // the event samples are divided into subsets
  //  for most exercises we divide the event samples into 30 subsets
  //  and use only the first subset
  int nSubset=30;
  int iSubset=0;

  // 200 replicas of the data are used to calculate the covariance matrix
#define NREPLICA 200

  // 100 iterations of the unfolding
#define NITERATION 100

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
  TH1D *hist_data_rec[NREPLICA];
  for(int iReplica=0;iReplica<NREPLICA;iReplica++) {
    hist_data_rec[iReplica]=new TH1D
      (TString::Format("hist_data_rec%d",iReplica),";r_{rec}",nBin,r0,r1);
  }
  TH1D *hist_data_truth=new TH1D("hist_data_truth",";r_{gen}",nBin,r0,r1);
  // MC1
  TH1D *hist_mc1_rec=new TH1D("hist_mc1_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBin,r0,r1);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  // derived histograms
  // matrix of probabilities
  TH2D *hist_mc1_prob=new TH2D("hist_mc1_prob",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  TH1D *hist_mc1_eff=new TH1D("hist_mc1_eff",";r_{gen}",nBin,r0,r1);

  // unfolding results for all replicas
  TH1D *hist_data_AGO[NREPLICA];
  for(int iReplica=0;iReplica<NREPLICA;iReplica++) {
    hist_data_AGO[iReplica]=new TH1D
      (TString::Format("hist_data_AGO%d",iReplica),";r_{AGO}",nBin,r0,r1);
  }

  TH1D *hist_data_ITER[NITERATION];
  TH2D *hist_data_ITERrho[NITERATION];
  for(int iIteration=0;iIteration<NITERATION;iIteration++) {
    hist_data_ITER[iIteration]=new TH1D
      (TString::Format("hist_data_ITER%d",iIteration+1),";r_{ITER}",nBin,r0,r1);
    hist_data_ITERrho[iIteration]=new TH2D
      (TString::Format("hist_data_ITERrho%d",iIteration+1),
       ";r_{ITER};r_{ITER}",nBin,r0,r1,nBin,r0,r1);
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
    tree->Project("hist_data_rec0","rRec",recWeight,"",nEvent,firstEvent);
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

  // produce data replicas with extra stat. fluctuations
  TRandom3 rnd;
  for(int iRec=1;iRec<=nBin;iRec++) {
    double c0=hist_data_rec[0]->GetBinContent(iRec);
    double e0=hist_data_rec[0]->GetBinError(iRec);
    double n0=c0*c0/(e0*e0);
    double a=e0*e0/c0;
    for(int iReplica=1;iReplica<NREPLICA;iReplica++) {
      double ni=rnd.Poisson(n0);
      double ci=a*ni;
      double ei=a*TMath::Sqrt(ni);
      hist_data_rec[iReplica]->SetBinContent(iRec,ci);
      hist_data_rec[iReplica]->SetBinError(iRec,ei);
    }
  }

  // iteration zero: use MC1 truth as start value
  for(int iReplica=0;iReplica<NREPLICA;iReplica++) {
    for(int jGen=1;jGen<=nBin;jGen++) {
      double nGenJ=hist_mc1_gen->GetBinContent(jGen);
      hist_data_AGO[iReplica]->SetBinContent(jGen,nGenJ);
    }
  }

  for(int iIteration=0;iIteration<NITERATION;iIteration++) {
    double *recRatio=new double[nBin+1];
    // repeat unfolding for all data subsets
    for(int iReplica=0;iReplica<NREPLICA;iReplica++) {
      for(int iRec=1;iRec<=nBin;iRec++) {
	double nRecMC=0.0;
	for(int jGen=1;jGen<=nBin;jGen++) {
	  nRecMC += hist_mc1_prob->GetBinContent(jGen,iRec)*
	    hist_data_AGO[iReplica]->GetBinContent(jGen);
	}
	recRatio[iRec]=hist_data_rec[iReplica]->GetBinContent(iRec)/nRecMC;
      }
      for(int jGen=1;jGen<=nBin;jGen++) {
	double fj=0.0;
	for(int iRec=1;iRec<=nBin;iRec++) {
	  fj += hist_mc1_prob->GetBinContent(jGen,iRec)*recRatio[iRec];
	}
	double xj=hist_data_AGO[iReplica]->GetBinContent(jGen)*fj/
	  hist_mc1_eff->GetBinContent(jGen);
	hist_data_AGO[iReplica]->SetBinContent(jGen,xj);
      }
    }
    delete [] recRatio;
    for(int jGen=1;jGen<=nBin;jGen++) {
      // unfolding central result: from first data subset
      hist_data_ITER[iIteration]->SetBinContent
	(jGen,hist_data_AGO[0]->GetBinContent(jGen));
      // unfolding covariance: estimate from data replicas
      for(int kGen=jGen;kGen>=1;kGen--) {
	double s[2][2]={{0.,0.},{0.,0.}};
	for(int iReplica=1;iReplica<NREPLICA;iReplica++) {
	  double yj=hist_data_AGO[iReplica]->GetBinContent(jGen);
	  double yk=hist_data_AGO[iReplica]->GetBinContent(kGen);
	  s[0][0]+=1.;
	  s[1][0]+=yj;
	  s[0][1]+=yk;
	  s[1][1]+=yj*yk;
	}
	double meanJ=s[1][0]/s[0][0];
	double meanK=s[0][1]/s[0][0];
	double rmsJK=s[1][1]/s[0][0]-meanJ*meanK;
	if(jGen==kGen) {
	  hist_data_ITER[iIteration]->SetBinError(jGen,TMath::Sqrt(rmsJK));
	}
	double rhoJK=rmsJK/
	  (hist_data_ITER[iIteration]->GetBinError(jGen)*
	   hist_data_ITER[iIteration]->GetBinError(kGen));
	hist_data_ITERrho[iIteration]->SetBinContent(jGen,kGen,rhoJK);
	hist_data_ITERrho[iIteration]->SetBinContent(kGen,jGen,rhoJK);
      }
    }
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("Iterative method","",900,600);
  canvas->Divide(3,2);

  int iter1=0;
  int iter2=9;
  int iter3=NITERATION-1;

  double ymax=3500;

  canvas->cd(1);
  hist_data_ITER[iter1]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_ITER[iter1]->SetLineColor(kRed);
  hist_data_ITER[iter1]->SetMarkerColor(kRed);
  hist_data_ITER[iter1]->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_ITER[iter1],
		    TString::Format("data iteration %d",iter1+1),"l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  TVirtualPad *p4=canvas->cd(4);
  p4->SetRightMargin(0.15);
  hist_data_ITERrho[iter1]->Draw("COLZ");
  hist_data_ITERrho[iter1]->Draw("BOX SAME");

  canvas->cd(2);
  hist_data_ITER[iter2]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_ITER[iter2]->SetLineColor(kRed);
  hist_data_ITER[iter2]->SetMarkerColor(kRed);
  hist_data_ITER[iter2]->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_data_ITER[iter2],
		    TString::Format("data iteration %d",iter2+1),"l");
  legend2->AddEntry(hist_data_truth,"data truth","l");
  legend2->AddEntry(mc1_truth,"MC1 truth","pl");
  legend2->Draw();

  TVirtualPad *p5=canvas->cd(5);
  p5->SetRightMargin(0.15);
  hist_data_ITERrho[iter2]->Draw("COLZ");
  hist_data_ITERrho[iter2]->Draw("BOX SAME");

  canvas->cd(3);
  //hist_data_ITER[iter3]->GetXaxis()->SetRangeUser(0.,0.89);
  //hist_data_ITER[iter3]->GetYaxis()->SetRangeUser(-500.,2000.);
  hist_data_ITER[iter3]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_ITER[iter3]->SetLineColor(kRed);
  hist_data_ITER[iter3]->SetMarkerColor(kRed);
  hist_data_ITER[iter3]->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9,"AGO method");
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_data_ITER[iter3],
		    TString::Format("data iteration %d",iter3+1),"l");
  legend3->AddEntry(hist_data_truth,"data truth","l");
  legend3->AddEntry(mc1_truth,"MC1 truth","pl");
  legend3->Draw();

  TVirtualPad *p6=canvas->cd(6);
  p6->SetRightMargin(0.15);
  hist_data_ITERrho[iter3]->Draw("COLZ");
  hist_data_ITERrho[iter3]->Draw("BOX SAME");

}
