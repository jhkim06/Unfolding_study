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
#include <TUnfold.h>

using namespace std;

void exercise8910(void) {
  //==============================================================
  // Unfolding exercise 3,7,8,9,10
  //  unfold (data only) using various methods
  //  and compare the results
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
  int nBinRec=100;
  double r0Rec=0.;
  double r1Rec=1.0;
  // data
  TH1D *hist_data_rec[NREPLICA];
  for(int iReplica=0;iReplica<NREPLICA;iReplica++) {
    hist_data_rec[iReplica]=new TH1D
      (TString::Format("hist_data_rec%d",iReplica),";r_{rec}",nBin,r0,r1);
  }
  TH1D *hist_data_recTUNF0=new TH1D("hist_data_recTUNF0",";r_{rec}",
				   nBinRec,r0Rec,r1Rec);
  TH1D *hist_data_truth=new TH1D("hist_data_truth",";r_{gen}",nBin,r0,r1);
  // MC1
  TH1D *hist_mc1_rec=new TH1D("hist_mc1_rec",";r_{rec}",nBin,r0,r1);
  TH1D *hist_mc1_gen=new TH1D("hist_mc1_gen",";r_{gen}",nBin,r0,r1);
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen",";r_{gen};r_{rec}",
				 nBin,r0,r1,nBin,r0,r1);
  TH2D *hist_mc1_recgenTUNF0=new TH2D("hist_mc1_recgenTUNF0",";r_{gen};r_{rec}",
				    nBin,r0,r1,nBinRec,r0Rec,r1Rec);
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

  // D'Agostini and iterated method

  TH1D *hist_data_ITER[NITERATION];
  TH2D *hist_data_ITERrho[NITERATION];
  for(int iIteration=0;iIteration<NITERATION;iIteration++) {
    hist_data_ITER[iIteration]=new TH1D
      (TString::Format("hist_data_ITER%d",iIteration+1),
       iIteration ? 
       TString::Format("Iterative n(iteration)=%d;r_{ITER}",iIteration+1) :
       "\"Bayesian\" (D'Agostini);r_{AGO}",
       nBin,r0,r1);
    hist_data_ITERrho[iIteration]=new TH2D
      (TString::Format("hist_data_ITERrho%d",iIteration+1),
       ";r_{ITER};r_{ITER}",nBin,r0,r1,nBin,r0,r1);
  }

    //  correction factor, simple Bin-By-Bin (BBB)
  TH1D *hist_mc1_fBBB=new TH1D("hist_mc1_BBB",";r",nBin,r0,r1);
  // results for BBB
  TH1D *hist_data_BBB=new TH1D("hist_data_BBB","Bin-by-bin correction;r_{BBB}",nBin,r0,r1);

  // matrix inversion
  TH1D *hist_data_INVERT=new TH1D("hist_data_INVERT",
				  "Matrix inversion;r_{INVERT}",nBin,r0,r1);
  // fraction fit
  TH1D *hist_data_TUNF0=new TH1D("hist_data_TUNF0",
				  "Fraction fit (TUnfold) nRec=100;r_{fit}",
				 nBin,r0,r1);

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
    tree->Project("hist_data_recTUNF0","rRec",recWeight,"",nEvent,firstEvent);
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
    tree->Project("hist_mc1_rec","rRec",recWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_gen","rGen",genWeight,"",nEvent,firstEvent);
    tree->Project("+hist_mc1_gen","0.95",bgrWeight,"",nEvent,firstEvent);
    tree->Project("hist_mc1_recgen","rRec:rGen",recsigWeight,
		  "",nEvent,firstEvent);
    tree->Project("+hist_mc1_recgen","rRec:0.95",recbgrWeight,
		  "",nEvent,firstEvent);

    // event with generator and reco information
    tree->Project("hist_mc1_recgenTUNF0",
		  "rRec:rGen",recsigWeight,"",nEvent,firstEvent);
    // background event with reco information
    tree->Project("+hist_mc1_recgenTUNF0",
		  "rRec:0.95",recbgrWeight,"",nEvent,firstEvent);
    // events with generator but no reco information
    tree->Project("+hist_mc1_recgenTUNF0",
		  "-1.0:rGen",norecsigWeight,"",nEvent,firstEvent);
    // background event without reco information
    tree->Project("+hist_mc1_recgenTUNF0",
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
    histogramDir->cd();

    // clean up
    delete tree;
  }

  // calculate correction factors BBB
  hist_mc1_fBBB->Divide(hist_mc1_gen,hist_mc1_rec);

  // results for BBB
  hist_data_BBB->Multiply(hist_mc1_fBBB,hist_data_rec[0]);

  // matrix inversion method
  TMatrixDSym dataRecEmatrix(nBin);
  for(int iRec=1;iRec<=nBin;iRec++) {
    int i=iRec-1;
    dataRecEmatrix[i][i]=TMath::Power(hist_data_rec[0]->GetBinError(iRec),2.);
  }

  // fill derived histograms
  // matrix of probabilities and 
  // histogram of efficiencies
  TMatrixD A(nBin,nBin);
  for(int iGen=1;iGen<=nBin;iGen++) {
    double nGen=hist_mc1_gen->GetBinContent(iGen);
    if(nGen>0.) {
      double recSum=0.0;
      for(int iRec=1;iRec<=nBin;iRec++) {
	double c=hist_mc1_recgen->GetBinContent(iGen,iRec);
	hist_mc1_prob->SetBinContent(iGen,iRec,c/nGen);
	A(iRec-1,iGen-1)=c/nGen;
	recSum +=c;
      }
      hist_mc1_eff->SetBinContent(iGen,recSum/nGen);
    }
  }
  TMatrixD Ainv(TMatrixD::kInverted,A);

  TVectorD data(nBin);
  for(int iRec=1;iRec<=nBin;iRec++) {
    data[iRec-1]=hist_data_rec[0]->GetBinContent(iRec);
  }


  // inverssion method
  // calculate unfolding results and covariance matrices
  TVectorD dataUNFOLD(Ainv*data);
  TMatrixD AinvT(TMatrixD::kTransposed,Ainv);
  TMatrixD dataUNFOLDcov(Ainv,TMatrixD::kMult,dataRecEmatrix*AinvT);

  // fill results into histograms
  for(int iGen=1;iGen<=nBin;iGen++) {
    int i=iGen-1;
    hist_data_INVERT->SetBinContent(iGen,dataUNFOLD[i]);
    hist_data_INVERT->SetBinError(iGen,TMath::Sqrt(dataUNFOLDcov[i][i]));
  }

  // Bayesian and iterative
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

  // TUnfold
  TUnfold unfold(hist_mc1_recgenTUNF0,
		 TUnfold::kHistMapOutputHoriz,
		 TUnfold::kRegModeSize,
		 TUnfold::kEConstraintNone);
  unfold.SetInput(hist_data_recTUNF0);
  unfold.DoUnfold(0.0);
  unfold.GetOutput(hist_data_TUNF0);

  // make plot
  TCanvas *canvas=new TCanvas("Comparison unfolding methods","",900,600);
  canvas->Divide(3,2);

  int iter1=0;
  int iter2=NITERATION-1;

  double ymax=3500;

  canvas->cd(1);
  hist_data_BBB->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_BBB->SetLineColor(kRed);
  hist_data_BBB->SetMarkerColor(kRed);
  hist_data_BBB->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  hist_mc1_gen->SetLineColor(kBlack);
  hist_mc1_gen->SetLineStyle(2);
  TH1 *mc1_truth=hist_mc1_gen->DrawCopy("HIST SAME");
  TLegend *legend1=new TLegend(0.25,0.7,0.85,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(hist_data_ITER[iter1],"data bin-by-bin","l");
  legend1->AddEntry(hist_data_truth,"data truth","l");
  legend1->AddEntry(mc1_truth,"MC1 truth","pl");
  legend1->Draw();

  canvas->cd(2);
  hist_data_ITER[iter1]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_ITER[iter1]->SetLineColor(kRed);
  hist_data_ITER[iter1]->SetMarkerColor(kRed);
  hist_data_ITER[iter1]->Draw();
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.25,0.7,0.85,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(hist_data_ITER[iter1],"data D'Agostini","l");
  legend2->AddEntry(hist_data_truth,"data truth","l");
  legend2->AddEntry(mc1_truth,"MC1 truth","pl");
  legend2->Draw();

  canvas->cd(3);
  hist_data_ITER[iter2]->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_ITER[iter2]->SetLineColor(kRed);
  hist_data_ITER[iter2]->SetMarkerColor(kRed);
  hist_data_ITER[iter2]->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend3=new TLegend(0.25,0.7,0.85,0.9);
  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->AddEntry(hist_data_ITER[iter2],
		    TString::Format("data iteration %d",iter2+1),"l");
  legend3->AddEntry(hist_data_truth,"data truth","l");
  legend3->AddEntry(mc1_truth,"MC1 truth","pl");
  legend3->Draw();

  canvas->cd(4);
  hist_data_INVERT->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_INVERT->SetLineColor(kRed);
  hist_data_INVERT->SetMarkerColor(kRed);
  hist_data_INVERT->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend4=new TLegend(0.25,0.7,0.85,0.9);
  legend4->SetBorderSize(0);
  legend4->SetFillStyle(0);
  legend4->AddEntry(hist_data_INVERT,
		    "data matrix inversion","l");
  legend4->AddEntry(hist_data_truth,"data truth","l");
  legend4->AddEntry(mc1_truth,"MC1 truth","pl");
  legend4->Draw();

  canvas->cd(5);
  hist_data_TUNF0->GetYaxis()->SetRangeUser(0.,ymax);
  hist_data_TUNF0->SetLineColor(kRed);
  hist_data_TUNF0->SetMarkerColor(kRed);
  hist_data_TUNF0->Draw();
  hist_data_truth->Draw("HIST SAME");
  mc1_truth->Draw("HIST SAME");
  TLegend *legend5=new TLegend(0.25,0.7,0.85,0.9);
  legend5->SetBorderSize(0);
  legend5->SetFillStyle(0);
  legend5->AddEntry(hist_data_TUNF0,
		    "data TUnfold #tau=0","l");
  legend5->AddEntry(hist_data_truth,"data truth","l");
  legend5->AddEntry(mc1_truth,"MC1 truth","pl");
  legend5->Draw();

}
