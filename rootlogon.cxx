void rootlogon(void) {

  // if(gROOT->GetVersionInt()<51800) return;

  TStyle *square=new TStyle("square","square");
  gROOT->SetStyle("Plain");

  gStyle->Copy(*square);
  gROOT->SetStyle("square");
  
  gStyle->SetPadBorderSize(0.0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadBottomMargin(0.20);
  gStyle->SetPadLeftMargin(0.23);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetHistLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.8);
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleSize(0.8,"P");
  gStyle->SetTitleSize(0.07,"xy");
  gStyle->SetLabelSize(0.07,"xy");
  gStyle->SetLabelOffset(0.01,"xy");
  gStyle->SetNdivisions(505,"xy");
  double titleBorder=0.005;
  gStyle->SetTitleX(gStyle->GetPadLeftMargin());
  gStyle->SetTitleY(1.0-titleBorder);
  gStyle->SetTitleBorderSize(0.0);
  gStyle->SetTitleW(1.0-gStyle->GetPadLeftMargin()
                    -gStyle->GetPadRightMargin());
  gStyle->SetTitleH(gStyle->GetPadTopMargin()-2.*titleBorder);
  gStyle->SetGridColor(kGray);

  int d=250;
  // size of one plot for 2x2 is 250
  // but standard root canvas should have size of about 500
  gStyle->SetCanvasDefW(2.0*d/(1.-gStyle->GetPadRightMargin()
                           -gStyle->GetPadLeftMargin()));
  gStyle->SetCanvasDefH(2.0*d/(1.-gStyle->GetPadBottomMargin()
                           -gStyle->GetPadTopMargin()));

  gROOT->ForceStyle();

}
