// Author: Stefan Schmitt
// DESY, 10/08/11

//  Version 17.2, with XMNL support, bug fix with bin map creation,
//                isPeriodic option for neighbour bins
//
//  History:
//    Version 17.1, in parallel to changes in TUnfold
//    Version 17.0, initial version, numbered in parallel to TUnfold

//////////////////////////////////////////////////////////////////////////
//
//  TUnfoldBinning
//
//  This class serves as a container of analysis bins
//  analysis bins are specified by defining the axes of a distribution.
//  It is also possible to have unconnected analysis bins without axis.
//  Multiple TUnfoldBinning objects may be arranged in a tree,
//  such that a full tree structure of histograms and bins is supported
//
//  If you use this software, please consider the following citation
//       S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]
//
//  More documentation and updates are available on
//      http://www.desy.de/~sschmitt
//
//  Functionality
//
//  The class gives access to all analysis bins numbered in sequence.
//  Such a sequence of bins may be stored in a 1-dimension histogram.
//  Correlations between two TUnfoldBinning objects may be stored in
//  a 2-dimensional histogram. This type of ordering is required for
//  the TUnfold class.
//
//  In addition, it is possible to have root histograms, using the
//  axes as defined with the distributions. Underflow/overflow bins
//  can be included or excluded when mapping bins on root histograms.
//  In addition, it is possible to collapse one of more axes when going
//  from a N-dimensional distribution to a root histogram.
//
//////////////////////////////////////////////////////////////////////////

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "TUnfoldBinningXML.h"
#include <TVectorD.h>
#include <TAxis.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TIterator.h>
#include <iomanip>

// #define DEBUG

using namespace std;

ClassImp(TUnfoldBinningV17)

/********************* setup **************************/

void TUnfoldBinningV17::Initialize(Int_t nBins)
{
  // initialize variables
   parentNode=0;
   childNode=0;
   nextNode=0;
   prevNode=0;
   fAxisList=new TObjArray();
   fAxisLabelList=new TObjArray();
   fAxisList->SetOwner();
   fAxisLabelList->SetOwner();
   fHasUnderflow=0;
   fHasOverflow=0;
   fDistributionSize=nBins;
   fBinFactorFunction=0;
   fBinFactorConstant=1.0;
}

Int_t TUnfoldBinningV17::UpdateFirstLastBin(Bool_t startWithRootNode)
{
   // update fFirstBin and fLastBin members of this node and its children
   //   startWithRootNode: if true, start the update with the root node
   if(startWithRootNode) {
      return GetRootNode()->UpdateFirstLastBin(kFALSE);
   }
   if(GetPrevNode()) {
      // if this is not the first node in a sequence,
      // start with the end bin of the previous node
      fFirstBin=GetPrevNode()->GetEndBin();
   } else if(GetParentNode()) {
      // if this is the first node in a sequence but has a parent,
      // start with the end bin of the parent's distribution
      fFirstBin=GetParentNode()->GetStartBin()+
         GetParentNode()->GetDistributionNumberOfBins();
   } else {
      // if this is the top level node, the first bin number is 1
     fFirstBin=1;
     //  ... unless the top level node is the only node
     //  ... with dimension=1
     //  ... and there are no child nodes
     //  ... and there is an underflow bin
     if((!GetChildNode())&&(GetDistributionDimension()==1)&&
        (fHasUnderflow==1)) {
        fFirstBin=0;
     }
   }
   fLastBin=fFirstBin+fDistributionSize;
   // now update count for all children
   for(TUnfoldBinningV17 *node=childNode;node;node=node->nextNode) {
      fLastBin=node->UpdateFirstLastBin(kFALSE);
   }
   return fLastBin;
}

TUnfoldBinningV17::TUnfoldBinningV17
(const char *name,Int_t nBins,const char *binNames)
   : TNamed(name ? name : "",name ? name : "")
{
   // initialize a node with bins but without axis
   //   name: name of the node
   //   nBin: number of extra bins (could be zero)
   //   binNames: (optionally) names of the bins sepatared by ';'
   Initialize(nBins);
   if(binNames) {
      TString nameString(binNames);
      delete fAxisLabelList;
      fAxisLabelList=nameString.Tokenize(";");
   }
   UpdateFirstLastBin();
}

TUnfoldBinningV17::TUnfoldBinningV17
(const TAxis &axis,Int_t includeUnderflow,Int_t includeOverflow)
   : TNamed(axis.GetName(),axis.GetTitle())
{
   // create binning containing a distribution with one axis
   //     axis: the axis to represent
   //     includeUnderflow: include underflow bin
   //     includeOverflow: include overflow bin
   Initialize(0);
   AddAxis(axis,includeUnderflow,includeOverflow);
   UpdateFirstLastBin();
}

TUnfoldBinningV17::~TUnfoldBinningV17(void)
{
   // delete all children
   while(childNode) delete childNode;
   // remove this node from the tree
   if(GetParentNode() && (GetParentNode()->GetChildNode()==this)) {
      parentNode->childNode=nextNode;
   }
   if(GetPrevNode()) prevNode->nextNode=nextNode;
   if(GetNextNode()) nextNode->prevNode=prevNode;
   delete fAxisList; 
   delete fAxisLabelList; 
   if(fBinFactorFunction) {
      if(!dynamic_cast<TF1 *>(fBinFactorFunction)) 
         delete fBinFactorFunction;
   }
}

TUnfoldBinningV17 *TUnfoldBinningV17::AddBinning
(const char *name,Int_t nBins,const char *binNames)
{
  // add a binning as last daughter to this tree
  //   name: name of the node
  //   nBin: number of bins not belonging to a distribution (usually zero)
  //   binNames: (optionally) names of these bins sepatared by ';'
  return AddBinning(new TUnfoldBinning(name,nBins,binNames));
}

TUnfoldBinningV17 *TUnfoldBinningV17::AddBinning(TUnfoldBinningV17 *binning)
{
  // add a binning as last daughter to this tree
  //   binning: pointer to the new binning
  // return value: if succeeded, return "binning"
  //               otherwise return 0
  TUnfoldBinning *r=0;
  if(binning->GetParentNode()) {
     Error("AddBinning",
           "binning \"%s\" already has parent \"%s\", can not be added to %s",
	  (char *)binning->GetName(),
	  (char *)binning->GetParentNode()->GetName(),
	  (char *)GetName());
  } else if(binning->GetPrevNode()) {
    Error("AddBinning",
          "binning \"%s\" has previous node \"%s\", can not be added to %s",
	  (char *)binning->GetName(),
	  (char *)binning->GetPrevNode()->GetName(),
	  (char *)GetName());
  } else if(binning->GetNextNode()) {
    Error("AddBinning",
          "binning \"%s\" has next node \"%s\", can not be added to %s",
	  (char *)binning->GetName(),
	  (char *)binning->GetNextNode()->GetName(),
	  (char *)GetName());
  } else {
    r=binning;
    binning->parentNode=this;
    if(childNode) {
      TUnfoldBinningV17 *child=childNode;
      // find last child
      while(child->nextNode) {
	child=child->nextNode;
      }
      // add as last child
      child->nextNode=r;
      r->prevNode=child;
    } else {
      childNode=r;
    }
    UpdateFirstLastBin();
    r=binning;
  }
  return r;
}

Bool_t TUnfoldBinningV17::AddAxis
(const char *name,Int_t nBin,Double_t xMin,Double_t xMax,
 Bool_t hasUnderflow,Bool_t hasOverflow)
{
  // add an axis with equidistant bins to the distribution
  //    name: name of the axis
  //    nBin: number of bins
  //    xMin: lower edge of the first bin
  //    xMax: upper edge of the last bin
  //    hasUnderflow: decide whether the axis has an underflow bin
  //    hasOverflow: decide whether the axis has an overflow bin
  // return: true if the axis has been added
  Bool_t r=kFALSE;
   if(nBin<=0) {
      Fatal("AddAxis","number of bins %d is not positive",
            nBin);
   } else if((!TMath::Finite(xMin))||(!TMath::Finite(xMax))||
      (xMin>=xMax)) {
      Fatal("AddAxis","xmin=%f required to be smaller than xmax=%f",
            xMin,xMax);
   } else {
     Double_t *binBorders=new Double_t[nBin+1];
     Double_t x=xMin;
     Double_t dx=(xMax-xMin)/nBin;
     for(Int_t i=0;i<=nBin;i++) {
       binBorders[i]=x+i*dx;
     }
     r=AddAxis(name,nBin,binBorders,hasUnderflow,hasOverflow);
     delete [] binBorders;
   }
   return r;
}

Bool_t TUnfoldBinningV17::AddAxis
(const TAxis &axis,Bool_t hasUnderflow,Bool_t hasOverflow)
{
  // add an axis to the distribution
  //    axis: the axis
  //    hasUnderflow: decide whether the underflow bin should be included
  //    hasOverflow: decide whether the overflow bin should be included
  // return: true if the axis has been added
  //
  // Note: axis labels are not imported
  Int_t nBin=axis.GetNbins();
  Double_t *binBorders=new Double_t[nBin+1];
  for(Int_t i=0;i<nBin;i++) {
    binBorders[i]=axis.GetBinLowEdge(i+1);
  }
  binBorders[nBin]=axis.GetBinUpEdge(nBin);
  Bool_t r=AddAxis(axis.GetTitle(),nBin,binBorders,hasUnderflow,hasOverflow);
  delete [] binBorders;
  return r;
}

Bool_t TUnfoldBinningV17::AddAxis
(const char *name,Int_t nBin,const Double_t *binBorders,
 Bool_t hasUnderflow,Bool_t hasOverflow)
{
  // add an axis with the specified bin borders to the distribution
  //    name: name of the axis
  //    nBin: number of bins
  //    binBorders: array of bin borders, with nBin+1 elements
  //    hasUnderflow: decide whether the axis has an underflow bin
  //    hasOverflow: decide whether the axis has an overflow bin
  Bool_t r=kFALSE;
  if(HasUnconnectedBins()) {
    Fatal("AddAxis","node already has %d bins without axis",
	  GetDistributionNumberOfBins());
  } else if(nBin<=0) {
    Fatal("AddAxis","number of bins %d is not positive",
	  nBin);
  } else {
    TVectorD *bins=new TVectorD(nBin+1);
    r=kTRUE;
    for(Int_t i=0;i<=nBin;i++) {
      (*bins)(i)=binBorders[i];
      if(!TMath::Finite((*bins)(i))) {
	Fatal("AddAxis","bin border %d is not finite",i);
	r=kFALSE;
      } else if((i>0)&&((*bins)(i)<=(*bins)(i-1))) {
	Fatal("AddAxis","bins not in order x[%d]=%f <= %f=x[%d]",
	      i,(*bins)(i),(*bins)(i-1),i-1);
	r=kFALSE;
      }
    }
    if(r) {
      Int_t axis=fAxisList->GetEntriesFast();
      Int_t bitMask=1<<axis;
      Int_t nBinUO=nBin;
      if(hasUnderflow) {
	fHasUnderflow |= bitMask;
	nBinUO++;
      } else {
	fHasUnderflow &= ~bitMask;
      }
      if(hasOverflow) {
	fHasOverflow |= bitMask;
	nBinUO++;
      } else {
	fHasOverflow &= ~bitMask;
      }
      fAxisList->AddLast(bins);
      fAxisLabelList->AddLast(new TObjString(name));
      if(!fDistributionSize) fDistributionSize=1;
      fDistributionSize *= nBinUO;
      UpdateFirstLastBin();
    }
  }
  return r;
}

void TUnfoldBinningV17::PrintStream(ostream &out,Int_t indent,int debug)
   const {
  // print some information about this binning tree
  //    out: stream to write to
  //    indent: initial indentation (sub-trees have indent+1)
   for(Int_t i=0;i<indent;i++) out<<"  ";
   out<<"TUnfoldBinning \""<<GetName()<<"\" has ";
   Int_t nBin=GetEndBin()-GetStartBin();
   if(nBin==1) {
      out<<"1 bin";
   } else {
      out<<nBin<<" bins";
   }
   out<<" ["
      <<GetStartBin()<<","<<GetEndBin()<<"] nTH1x="
      <<GetTH1xNumberOfBins()
      <<"\n";
   if(GetDistributionNumberOfBins()) {
      for(Int_t i=0;i<indent;i++) out<<"  ";
      out<<" distribution: "<<GetDistributionNumberOfBins()<<" bins\n";
      if(fAxisList->GetEntriesFast()) {
         /* for(Int_t i=0;i<indent;i++) out<<"  ";
            out<<" axes:\n"; */
          for(Int_t axis=0;axis<GetDistributionDimension();axis++) {
             for(Int_t i=0;i<indent;i++) out<<"  ";
             out<<"  \""
                 <<GetDistributionAxisLabel(axis)
                 <<"\" nbin="<<GetDistributionBinning(axis)->GetNrows()-1;
             if(HasUnderflow(axis)) out<<" plus underflow";
             if(HasOverflow(axis)) out<<" plus overflow";
             out<<"\n";
          }
      } else {
         for(Int_t i=0;i<indent;i++) out<<"  ";
         out<<" no axis\n";
         for(Int_t i=0;i<indent;i++) out<<"  ";
         out<<" names: ";
         for(Int_t ibin=0;(ibin<GetDistributionNumberOfBins())&&
                (ibin<fAxisLabelList->GetEntriesFast());ibin++) {
            if(ibin) out<<";";
            if(GetDistributionAxisLabel(ibin)) {
               out<<GetDistributionAxisLabel(ibin);
            }
         }
         out<<"\n";
      }
      if(debug>0) {
         // print all bins with full name, size, status, user factor
         for(int iBin=GetStartBin();iBin<GetEndBin();iBin++) {
            for(Int_t i=0;i<indent;i++) out<<"  ";
            out<<GetBinName(iBin)
               <<" size="<<GetBinSize(iBin)
               <<" factor="<<GetBinFactor(iBin);
            out<<"\n";
         }
      }
   }
   TUnfoldBinningV17 const *child=GetChildNode();
   if(child) {
      while(child) {
         child->PrintStream(out,indent+1,debug);
         child=child->GetNextNode();
      }
   }
}

void TUnfoldBinningV17::SetBinFactor
(Double_t normalisation,TObject *binfactor) {
   fBinFactorConstant=normalisation;
   if(fBinFactorFunction) {
      if(!dynamic_cast<TF1 *>(fBinFactorFunction)) 
         delete fBinFactorFunction;
   }
   fBinFactorFunction=binfactor;
}

void TUnfoldBinningV17::SetBinFactorFunction
(Double_t normalisation,TF1 *userFunc) {
   SetBinFactor(normalisation,userFunc);
}

/********************* Navigation **********************/

TUnfoldBinningV17 const *TUnfoldBinningV17::FindNode(char const *name) const
{
   // parse the tree and return a node with the given name
   //   name: the name of the node to find
   TUnfoldBinningV17 const *r=0;
   if((!name)||(!TString(GetName()).CompareTo(name))) {
      r=this;
   }
   for(TUnfoldBinningV17 const *child=GetChildNode();
       (!r) && child;child=child->GetNextNode()) {
      r=child->FindNode(name);
   }
   return r;
}

TUnfoldBinningV17 *TUnfoldBinningV17::GetRootNode(void)
{
   // return root node
   TUnfoldBinningV17 *node=this;
   while(node->GetParentNode()) node=node->parentNode;
   return node;
}

TUnfoldBinningV17 const *TUnfoldBinningV17::GetRootNode(void) const
{
   // return root node
   TUnfoldBinningV17 const *node=this;
   while(node->GetParentNode()) node=node->GetParentNode();
   return node;
}

/********************* Create THxx histograms **********/

TString TUnfoldBinningV17::BuildHistogramTitle
(const char *histogramName,const char *histogramTitle,Int_t const *axisList)
   const
{
   // build a title
   // input:
   //   histogramTitle : if this is non-zero, use that title
   //  otherwise:
   //   title=histogramName[;x[;y[;z]]]
   //   ?Axis : -2 stop adding text to the title
   //           -1 add name of this node
   //          >=0 use name of the corresponding axis
   TString r;
   if(histogramTitle) {
      r=histogramTitle;
   } else {
      r=histogramName;
      Int_t iEnd;
      for(iEnd=2;iEnd>0;iEnd--) {
         if(axisList[iEnd]>=0) break;
      }
      for(Int_t i=0;i<=iEnd;i++) {
         r += ";";
         if(axisList[i]<0) {
            r += GetName();
         } else {
            r += GetNonemptyNode()->GetDistributionAxisLabel(axisList[i]);
         }
      }
   }
   return r;
}

TString TUnfoldBinningV17::BuildHistogramTitle2D
(const char *histogramName,const char *histogramTitle,
 Int_t xAxis,const TUnfoldBinningV17 *yAxisBinning,Int_t yAxis) const
{
   // build a title
   // input:
   //   histogramTitle : if this is non-zero, use that title
   //  otherwise:
   //   title=histogramName;x;y
   //   xAxis : -1 no title for this axis
   //          >=0 use name of the corresponding axis
   TString r;
   if(histogramTitle) {
      r=histogramTitle;
   } else {
      r=histogramName;
      r += ";";
      if(xAxis==-1) {
         r += GetName();
      } else if(xAxis>=0) {
         r += GetNonemptyNode()->GetDistributionAxisLabel(xAxis);
      }
      r+= ";";
      if(yAxis==-1) {
         r += yAxisBinning->GetName();
      } else if(yAxis>=0) {
         r += yAxisBinning->GetNonemptyNode()->GetDistributionAxisLabel(yAxis);
      }

   }
   return r;
}

Int_t TUnfoldBinningV17::GetTH1xNumberOfBins
(Bool_t originalAxisBinning,const char *axisSteering) const
{
  // return the number of histogram bins required when storing
  // this binning in a one-dimensional histogram
  //     originalAxisBinning : try to preserve the axis binning,
  //                           if this binning is 1-dimensional then the
  //                           underflow/overflow are stored in the
  //                           corresponding bins of the TH1 and 
  //   axisSteering:
  //       "pattern1;pattern2;...;patternN"
  //       patternI = axis[mode]
  //       axis = name or *
  //       mode = C|U|O 
  //        C: collapse axis into one bin
  //        U: discard underflow bin
  //        O: discard overflow bin
  //  return: number of bins of the TH1 (underflow/overflow are not counted)
   Int_t axisBins[3],axisList[3];
   GetTHxxBinning(originalAxisBinning ? 1 : 0,axisBins,axisList,
                  axisSteering);
   return axisBins[0];
}

TH1 *TUnfoldBinningV17::CreateHistogram
(const char *histogramName,Bool_t originalAxisBinning,Int_t **binMap,
 const char *histogramTitle,const char *axisSteering) const
{
  // create a THxx histogram capable to hold the bins of this binning
  // scheme and its children
  //  input:
  //     histogramName: name of the histogram which is created
  //     originalAxisBinning : try to preserve the axis binning
  //                           if this parameter is true, the resulting
  //                           histogram has bin widths and histogram
  //                           dimension (TH1D, TH2D, TH3D)
  //                           in parallel to the binning scheme
  //                           (if possible)
  //     binMap : mapping of global bins to histogram bins 
  //               see method CreateBinMap()
  //               if(binMap==0), no binMap is created
  //     histogramTitle: if this is non-zero, it is taken as histogram title
  //                     otherwise, the title is created automatically
  //     axisSteering:
  //       "pattern1;pattern2;...;patternN"
  //       patternI = axis[mode]
  //       axis = name or *
  //       mode = C|U|O 
  //        C: collapse axis into one bin
  //        U: discard underflow bin
  //        O: discard overflow bin
  // returns: a new histogram (TH1D, TH2D or TH3D)
   Int_t nBin[3],axisList[3];
   Int_t nDim=GetTHxxBinning(originalAxisBinning ? 3 : 0,nBin,axisList,
                             axisSteering);
   const TUnfoldBinningV17 *neNode=GetNonemptyNode();
   TString title=BuildHistogramTitle(histogramName,histogramTitle,axisList);
   TH1 *r=0;
   if(nDim>0) {
      const TVectorD *axisBinsX=
         neNode->GetDistributionBinning(axisList[0]);
      if(nDim>1) {
         const TVectorD *axisBinsY=
            neNode->GetDistributionBinning(axisList[1]);
         if(nDim>2) {
            const TVectorD *axisBinsZ=
               neNode->GetDistributionBinning(axisList[2]);
            r=new TH3D(histogramName,title,
                       nBin[0],axisBinsX->GetMatrixArray(),
                       nBin[1],axisBinsY->GetMatrixArray(),
                       nBin[2],axisBinsZ->GetMatrixArray());
         } else {
            r=new TH2D(histogramName,title,
                       nBin[0],axisBinsX->GetMatrixArray(),
                       nBin[1],axisBinsY->GetMatrixArray());
         }
      } else {
         r=new TH1D(histogramName,title,nBin[0],axisBinsX->GetMatrixArray());
      }
   } else {
      if(originalAxisBinning) {
         Warning("CreateHistogram",
		 "Original binning can not be represented as THxx");
      }
      r=new TH1D(histogramName,title,nBin[0],0.5,nBin[0]+0.5);
      nDim=0;
   }
   if(binMap) {
      *binMap=CreateBinMap(r,nDim,axisList,axisSteering);
   }
   return r;
}

TH2D *TUnfoldBinningV17::CreateErrorMatrixHistogram
(const char *histogramName,Bool_t originalAxisBinning,Int_t **binMap,
 const char *histogramTitle,const char *axisSteering) const
{
  // create a TH2D histogram capable to hold an error matrix
  // coresponding to the bins of this binning scheme and its children
  //  input:
  //     histogramName: name of the histogram which is created
  //     originalAxisBinning : try to preserve the axis binning
  //                           if this parameter is true, the resulting
  //                           histogram has bin widths
  //                           in parallel to this binning scheme
  //                           (if possible)
  //     binMap : mapping of global bins to histogram bins 
  //               see method CreateBinMap()
  //               if(binMap==0), no binMap is created
  //     histogramTitle: if this is non-zero, it is taken as histogram title
  //                     otherwise, the title is created automatically
  //     axisSteering:
  //       "pattern1;pattern2;...;patternN"
  //       patternI = axis[mode]
  //       axis = name or *
  //       mode = C|U|O 
  //        C: collapse axis into one bin
  //        U: discard underflow bin
  //        O: discard overflow bin
  // returns: a new TH2D

   Int_t nBin[3],axisList[3];
   Int_t nDim=GetTHxxBinning(originalAxisBinning ? 1 : 0,nBin,axisList,
                             axisSteering);
   TString title=BuildHistogramTitle(histogramName,histogramTitle,axisList);
   TH2D *r=0;
   if(nDim==1) {
      const TVectorD *axisBinsX=(TVectorD const *)
         GetNonemptyNode()->fAxisList->At(axisList[0]);
      r=new TH2D(histogramName,title,nBin[0],axisBinsX->GetMatrixArray(),
                 nBin[0],axisBinsX->GetMatrixArray());
   } else {
      if(originalAxisBinning) {
         Info("CreateErrorMatrixHistogram",
              "Original binning can not be represented on one axis");
      }
      r=new TH2D(histogramName,title,nBin[0],0.5,nBin[0]+0.5,
                 nBin[0],0.5,nBin[0]+0.5);
      nDim=0;
   }
   if(binMap) {
      *binMap=CreateBinMap(r,nDim,axisList,axisSteering);
   }
   return r;
}

TH2D *TUnfoldBinningV17::CreateHistogramOfMigrations
(TUnfoldBinningV17 const *xAxis,TUnfoldBinningV17 const *yAxis,
 char const *histogramName,Bool_t originalXAxisBinning,
 Bool_t originalYAxisBinning,char const *histogramTitle)
{
  // create a TH2D histogram capable to hold the bins of the two
  // input binning schemes on the x and y axes, respectively
  //  input:
  //     histogramName: name of the histogram which is created
  //     xAxis: binning scheme for the x axis
  //     yAxis: binning scheme for the y axis
  //     originalXAxisBinning: preserve x-axis bin widths if possible
  //     originalXAxisBinning: preserve y-axis bin widths if possible
  //     histogramTitle: if this is non-zero, it is taken as histogram title
  //                     otherwise, the title is created automatically
  // returns: a new TH2D

   Int_t nBinX[3],axisListX[3];
   Int_t nDimX=
      xAxis->GetTHxxBinning(originalXAxisBinning ? 1 : 0,nBinX,axisListX,0);
   Int_t nBinY[3],axisListY[3];
   Int_t nDimY=
      yAxis->GetTHxxBinning(originalYAxisBinning ? 1 : 0,nBinY,axisListY,0);
   TString title=xAxis->BuildHistogramTitle2D
      (histogramName,histogramTitle,axisListX[0],yAxis,axisListY[0]);
   if(nDimX==1) {
      const TVectorD *axisBinsX=(TVectorD const *)
         xAxis->fAxisList->At(axisListX[0]);
      if(nDimY==1) {
         const TVectorD *axisBinsY=(TVectorD const *)
            yAxis->fAxisList->At(axisListY[0]);
         return new TH2D(histogramName,title,
                         nBinX[0],axisBinsX->GetMatrixArray(),
                         nBinY[0],axisBinsY->GetMatrixArray());
      } else {
         return new TH2D(histogramName,title,
                         nBinX[0],axisBinsX->GetMatrixArray(),
                         nBinY[0],0.5,0.5+nBinY[0]);
      }
   } else {
      if(nDimY==1) {
         const TVectorD *axisBinsY=(TVectorD const *)
            yAxis->fAxisList->At(axisListY[0]);
         return new TH2D(histogramName,title,
                         nBinX[0],0.5,0.5+nBinX[0],
                         nBinY[0],axisBinsY->GetMatrixArray());
      } else {
         return new TH2D(histogramName,title,
                         nBinX[0],0.5,0.5+nBinX[0],
                         nBinY[0],0.5,0.5+nBinY[0]);
      }
   }
}

Int_t TUnfoldBinningV17::GetTHxxBinning
(Int_t maxDim,Int_t *axisBins,Int_t *axisList,
 const char *axisSteering) const
{
   // calculate properties of a THxx histogram to store this binning
   // input:
   //   maxDim : maximum dimension of the THxx (0 or 1..3)
   //              maxDim==0 is used to indicate that the histogram should be
   //              1-dimensional with all bins mapped on one axis,
   //               bin centers equal to bin numbers
   //   axisSteering:
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O 
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   // output;
   //   axisBins[0..2] : number of bins on the THxx axes [0]:x [1]:y [2]:z
   //   axisList[0..2] : TUnfoldBinning axis number corresponding to TH1 axis
   //                        [0]:x [1]:y [2]:z
   // return value :  1-3: dimension of THxx
   //                 0 : use 1-dim THxx, binning structure is not preserved
   for(Int_t i=0;i<3;i++) {
      axisBins[i]=0;
      axisList[i]=-1;
   }
   const TUnfoldBinningV17 *theNode=GetNonemptyNode();
   if(theNode) {
      Int_t r=theNode->GetTHxxBinningSingleNode
         (maxDim,axisBins,axisList,axisSteering);
      return r;
   } else {
      axisBins[0]=GetTHxxBinsRecursive(axisSteering);
      return 0;
   }
}

const TUnfoldBinningV17 *TUnfoldBinningV17::GetNonemptyNode(void) const
{
   // get the node which has non-empty distributions
   // if there is none or if there are many, return zero
   const TUnfoldBinningV17 *r=GetDistributionNumberOfBins()>0 ? this : 0;
   for(TUnfoldBinningV17 const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      const TUnfoldBinningV17 *c=child->GetNonemptyNode();
      if(!r) {
         // new candidate found
         r=c;
      } else {
         if(c) {
            // multiple nodes found
            r=0;
            break;
         }
      }
   }
   return r;
}

Int_t TUnfoldBinningV17::GetTHxxBinningSingleNode
(Int_t maxDim,Int_t *axisBins,Int_t *axisList,const char *axisSteering) const
{
   // get the properties of a histogram capable to hold the distribution
   // of this node
   // input:
   //   maxDim : maximum dimension of the THxx (0 or 1..3)
   //              maxDim==0 is used to indicate that the histogram should be
   //              1-dimensional with all bins mapped on one axis
   //   axisSteering :
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   //   input/output:
   //    axisBins[0..2] : cumulated number of bins on the THxx axes
   //                         [0]:x [1]:y [2]:z
   //    axisList[0..2] : TUnfoldBinning axis number corresponding to TH1 axis
   //                         [0]:x [1]:y [2]:z
   // return value :  1-3: dimension of THxx
   //                 0 : use 1-dim THxx, binning structure is not preserved


   // decode axisSteering
   //   isOptionGiven[0] ('C'): bit vector which axes to collapse
   //   isOptionGiven[1] ('U'): bit vector to discard underflow bins
   //   isOptionGiven[2] ('O'): bit vector to discard overflow bins
   Int_t isOptionGiven[3];
   DecodeAxisSteering(axisSteering,"CUO",isOptionGiven);
   // count number of axes after projecting
   Int_t numDimension=GetDistributionDimension();
   Int_t r=0;
   for(Int_t i=0;i<numDimension;i++) {
      if(isOptionGiven[0] & (1<<i)) continue;
      r++;
   }
   if((r>0)&&(r<=maxDim)) {
      // 0<r<=maxDim
      //
      // -> preserve the original binning
      //    axisList[] and axisBins[] are overwritten
      r=0;
      for(Int_t i=0;i<numDimension;i++) {
         if(isOptionGiven[0] & (1<<i)) continue;
         axisList[r]=i;
         axisBins[r]=GetDistributionBinning(i)->GetNrows()-1;
         r++;
      }
   } else {
      // map everything on one axis
      //  axisBins[0] is the number of bins
      if(HasUnconnectedBins() || (GetDistributionNumberOfBins()<=0)) {
         axisBins[0] = GetDistributionNumberOfBins();
      } else {
         Int_t nBin=1;
         for(Int_t i=0;i<numDimension;i++) {
            Int_t mask=(1<<i);
            if(isOptionGiven[0] & mask) continue;
            Int_t nBinI=GetDistributionBinning(i)->GetNrows()-1;
            if((fHasUnderflow & mask)&& !(isOptionGiven[1] & mask)) nBinI++;
            if((fHasOverflow & mask)&& !(isOptionGiven[2] & mask)) nBinI++;
            nBin *= nBinI;
         }
         axisBins[0] = nBin;
      }
      r=0;
   }
   return r;
}

Int_t TUnfoldBinningV17::GetTHxxBinsRecursive(const char *axisSteering) const
{
   // input:
   //   axisSteering :
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O 
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   // output:
   //   binMap[] : map global bin numbers to histogram bins
   // return value :  number of bins

   Int_t r=0;
   for(TUnfoldBinningV17 const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      r +=child->GetTHxxBinsRecursive(axisSteering);
   }
   // here: process distribution of this node
   Int_t axisBins[3],axisList[3];
   GetTHxxBinningSingleNode(0,axisBins,axisList,axisSteering);
   r += axisBins[0];
   return r;
}

Int_t *TUnfoldBinningV17::CreateEmptyBinMap(void) const {
   // create empty bin map which can be manipulated by
   //  MapGlobalBin()
   Int_t nMax=GetRootNode()->GetEndBin()+1;
   Int_t *r=new Int_t[nMax];
   for(Int_t i=0;i<nMax;i++) {
         r[i]=-1;
   }
   return r;
}

void TUnfoldBinningV17::SetBinMapEntry
(Int_t *binMap,Int_t globalBin,Int_t destBin) const {
   // set one entry in the bin map
   //   binMap : bin map, to be used with TUnfoldSys::GetOutput() etc
   //   globalBin : source bin, global bin number in this binning scheme
   //   destBin : destination bin.
   //        global bins with the same destination bin will be added
   //        if GetOutput is called with a bin map set up this way
   Int_t nMax=GetRootNode()->GetEndBin()+1;
   if((globalBin<0)||(globalBin>=nMax)) {
      Error("SetBinMapEntry","global bin number %d outside range (max=%d)",
            globalBin,nMax);
   } else {
      binMap[globalBin]=destBin;
   }
}

Int_t TUnfoldBinningV17::FillBinMap1D
(Int_t *binMap,const char *axisSteering,TH1 *destHist,
 Int_t firstBinX) const {
   // map all global bins referenced by this node to the one-dimensional
   // histogram destHist, starting with bin firstBinX
   //   binMap : the bin map to be filled
   //  axisSteering:
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O 
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   //  destHist: destination histogram
   //  firstBinX: first bin in destination histogram which should be filled
   // returns: first bin number of destHist which has not been mapped
   if(destHist->GetDimension()!=1) {
      Error("FillBinMap1D","histogram %s is not 1D",
            (char *)destHist->GetName());
   }
   Int_t r=firstBinX;
   Int_t axisBins[3],axisList[3];
   Int_t nDim=GetTHxxBinningSingleNode(3,axisBins,axisList,axisSteering);
   if((nDim==1)|| !GetDistributionDimension()) {
      r+=FillBinMapSingleNode(0,r,0,0,axisSteering,binMap);
   } else {
      Error("FillBinMap1D","distribution %s with steering=%s is not 1D",
            (char *)GetName(),axisSteering);
   }
   for(TUnfoldBinningV17 const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      r =child->FillBinMap1D(binMap,axisSteering,destHist,r);
   }
   return r;
}

Int_t *TUnfoldBinningV17::CreateBinMap
(const TH1 *hist,Int_t nDim,const Int_t *axisList,const char *axisSteering)
   const
{
   // create mapping from global bin number to a histogram for this node
   // global bins are the bins of the root node binning scheme
   // when projecting them on a TH1 histogram "hRootNode" without special
   // axis steering and without attempting to preserve the axis binnings
   // 
   // The bin map is an array of size hRootNode->GetNbinsX()+2
   // For each bin of the "hRootNode" histogram it holds the target bin in
   // "hist" or the number -1 if the corresponding "hRootNode" bin is not
   // represented in "hist"
   //
   // input
   //   hist : the histogram (to calculate root bin numbers)
   //   nDim : target dimension of the TUnfoldBinning
   //          if(nDim==0) all bins are mapped linearly
   //   axisSteering:
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O 
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   //
   // input used only if nDim>0:
   //    axisList : for each THxx axis give the TUnfoldBinning axis number
   //
   // return value:
   //   an new array which holds the bin mapping 
   //    r[0] : to which THxx bin to map global bin number 0
   //    r[1] : to which THxx bin to map global bin number 1
   //      ...
   //    r[nmax]
   //  where nmax=GetRootNode()->GetEndBin()+1
   Int_t *r=CreateEmptyBinMap();
   Int_t startBin=GetRootNode()->GetStartBin();
   if(nDim>0) {
     const TUnfoldBinning *nonemptyNode=GetNonemptyNode();
     if(nonemptyNode) {
       nonemptyNode->
          FillBinMapSingleNode(hist,startBin,nDim,axisList,axisSteering,r);
     } else {
       Fatal("CreateBinMap","called with nDim=%d but GetNonemptyNode()=0",
	     nDim);
     }
   } else {
     FillBinMapRecursive(startBin,axisSteering,r);
   }
   return r;
}

Int_t TUnfoldBinningV17::FillBinMapRecursive
(Int_t startBin,const char *axisSteering,Int_t *binMap) const
{
   // fill bin map recursively
   // input
   //   startBin : first histogram bin
   //   axisSteering : specify how to project the axes
   //   binMap : the bin mapping which is to be filled
   //
   // the positions
   //       binMap[GetStartBin()]
   //         ..
   //       binMap[GetEndBin()-1]
   // are filled
   //
   Int_t nbin=0;
   nbin = FillBinMapSingleNode(0,startBin,0,0,axisSteering,binMap);
   for(TUnfoldBinningV17 const *child=GetChildNode();child;
       child=child->GetNextNode()) {
      nbin += child->FillBinMapRecursive(startBin+nbin,axisSteering,binMap);
   }
   return nbin;
}

Int_t TUnfoldBinningV17::FillBinMapSingleNode
(const TH1 *hist,Int_t startBin,Int_t nDim,const Int_t *axisList,
 const char *axisSteering,Int_t *binMap) const
{
  // fill bin map for a single node
  //  input:
  //     hist: the histogram representing this node (used if nDim>0)
  //     startBin: start bin in the bin map
  //     nDim:
  //        0: bins are mapped in linear order, ignore hist and axisList
  //        nDim=hist->GetDimension():
  //           bins are mapped to "hist" bin numbers
  //           the corresponding TUnfoldBinning axes are taken from axisList[]
  //        nDim=1 and hist->GetDimension()>1:
  //           bins are mapped to th x-axis of "hist"
  //           the corresponding TUnfoldBinning axis is taken from axisList[0]
  //     axisList[]:
  //           TUnfoldBinning axis numbers corresponding to the
  //           x[0], y[1], z[2] axes of "hist"
  //     axisSteering:
  //       "pattern1;pattern2;...;patternN"
  //       patternI = axis[mode]
  //       axis = name or *
  //       mode = C|U|O|0..9
  //        C: collapse axis into one bin
  //        U: discard underflow bin
  //        O: discard overflow bin
  //        0,1,2,3,4,5,6,7,8,9: define which bins
  //                             of the collapsed axis to include
  //                             only works for axes with up to 10 bins
  //     binMap: the bin map to fill
  // return value: 
  //     the number of bins mapped
  //     (only relevant if nDim==0)
   // first, decode axisSteering
   //   isOptionGiven[0] ('C'): bit vector which axes to collapse
   //   isOptionGiven[1] ('U'): bit vector to discard underflow bins
   //   isOptionGiven[2] ('O'): bit vector to discard overflow bins
   Int_t isOptionGiven[3+10];
   DecodeAxisSteering(axisSteering,"CUO0123456789",isOptionGiven);
   Int_t haveSelectedBin=0; 
   for(Int_t i=3;i<3+10;i++) {
      haveSelectedBin |= isOptionGiven[i];
   }

   Int_t axisBins[MAXDIM];
   Int_t dimension=GetDistributionDimension();
   Int_t axisNbin[MAXDIM];
   for(Int_t i=0;i<dimension;i++) {
      const TVectorD *binning=GetDistributionBinning(i);
      axisNbin[i]=binning->GetNrows()-1;
   };
   for(Int_t i=0;i<GetDistributionNumberOfBins();i++) {
      Int_t globalBin=GetStartBin()+i;
      const TUnfoldBinningV17 *dest=ToAxisBins(globalBin,axisBins);
      if(dest!=this) {
         if(!dest) {
            Fatal("FillBinMapSingleNode",
                  "bin %d outside binning scheme",
                  globalBin);
         } else {
            Fatal("FillBinMapSingleNode",
                  "bin %d located in %s %d-%d rather than %s %d=%d",
                  i,(const char *)dest->GetName(),
                  dest->GetStartBin(),dest->GetEndBin(),
                  (const char *)GetName(),GetStartBin(),GetEndBin());
         }
      }
      // check whether this bin has to be skipped
      Bool_t skip=kFALSE;
      for(Int_t axis=0;axis<dimension;axis++) {
         Int_t mask=(1<<axis);
         // underflow/overflow excluded by steering
         if(((axisBins[axis]<0)&&(isOptionGiven[1] & mask))||
            ((axisBins[axis]>=axisNbin[axis])&&(isOptionGiven[2] & mask)))
            skip=kTRUE;
         // only certain bins selected by steering
         if((axisBins[axis]>=0)&&(axisBins[axis]<axisNbin[axis])&&
            (haveSelectedBin & mask)) {
            if(!(isOptionGiven[3+axisBins[axis]] & mask)) skip=kTRUE;
         }
      }
      if(skip) {
         continue;
      }

      if(nDim>0) {
         // get bin number from THxx function(s)
         if(nDim==hist->GetDimension()) {
            Int_t ibin[3];
            ibin[0]=ibin[1]=ibin[2]=0;
            for(Int_t hdim=0;hdim<nDim;hdim++) {
               Int_t axis=axisList[hdim];
               ibin[hdim]=axisBins[axis]+1;
            }
            binMap[globalBin]=hist->GetBin(ibin[0],ibin[1],ibin[2]);
         } else if(nDim==1) {
            // histogram has more dimensions than the binning scheme
            // and the binning scheme has one axis only
            // -> use the first valid axis only
            Error("FillBinMapSingleNode","inconsistent dimensions %d %d",nDim,
                  hist->GetDimension());
            for(Int_t ii=0;ii<hist->GetDimension();ii++) {
               if(axisList[ii]>=0) {
                  binMap[globalBin]=axisBins[axisList[ii]]+1;
                  break;
               }
            }
         } else {
            Fatal("FillBinMapSingleNode","inconsistent dimensions %d %d",nDim,
                  hist->GetDimension());
         }
      } else {
         // order all bins in sequence
         // calculation in parallel to ToGlobalBin()
         // but take care of
         //   startBin,collapseAxis,discardeUnderflow,discardeOverflow
         if(dimension>0) {
            Int_t r=0;
            for(Int_t axis=dimension-1;axis>=0;axis--) {
               Int_t mask=(1<<axis);
               if(isOptionGiven[0] & mask) {
                  // bins on this axis are integrated over
                  continue;
               }
               Int_t iBin=axisBins[axis];
               Int_t nMax=axisNbin[axis];
               if((fHasUnderflow & ~isOptionGiven[1]) & mask) {
                  nMax +=1;
                  iBin +=1;
               }
               if((fHasOverflow & ~isOptionGiven[2]) & mask) {
                  nMax += 1;
               }
               r = r*nMax +iBin;
            }
            binMap[globalBin] = startBin + r;
         } else {
	   binMap[globalBin] = startBin + axisBins[0];
         }
      }
   }
   Int_t nbin;
   if(dimension>0) {
     nbin=1;
     for(Int_t axis=dimension-1;axis>=0;axis--) {
       Int_t mask=(1<<axis);
       if(isOptionGiven[0] & mask) {
	 // bins on this axis are integrated over
	 continue;
       }
       Int_t nMax=axisNbin[axis];
       if((fHasUnderflow & ~isOptionGiven[1]) & mask) {
	 nMax +=1;
       }
       if((fHasOverflow & ~isOptionGiven[2]) & mask) {
	 nMax += 1;
       }
       nbin = nbin*nMax;
     }
   } else {
     nbin=GetDistributionNumberOfBins();
   }
   return nbin;
}

TH1 *TUnfoldBinningV17::ExtractHistogram
(const char *histogramName,const TH1 *globalBins,
 const TH2 *globalBinsEmatrix,Bool_t originalAxisBinning,
 const char *axisSteering) const
{
   // extract a distribution from the given set of global bins
   // input:
   //    histogramName : name of the histogram which ic created
   //    globalBins : histogram with all bins
   //    globalBinsEmatrix : corresponding error matrix
   //                 if this pointer is zero, only diagonal errors
   //                 are considered
   //    originalAxisBinning :  extract  histogram with proper binning
   //                          (if possible)
   //    axisSteering
   //       "pattern1;pattern2;...;patternN"
   //       patternI = axis[mode]
   //       axis = name or *
   //       mode = C|U|O 
   //        C: collapse axis into one bin
   //        U: discard underflow bin
   //        O: discard overflow bin
   Int_t *binMap=0;
   TH1 *r=CreateHistogram(histogramName,originalAxisBinning,&binMap,0,
                          axisSteering);
   if(!r) return 0;
   TUnfoldBinningV17 const *root=GetRootNode();
   Int_t nMax=-1;
   for(Int_t iSrc=root->GetStartBin();iSrc<root->GetEndBin();iSrc++) {
      if(binMap[iSrc]>nMax) nMax=binMap[iSrc];
   }
   if(nMax<0) {
      delete r;
      r=0;
   } else {
      TVectorD eSquared(nMax+1);
      for(Int_t iSrc=root->GetStartBin();iSrc<root->GetEndBin();iSrc++) {
         Int_t iDest=binMap[iSrc];
         if(iDest>=0) {
            Double_t c=r->GetBinContent(iDest);
            r->SetBinContent(iDest,c+globalBins->GetBinContent(iSrc));
            if(!globalBinsEmatrix) {
               eSquared(iDest)+=TMath::Power(globalBins->GetBinError(iSrc),2.);
            } else {
               for(Int_t jSrc=root->GetStartBin();jSrc<root->GetEndBin();
                   jSrc++) {
                  if(binMap[jSrc]==iDest) {
                     eSquared(iDest) += 
                        TMath::Power(globalBins->GetBinError(jSrc),2.);
                  }
               }
            }
         }
      }
      for(Int_t i=0;i<nMax;i++) {
         Double_t e2=eSquared(i);
         if(e2>0.0) {
            r->SetBinError(i,TMath::Sqrt(e2));
         }
      }
   }
   delete binMap;
   return r;
}

/********************* Calculate global bin number ******/

Int_t TUnfoldBinningV17::GetGlobalBinNumber(Double_t x) const
{
  // locate bin on a one-dimensional distribution
  // input
  //    x: coordinate to locate
   if(GetDistributionDimension()!=1) {
      Fatal("GetBinNumber",
            "called with 1 argument for %d dimensional distribution",
            GetDistributionDimension());
   }
   return GetGlobalBinNumber(&x);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber(Double_t x,Double_t y) const
{
  // locate bin on a two-dimensional distribution
  // input
  //    x,y: coordinates to locate
   if(GetDistributionDimension()!=2) {
      Fatal("GetBinNumber",
            "called with 2 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[2];
   xx[0]=x;
   xx[1]=y;
   return GetGlobalBinNumber(xx);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber
(Double_t x,Double_t y,Double_t z) const
{
  // locate bin on a three-dimensional distribution
  // input
  //    x,y,z: coordinates to locate
   if(GetDistributionDimension()!=3) {
      Fatal("GetBinNumber",
            "called with 3 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   return GetGlobalBinNumber(xx);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3) const
{
  // locate bin on a four-dimensional distribution
  // input
  //    x0,x1,x2,x3: coordinates to locate
   if(GetDistributionDimension()!=4) {
      Fatal("GetBinNumber",
            "called with 4 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[4];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   return GetGlobalBinNumber(xx);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3,Double_t x4) const
{
  // locate bin on a five-dimensional distribution
  // input
  //    x0,x1,x2,x3,x4: coordinates to locate
   if(GetDistributionDimension()!=5) {
      Fatal("GetBinNumber",
            "called with 5 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[5];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   xx[4]=x4;
   return GetGlobalBinNumber(xx);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber
(Double_t x0,Double_t x1,Double_t x2,Double_t x3,Double_t x4,Double_t x5) const
{
  // locate bin on a five-dimensional distribution
  // input
  //    x0,x1,x2,x3,x4,x5: coordinates to locate
   if(GetDistributionDimension()!=6) {
      Fatal("GetBinNumber",
            "called with 5 arguments for %d dimensional distribution",
            GetDistributionDimension());
   }
   Double_t xx[6];
   xx[0]=x0;
   xx[1]=x1;
   xx[2]=x2;
   xx[3]=x3;
   xx[4]=x4;
   xx[5]=x5;
   return GetGlobalBinNumber(xx);
}

Int_t TUnfoldBinningV17::GetGlobalBinNumber
(const Double_t *x,Int_t *isBelow,Int_t *isAbove) const
{
   // locate bin on a n-dimensional distribution
   // input
   //    x[]: coordinates to locate
   // output:
   //    isBelow,isAbove: bit vectors,
   //       indicating which cut on which axis failed
   if(!GetDistributionDimension()) {
      Fatal("GetBinNumber",
            "no axes are defined for node %s",
            (char const *)GetName());
   }
   Int_t iAxisBins[MAXDIM];
   for(Int_t dim=0;dim<GetDistributionDimension();dim++) {
      TVectorD const *bins=(TVectorD const *) fAxisList->At(dim);
      Int_t i0=0;
      Int_t i1=bins->GetNrows()-1;
      Int_t iBin= 0;
      if(x[dim]<(*bins)[i0]) {
         iBin += i0-1;
         // underflow
      } else if(x[dim]>=(*bins)[i1]) {
         // overflow
         iBin += i1;
      } else {
         while(i1-i0>1) {
            Int_t i2=(i0+i1)/2;
            if(x[dim]<(*bins)[i2]) {
               i1=i2;
            } else {
               i0=i2;
            }
         }
         iBin += i0;
      }
      iAxisBins[dim]=iBin;
   }
   Int_t r=ToGlobalBin(iAxisBins,isBelow,isAbove);
   if(r<0) r=0;
   return r;
}

/********************* access by global bin number ******/

TString TUnfoldBinningV17::GetBinName(Int_t iBin) const
{
  // get the name of a bin in the given tree
  //    iBin: bin number
   Int_t axisBins[MAXDIM];
   TString r=TString::Format("#%d",iBin);
   TUnfoldBinningV17 const *distribution=ToAxisBins(iBin,axisBins);
   if(distribution) {
      r +=" (";
      r += distribution->GetName();
      Int_t dimension=distribution->GetDistributionDimension();
      if(dimension>0) {
         TString axisString;
         for(Int_t axis=0;axis<dimension;axis++) {
            TString thisAxisString=
               distribution->GetDistributionAxisLabel(axis);
            TVectorD const *bins=distribution->GetDistributionBinning(axis);
            Int_t i=axisBins[axis];
            if(i<0) thisAxisString += "[ufl]";
            else if(i>=bins->GetNrows()-1) thisAxisString += "[ofl]";
            else {
               thisAxisString +=
                  TString::Format("[%.3g,%.3g]",(*bins)[i],(*bins)[i+1]);
            }
            axisString = ":"+thisAxisString+axisString;
         }
         r += axisString;
      } else {
         // extra bins
         Int_t i=axisBins[0];
         if((i>=0)&&(i<distribution->fAxisLabelList->GetEntriesFast())) {
            r += distribution->GetDistributionAxisLabel(i);
         } else {
            r += TString::Format(" %d",i);
         }
      }
      r +=")";
   }
   return r;
}

Double_t TUnfoldBinningV17::GetBinSize(Int_t iBin) const
{
  // get N-dimensional bin size
  // input:
  //    iBin : global bin number
  //    includeUO : include underflow/overflow bins or not
   Int_t axisBins[MAXDIM];
   TUnfoldBinningV17 const *distribution=ToAxisBins(iBin,axisBins);
   Double_t r=0.0;
   if(distribution) {
      if(distribution->GetDistributionDimension()>0) r=1.0;
      for(Int_t axis=0;axis<distribution->GetDistributionDimension();axis++) {
         TVectorD const *bins=distribution->GetDistributionBinning(axis);
         Int_t pos=axisBins[axis];
         if(pos<0) {
	    r *= distribution->GetDistributionUnderflowBinWidth(axis);
         } else if(pos>=bins->GetNrows()-1) {
	    r *= distribution->GetDistributionOverflowBinWidth(axis);
         } else {
            r *= (*bins)(pos+1)-(*bins)(pos);
         }
         if(r<=0.) break;
      }
   }
   return r;
}

Bool_t TUnfoldBinningV17::IsBinFactorGlobal(void) const {
   return fBinFactorFunction ? kFALSE : kTRUE;
}

Double_t TUnfoldBinningV17::GetGlobalFactor(void) const {
   return fBinFactorConstant;
}

Double_t TUnfoldBinningV17::GetBinFactor(Int_t iBin) const
{
  // return user factor for a bin
  //    iBin : global bin number
   Int_t axisBins[MAXDIM];
   TUnfoldBinningV17 const *distribution=ToAxisBins(iBin,axisBins);   
   Double_t r=distribution->fBinFactorConstant;
   if((r!=0.0) && distribution->fBinFactorFunction) {
      TF1 *function=dynamic_cast<TF1 *>(distribution->fBinFactorFunction);
      if(function) {
         Double_t x[MAXDIM];
         Int_t dimension=distribution->GetDistributionDimension();
         if(dimension>0) {
            for(Int_t  axis=0;axis<dimension;axis++) {
               x[axis]=distribution->GetDistributionBinCenter
                  (axis,axisBins[axis]);
            }
            r *= function->EvalPar(x,function->GetParameters());
         } else {
            x[0]=axisBins[0];
            r *= function->Eval(x[0]);
         }
      } else {
         TVectorD *vect=dynamic_cast<TVectorD *>
            (distribution->fBinFactorFunction);
         if(vect) {
            r=(*vect)[iBin-GetStartBin()];
         } else {
            Error("GetBinFactor",
                  "internal error: user function is neiter TF1 or TVectorD");
         }
      }
   }
   return r;
}

Int_t TUnfoldBinningV17::GetBinNeighbours
(Int_t bin,Int_t axis,Int_t *prev,Double_t *distPrev,
 Int_t *next,Double_t *distNext,Bool_t isPeriodic) const
{
   // get neighbour bins along the specified axis
   // input:
   //    bin,axis : bin number and axis number
   //    isPeriodic : the first bin is a neighbour of the last bin
   // output:
   //    prev: bin number of previous bin (or -1 if not existing)
   //    distPrev: distance to previous bin
   //    next: bin number of next bin (or -1 if not existing)
   //    distNext: distance to next bin
   // return code:
   //    0: everything is fine
   //   >0: isPeriodic option was ignored
   //   +1: isPeriodic option was specified with underflow bin
   //   +2: isPeriodic option was specified with overflow bin

   Int_t axisBins[MAXDIM];
   TUnfoldBinningV17 const *distribution=ToAxisBins(bin,axisBins);
   Int_t dimension=distribution->GetDistributionDimension();
   *prev=-1;
   *next=-1;
   *distPrev=0.;
   *distNext=0.;
   Int_t r=0;
   if((axis>=0)&&(axis<dimension)) {
     //TVectorD const *bins=distribution->GetDistributionBinning(axis);
      //Int_t nBin=bins->GetNrows()-1;
      Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
      Int_t centerBin= axisBins[axis];
      axisBins[axis] =centerBin-1;
      if(isPeriodic) {
         if(HasUnderflow(axis)) {
            r +=1;
         } else if((axisBins[axis]<0)&&(nMax>=3)) {
            axisBins[axis]=nMax-1;
         }
      }
      *prev=ToGlobalBin(axisBins);
      if(*prev>=0) {
	*distPrev=distribution->GetDistributionBinCenter(axis,axisBins[axis])-
	  distribution->GetDistributionBinCenter(axis,centerBin);
      }
      axisBins[axis] =centerBin+1;
      if(isPeriodic) {
         if(HasOverflow(axis)) {
            r +=2;
         } else if((axisBins[axis]==nMax)&&(nMax>=3)) {
            axisBins[axis]=0;
         }
      }
      *next=ToGlobalBin(axisBins);
      if(*next>=0) {
	*distNext=distribution->GetDistributionBinCenter(axis,axisBins[axis])-
	  distribution->GetDistributionBinCenter(axis,centerBin);
      }
   }
   return r;
}

void TUnfoldBinningV17::GetBinUnderflowOverflowStatus
(Int_t iBin,Int_t *uStatus,Int_t *oStatus) const
{
  // return bit map indicating underflow and overflow status
  //  iBin: global bin number
  // output
  //  uStatus: bin map indicating whether the bin is in underflow
  //  oStatus: bin map indicating whether the bin is in overflow
  Int_t axisBins[MAXDIM];
  TUnfoldBinningV17 const *distribution=ToAxisBins(iBin,axisBins);
  Int_t dimension=distribution->GetDistributionDimension();
  *uStatus=0;
  *oStatus=0;
  for(Int_t axis=0;axis<dimension;axis++) {
    TVectorD const *bins=distribution->GetDistributionBinning(axis);
    Int_t nBin=bins->GetNrows()-1;
    if(axisBins[axis]<0) *uStatus |= (1<<axis);
    if(axisBins[axis]>=nBin) *oStatus |= (1<<axis);
  }
}

Bool_t TUnfoldBinningV17::HasUnconnectedBins(void) const
{
   // check whether there are bins but no axis
   return (!GetDistributionDimension())&&(GetDistributionNumberOfBins()>0);
}

const TObjString *TUnfoldBinningV17::GetUnconnectedBinName(Int_t bin) const {
   TObjString *r=0;
   if(HasUnconnectedBins()) {
      if(bin<fAxisLabelList->GetEntriesFast()) {
         r=((TObjString * const)fAxisLabelList->At(bin));
      }
   }
   return r;
}


Double_t TUnfoldBinningV17::GetDistributionAverageBinSize
(Int_t axis,Bool_t includeUnderflow,Bool_t includeOverflow) const
{
  // get average bin size of the specified axis
  //   axis : axis number
  //   includeUnderflow : include underflow bin
  //   includeOverflow : include overflow bin
   Double_t r=0.0;
   if((axis>=0)&&(axis<GetDistributionDimension())) {
      TVectorD const *bins=GetDistributionBinning(axis);
      Double_t d=(*bins)[bins->GetNrows()-1]-(*bins)[0];
      Double_t nBins=bins->GetNrows()-1;
      if(includeUnderflow && HasUnderflow(axis)) {
         Double_t w=GetDistributionUnderflowBinWidth(axis);
         if(w>0) {
            nBins++;
            d += w;
         }
      }
      if(includeOverflow && HasOverflow(axis)) {
         Double_t w=GetDistributionOverflowBinWidth(axis);
         if(w>0.0) {
            nBins++;
            d += w;
         }
      }
      if(nBins>0) {
         r=d/nBins;
      }
   } else {
      Error("GetDistributionAverageBinSize","axis %d does not exist",axis);
   }
   return r;
}

Double_t TUnfoldBinningV17::GetDistributionUnderflowBinWidth(Int_t axis) const
{
   // return width of the underflow bin
   //   axis: axis number
   TVectorD const *bins=GetDistributionBinning(axis);
   return (*bins)[1]-(*bins)[0];
}

Double_t TUnfoldBinningV17::GetDistributionOverflowBinWidth(Int_t axis) const
{
   // return width of the underflow bin
   //   axis: axis number
   TVectorD const *bins=GetDistributionBinning(axis);
   return (*bins)[bins->GetNrows()-1]-(*bins)[bins->GetNrows()-2];
}

Double_t TUnfoldBinningV17::GetDistributionBinCenter
(Int_t axis,Int_t bin) const
{
   // position of the bin center
   // input
   //   axis : axis number
   //   bin : bin number on the axis
   TVectorD const *bins=GetDistributionBinning(axis);
   Double_t r=0.0;
   if(bin<0) {
      // underflow bin
      r=(*bins)[0]-0.5*GetDistributionUnderflowBinWidth(axis);
   } else if(bin>=bins->GetNrows()-1) {
      // overflow bin
      r=(*bins)[bins->GetNrows()-1]+0.5*GetDistributionOverflowBinWidth(axis);
   } else {
      r=0.5*((*bins)[bin+1]+(*bins)[bin]);
   }
   return r;
}

Int_t TUnfoldBinningV17::ToGlobalBin
(Int_t const *axisBins,Int_t *isBelow,Int_t *isAbove) const
{
  // get global bin number, given axis bin numbers
  //   axisBins[]: bin numbers on each axis
  // return: global bin nmber or -1 if not inside distribution
   //  isBelow: bit vector indicating axes below lowest bin
   //  isAbove: bit vector indicating axes above highest bin
  Int_t dimension=GetDistributionDimension();
  Int_t r=0;
  if(isBelow) *isBelow=0;
  if(isAbove) *isAbove=0;
  if(dimension>0) {
    for(Int_t axis=dimension-1;axis>=0;axis--) {
      Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
      Int_t i=axisBins[axis];
      if(HasUnderflow(axis)) {
	nMax +=1;
	i +=1;
      }
      if(HasOverflow(axis)) nMax +=1;
      if((i>=0)&&(i<nMax)) {
         if(r>=0) r = r*nMax +i;
      } else {
         r=-1;
         if((i<0)&&(isBelow)) *isBelow |= 1<<axis;
         if((i>=nMax)&&(isAbove)) *isAbove |= 1<<axis;
      }
    }
    if(r>=0) {
      r += GetStartBin();
    }
  } else {
    if((axisBins[0]>=0)&&(axisBins[0]<GetDistributionNumberOfBins()))
      r=GetStartBin()+axisBins[0];
    else
       Fatal("ToGlobalBin","bad input %d for dimensionless binning %s %d",
             axisBins[0],(const char *)GetName(),
             GetDistributionNumberOfBins());
  }
  return r;
}

TUnfoldBinningV17 const *TUnfoldBinningV17::ToAxisBins
(Int_t globalBin,Int_t *axisBins) const
{
   // return distribution in which the bin is located
   // and bin numbers on the corresponding axes
   // input
   //   globalBin : global bin number
   // output
   //   if(GetDistributionDimension()>0) 
   //     axisBin[] : bin number on the individual axes
   //                 bin numbers are starting with -1 (underflow)
   //   else
   //     axisBin[0] : bin number, counted from zero
   // return value
   //   the distribution in which the globalBin is located
   //   or 0 if the globalBin is outside the tree
   TUnfoldBinningV17 const *r=0;
   if((globalBin>=GetStartBin())&&(globalBin<GetEndBin())) {
      TUnfoldBinningV17 const *node;
      for(node=GetChildNode();node && !r; node=node->GetNextNode()) {
         r=node->ToAxisBins(globalBin,axisBins);
      }
      if(!r) {
         r=this;
         Int_t i=globalBin-GetStartBin();
         Int_t dimension=GetDistributionDimension();
         if(dimension>0) {
            for(int axis=0;axis<dimension;axis++) {
               Int_t nMax=GetDistributionBinning(axis)->GetNrows()-1;
               axisBins[axis]=0;
               if(HasUnderflow(axis)) {
                  axisBins[axis] =-1;
                  nMax += 1;
               }
               if(HasOverflow(axis)) nMax +=1;
               axisBins[axis] += i % nMax;
               i /= nMax;
            }
         } else {
            axisBins[0]=i;
         }
      }
   }
   return r;
}

void TUnfoldBinningV17::DecodeAxisSteering
(const char *axisSteering,const char *options,Int_t *isOptionGiven) const
{
  // decode axis steering
  // input
  //   axisSteering: the steering to decode
  //   options: the allowed options to extract
  // output
  //   isOptionGiven[] : array of decoded steering options
  //            the dimension of isOptionGiven[] has to match the number of
  //            characters in "option"
  //
  // the axis steering is given in the form
  //   steering1;steering2;...;steeringN
  // all parts of the steering, sepatared by ';' are parsed
  // each part must have the form
  //   axisName[optionlist]
  // where:
  //   axisName : the name of the axis for which the optionlist is relevant
  //               if the name is the character '*', all axes are matched
  //   optionlist : for each option the corresonding character
  //
  // example:
  //   imagine this node has two axes, named "xgen" and "ygen"
  //  then
  //   DecodeAxisSteering("*[U];xgen[C]","CUO",options);
  //  will set the "U" option for all axes and the "C" option for the "xgen"
  //  axis, so:
  //    options[0]=0x1; // option 'C' is true for axis 0 (bit #0 is set)
  //    options[1]=0x3; // option 'U' is true for both axes (bit #0,#1 are set)
  //    options[2]=0x0; // option 'O' is not given (no bits are set)
  Int_t nOpt=TString(options).Length();
  for(Int_t i=0;i<nOpt;i++) isOptionGiven[i]=0;
  if(axisSteering) {
     TObjArray *patterns=TString(axisSteering).Tokenize(";");
     Int_t nPattern=patterns->GetEntries();
     Int_t nAxis=fAxisLabelList->GetEntries();
     for(Int_t i=0;i<nPattern;i++) {
        TString const &pattern=((TObjString * const)patterns->At(i))
           ->GetString();
        Int_t bracketBegin=pattern.Last('[');
        Int_t len=pattern.Length();
        if((bracketBegin>0)&&(pattern[len-1]==']')) {
	  TString axisId=pattern(0,bracketBegin);
	  Int_t mask=0;
	  if((axisId[0]=='*')&&(axisId.Length()==1)) {
	    // turn all bins on 
	    mask=(1<<nAxis)-1;
	  } else {
	    // if axis is there, turn its bit on
	    for(Int_t j=0;j<nAxis;j++) {
	      if(!axisId.CompareTo(GetDistributionAxisLabel(j))) {
		mask|= (1<<j);
	      }
	    }
	  }
	  for(Int_t o=0;o<nOpt;o++) {
	    if(pattern.Last(options[o])>bracketBegin) {
	      isOptionGiven[o] |= mask;
	    }
	  }
        } else {
           Error("DecodeAxisSteering",
                 "steering \"%s\" does not end with [options]",
		 (const char *)pattern);
        }
     }
  }
}

