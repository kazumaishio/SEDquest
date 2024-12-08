
using namespace std;
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "MStatusArray.h"
#include <iostream>
#include <fstream>

int CollareaExtractor(){

  TString inputpath ("/Users/kazuma/Workspace/MAGICana/Analysis/20211111_CrabNebula/Analysis_for_EestPaper/CrabNebula_ST0307_Zd05-35_RFv5/Thsq75_Had90_40bins/Status_flute.root");

  TH1D * hAeff;

  TFile * fflute=new TFile(inputpath);
  fflute->ls();
  MStatusArray* msaflute = (MStatusArray*)fflute->Get("MStatusDisplay");
  TCanvas *tcAeff = (TCanvas*)msaflute->FindCanvas("Coll. Area Eest");
  for(int i=0;i<(tcAeff->GetListOfPrimitives()->GetSize());i++)
  {
    cout<<"object"<< i<<" : " << tcAeff->GetListOfPrimitives()->At(i)->GetName()<<endl;
    if(tcAeff->GetListOfPrimitives()->At(i)->IsA()->InheritsFrom(TH1D::Class()))
    {
      hAeff = (TH1D*)tcAeff->GetListOfPrimitives()->At(i);

      continue;
    }
    //   if(!strcmp(tcAeff->GetListOfPrimitives()->At(i)->GetName(),"fAbsorbedSED"))
    //     fsed =(TF1*)tcAeff->GetListOfPrimitives()->At(i);
  }

  // TCanvas * can = new TCanvas("can","can",800,600);
  // hAeff->Draw();
  tcAeff->Draw();
  hAeff->Print("all");
  cout<<hAeff->GetBinContent(38)<<endl;
  cout<<"          lowedge, center, content, error"<<endl;
  cout<<"Bin0  :"<<0 <<": "<<hAeff->GetBinLowEdge(0)   <<" " <<hAeff->GetBinCenter(0) <<" " <<hAeff->GetBinLowEdge(0) +hAeff->GetBinWidth(0) <<" " <<hAeff->GetBinContent(0) <<" " <<hAeff->GetBinError(0) <<endl;
  cout<<"Bin1  :"<<1 <<": "<<hAeff->GetBinLowEdge(1)   <<" " <<hAeff->GetBinCenter(1) <<" " <<hAeff->GetBinLowEdge(1) +hAeff->GetBinWidth(1) <<" " <<hAeff->GetBinContent(1) <<" " <<hAeff->GetBinError(1) <<endl;
  cout<<"Bin2  :"<<2 <<": "<<hAeff->GetBinLowEdge(2)   <<" " <<hAeff->GetBinCenter(2) <<" " <<hAeff->GetBinLowEdge(2) +hAeff->GetBinWidth(2) <<" " <<hAeff->GetBinContent(2) <<" " <<hAeff->GetBinError(2) <<endl;
  cout<<"Bin3  :"<<3 <<": "<<hAeff->GetBinLowEdge(3)   <<" " <<hAeff->GetBinCenter(3) <<" " <<hAeff->GetBinLowEdge(3) +hAeff->GetBinWidth(3) <<" " <<hAeff->GetBinContent(3) <<" " <<hAeff->GetBinError(3) <<endl;
  cout<<"Bin11 :"<<11<<": "<<hAeff->GetBinLowEdge(11)<<" " <<hAeff->GetBinCenter(11)<<" " <<hAeff->GetBinLowEdge(11)+hAeff->GetBinWidth(11)<<" " <<hAeff->GetBinContent(11)<<" " <<hAeff->GetBinError(11)<<endl;
  cout<<"Bin12 :"<<12<<": "<<hAeff->GetBinLowEdge(12)<<" " <<hAeff->GetBinCenter(12)<<" " <<hAeff->GetBinLowEdge(12)+hAeff->GetBinWidth(12)<<" " <<hAeff->GetBinContent(12)<<" " <<hAeff->GetBinError(12)<<endl;
  cout<<"Bin13 :"<<13<<": "<<hAeff->GetBinLowEdge(13)<<" " <<hAeff->GetBinCenter(13)<<" " <<hAeff->GetBinLowEdge(13)+hAeff->GetBinWidth(13)<<" " <<hAeff->GetBinContent(13)<<" " <<hAeff->GetBinError(13)<<endl;


  ofstream fout;
  // fout.open("collarea_magic.txt");
  // fout<<" lowedge, center, content, error"<<endl;
  // for (int i = 11; i<39;i++)
  // {
  //   fout<<hAeff->GetBinLowEdge(i)<<", " <<hAeff->GetBinCenter(i)<<", " <<hAeff->GetBinContent(i)<<", " <<hAeff->GetBinError(i)<<endl;
  // }
  // fout<<hAeff->GetBinLowEdge(39)+hAeff->GetBinWidth(39)<<",,,"<<endl;
  

  fout.open("collarea_magic.ecsv");
  fout<<"# %ECSV 0.9"<<endl;
  fout<<"# ---"<<endl;
  fout<<"# datatype:"<<endl;
  fout<<"# - {name: e_min, unit: GeV, datatype: float64}"<<endl;
  fout<<"# - {name: e_ref, unit: GeV, datatype: float64}"<<endl;
  fout<<"# - {name: e_max, unit: GeV, datatype: float64}"<<endl;
  fout<<"# - {name: Aeff, unit: cm2, datatype: float64}"<<endl;
  fout<<"# - {name: Aeff_err, unit: cm2, datatype: float64}  "<<endl;
  fout<<"e_min e_ref e_max Aeff Aeff_err"<<endl;
  for (int i = 11; i<39;i++)
  {
    fout
      <<hAeff->GetBinLowEdge(i)<<" " 
      <<hAeff->GetBinCenter(i)<<" " 
      <<hAeff->GetBinLowEdge(i)+hAeff->GetBinWidth(i)<<" "
      <<hAeff->GetBinContent(i)<<" " 
      <<hAeff->GetBinError(i)<<endl;
  }

  fout.close();
  return 0;
  

}



