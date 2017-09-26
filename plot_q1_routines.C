#include <iostream>
#include <fstream>
#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <iterator>
#include <vector>
#include <algorithm>
using namespace std;

TF1 *neg_field_corr;
TF1 *pos_field_corr;
vector<double> angle;
Double_t fit_function(Double_t *x,Double_t *par);
void Get_two_col_data_from_file(TString filename,vector<double> &,vector<double> & );
void Get_three_col_data_from_file(TString filename,vector<double> &,vector<double> & ,vector<double> & );
void Plot_q1_zscan(TString fname);
void Set_probe();
void Read_zpos_file();
Double_t Get_zpos(Double_t num);
void Read_angle_file();
Double_t field_corr(Double_t field );
void Plot_q1_bi(vector<double> & ,vector<double> &,vector<double> &);
void run_plot_q1_bi();
vector<double> Q1_central_current_data;
vector<double> Q1_central_field_data;
vector<double> zpos_num;
vector<double> zpos_val;
void Plot_z0_data(TString fname);
void run_plot_q1_zscan(TString fname);
//
void run_plot_q1_zscan(TString fname) {
  Set_probe();
  Get_two_col_data_from_file("table_zpos_num_versus_dis.dat",zpos_num,zpos_val);
  Plot_q1_zscan(fname);
}
//
Double_t Get_zpos(Double_t num) {
  Double_t ret_val=0.;
 for (UInt_t i=0; i<zpos_num.size(); i++) {
   if (num==zpos_num[i]) ret_val=zpos_val[i];
 }
 return ret_val;
} 
//
void run_plot_q1_bi() {
  Set_probe();
  Read_angle_file();
  Double_t measure_rad=.1876 ; // Radius of probe for measurement
  vector<double> Q1_curratio_data;
  vector<double> Q1_radratio_data;
  const int nftot=4;
  TString cname[nftot]={"1228","1672","2169","2454"};
  Double_t Q1_tosca_curratio[nftot]={0.00353,0.00352,0.00342,0.00330};
  vector<double> Q1_tosca;
  for (int nf=0;nf<nftot;nf++) {
    Q1_central_current_data.push_back(cname[nf].Atof());
    Plot_z0_data(cname[nf]);
    Q1_central_field_data[nf]= abs(Q1_central_field_data[nf])/10.;
    Q1_radratio_data.push_back(Q1_central_field_data[nf]/measure_rad);
    Q1_curratio_data.push_back(Q1_radratio_data[nf]/Q1_central_current_data[nf]);
    Q1_tosca.push_back(Q1_tosca_curratio[nf]);
    cout << Q1_central_current_data[nf] << " " << Q1_central_field_data[nf] <<  " " << Q1_radratio_data[nf] <<  " " << Q1_curratio_data[nf] << " " << Q1_tosca[nf] << " " << Q1_curratio_data[nf]/Q1_tosca[nf] << endl;
    Q1_tosca[nf]=Q1_tosca[nf]*1.0;
  }
  Plot_q1_bi(Q1_central_current_data,Q1_curratio_data,Q1_tosca);
 //
}
//
//
 void Plot_q1_bi(vector<double> &x,vector<double> &y,vector<double> &ty) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(0.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.2);
 //
 TF1 *cur_field;
 //
vector<double> bvec;
 Int_t vecsize=x.size();
 for (UInt_t i=0; i<x.size(); i++) {
   bvec.push_back(y[i]*x[i]);
   cout << bvec[i] << endl;
 }
TCanvas *can4 = new TCanvas("can4","Q1 I versus B",800,800);
 can4->Divide(1,1);
 can4->cd(1);
 TGraph *grcur = new TGraph(vecsize,&(bvec[0]),&(x[0]));
 grcur->SetMarkerStyle(21);
 grcur->SetTitle("Q1 I versus B");
 grcur->GetYaxis()->SetTitle("Current [A]");
 grcur->GetXaxis()->SetTitle("B [T]");
 grcur->SetMarkerSize(1.0);
 grcur->SetMaximum(4000);
 grcur->SetMinimum(0);
 grcur->Draw("AP");
   grcur->Fit("pol3","QO","",0.,10.);
  cur_field = grcur->GetFunction("pol3");
  Double_t ftest=19.9049*6.4/11./1.844*18./20.;
  cout << " current for b field = " << ftest << " " << cur_field->Eval(ftest,0.,0.) << endl;
   //
TCanvas *can = new TCanvas("can","Q1 (T/m)/Current versus Current",800,800);
 TPad *pad1=new TPad("pad1","",0.05,0.05,0.98,0.98);
 can->Divide(1,1);
 can->cd(1);
 pad1->Draw();
 pad1->cd();
 pad1->SetLeftMargin(0.2);
 TGraph *grposup = new TGraph(vecsize,&(x[0]),&(y[0]));
 TGraph *grtosca = new TGraph(vecsize,&(x[0]),&(ty[0]));
 grposup->SetName("grposup");
 grtosca->SetName("grtosca");
 grposup->Draw("AP");
 // grtosca->Draw("L same");
 grposup->SetMarkerStyle(21);
 grposup->SetTitle("Q1 (B/r)/I versus I");
 grposup->GetXaxis()->SetTitle("Current (A)");
 grposup->GetYaxis()->SetTitle("(B/r)/I  [(T/m)/A]");
 grposup->SetMarkerSize(1.5);
 grposup->SetMaximum(.0037);
 grposup->SetMinimum(.0027);
 // grposup->Fit("pol2","O","",x[0],4000.);
 
 // can->Close();
 TLegend *leg = new TLegend(.5,.5,.8,.6);
TCanvas *can2 = new TCanvas("can2","Q1 (T/m)/Current versus Current",800,800);
 can2->Divide(1,1);
 can2->cd(1);
 grposup->Draw("AP");
 grtosca->Draw("P same");
  grtosca->SetMarkerStyle(22);
 grtosca->SetMarkerSize(1.5);
 leg->AddEntry("grposup","Data","p");
 leg->AddEntry("grtosca","Tosca","p");
 leg->Draw();

}
//
void Plot_z0_data(TString fname) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
//
  //
 //
vector<double> q1_ang_index;
vector<double> q1_current;
vector<double> q1_field;
vector<double> q1_angle;
 TString dname;
 dname="q1-"+fname+"-z-24-ang_scan.dat";
  Get_three_col_data_from_file(dname,q1_ang_index,q1_field,q1_current);
  //
 for (UInt_t i=0; i<q1_ang_index.size(); i++) {
   q1_angle.push_back(angle[int(q1_ang_index[i])-1]);
   q1_field[i]=q1_field[i]+field_corr(q1_field[i]);
 }
 //
 TString title;
 title="Q1 I="+fname+"A Z = 0 ";
  TCanvas *can1 = new TCanvas("can1",title,800,800);
 can1->Divide(1,1);
 can1->cd(1);
 TGraph *grcir1 = new TGraph(q1_angle.size(),&(q1_angle[0]),&(q1_field[0]));
 grcir1->Draw("AP");
 grcir1->SetMarkerStyle(21);
 grcir1->SetTitle(title);
 grcir1->GetXaxis()->SetTitle("Angle (radians)");
 grcir1->GetYaxis()->SetTitle("Field (kG)");
 grcir1->SetMarkerSize(.5);
 grcir1->SetMaximum(12.);
 grcir1->SetMinimum(-12.);
 TF1 *f1 = new TF1("fitfunc",fit_function,-3.,3.,10);
 grcir1->Fit("fitfunc","Q");
 f1->Draw("same");
 Q1_central_field_data.push_back(f1->GetParameter(0));
 can1->WaitPrimitive();
 can1->Close();
  //
}
//
Double_t field_corr(Double_t field )
{
  Double_t fieldcorr;
  fieldcorr=0.;
   if (field<0) fieldcorr=neg_field_corr->Eval(field,0.,0.);
   if (field >0) fieldcorr=-pos_field_corr->Eval(field,0.,0.);
  return fieldcorr;
}
//
void Get_two_col_data_from_file(TString filename,vector<double> &x,vector<double> &y )
{
 ifstream corfile;
 corfile.open(filename);
 TString curline;
 Int_t nline=0;
 while (corfile.good())
   {
     curline.ReadLine(corfile,kTRUE);
    TString sc=curline.Data();
    Int_t chi,ncomma=sc.CountChar(',');
    if (ncomma ==1) {
        chi=sc.Index(',');
        TString temp(sc(0,chi));
	x.push_back(temp.Atof());
        temp=sc(chi+1,sc.Sizeof()+1);
	y.push_back(temp.Atof());
	//  cout << " data " << x[nline] << " "<< y[nline] <<endl;
         nline++;
    }
   }
}
//
void Get_three_col_data_from_file(TString filename,vector<double> &x,vector<double> &y,vector<double> &dy )
{
 ifstream corfile;
 corfile.open(filename);
 TString curline;
 Int_t nline=0;
 while (corfile.good())
   {
     curline.ReadLine(corfile,kTRUE);
    TString sc=curline.Data();
    Int_t chi,ncomma=sc.CountChar(',');
    if (ncomma ==2) {
        chi=sc.Index(',');
        TString temp(sc(0,chi));
	x.push_back(temp.Atof());
        temp=sc(chi+1,sc.Sizeof()-chi);
	//cout << temp  << endl;
        chi=temp.Index(',');
        TString temp1(temp(0,chi));
        TString temp2(temp(chi+1,temp.Sizeof()-chi));
	//	cout << temp << " " << chi << endl;
	y.push_back(temp1.Atof());
        temp=sc(chi+1,temp.Sizeof()-chi);
	dy.push_back(temp2.Atof());
	//cout << " data " << x[nline] << " "<< y[nline]<< " " <<  dy[nline]<<endl;
         nline++;
    }
   }
 corfile.close();
}
//
Double_t fit_function(Double_t *x,Double_t *par) {
  Double_t xx=x[0];
  Double_t f=0;
  UInt_t norder=10;
   for (UInt_t i=0; i<5; i++) {
     f= f+ par[i]*TMath::Sin(2*(i+1)*xx) ;
    }  
   for (UInt_t i=5; i<norder; i++) {
     f= f+ par[i]*TMath::Cos(2*(i-5+1)*xx) ;
    }  
  return f;
}
//
void Read_angle_file() {
vector<double> ang_index;
vector<double> xpos;
vector<double> ypos;
  Get_three_col_data_from_file("ang_index_versus_pos.dat",ang_index,xpos,ypos);
 for (UInt_t i=0; i<ang_index.size(); i++) {
   angle.push_back(TMath::ATan2(xpos[i],ypos[i]));
 }
}
//
//
void Read_zpos_file() {
}
//
//
void Plot_q1_zscan(TString fname)  {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //
vector<double> q1_zpos;
vector<double> q1_current;
vector<double> q1_field;
//
  // Assume that files has zpos as number,
  Get_three_col_data_from_file(fname,q1_zpos,q1_field,q1_current);
 //
  Double_t central_field=0.;
  Double_t zpos_min=100.;
  Double_t IntBdl=0.;
 for (UInt_t i=0; i<q1_zpos.size(); i++) {
   //cout << q1_zpos[i] << " " << Get_zpos(q1_zpos[i]) << endl;
   q1_zpos[i]=(Get_zpos(q1_zpos[i]))/10.;
   q1_field[i]=q1_field[i]+field_corr(q1_field[i]);
 }
 for (UInt_t i=1; i<q1_zpos.size()-1; i++) {
   IntBdl = IntBdl + TMath::Abs(q1_field[i])*(q1_zpos[i+1]-q1_zpos[i]);
   if ( TMath::Abs(q1_zpos[i]) < zpos_min) {
      central_field=TMath::Abs(q1_field[i]);
      zpos_min = TMath::Abs(q1_zpos[i]);
   }
 }
 Double_t leff=IntBdl/central_field;
 cout << " Leff = " << leff << " Central field = " << central_field << " Bdl (kG-cm) = " << IntBdl << endl;
 //
TCanvas *can = new TCanvas("can"," ",800,800);
 can->Divide(1,1);
 can->cd(1);
 TGraph *grq1 = new TGraph(q1_zpos.size(),&(q1_zpos[0]),&(q1_field[0]));
 grq1->Draw("AP");
 grq1->SetMarkerStyle(21);
 grq1->SetTitle("");
 grq1->GetXaxis()->SetTitle("Z pos (cm)");
 grq1->GetYaxis()->SetTitle("Field (kG)");
 grq1->SetMarkerSize(.5);
}
//
void Plot_q1_200A_zscan() {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //
vector<double> q1_zpos;
vector<double> q1_field;
//
  Get_two_col_data_from_file("q1-200A-z-scan-ang-9.dat",q1_zpos,q1_field);
 //
  Double_t central_field=0.;
  Double_t zpos_min=100.;
  Double_t IntBdl=0.;
 for (UInt_t i=0; i<q1_zpos.size(); i++) {
   q1_zpos[i]=(q1_zpos[i])/10.;
   q1_field[i]=(q1_field[i])/1000.+field_corr(q1_field[i]/1000.);
 }
 for (UInt_t i=1; i<q1_zpos.size()-1; i++) {
   IntBdl = IntBdl + TMath::Abs(q1_field[i])*(q1_zpos[i+1]-q1_zpos[i]);
   if ( TMath::Abs(q1_zpos[i]) < zpos_min) {
      central_field=TMath::Abs(q1_field[i]);
      zpos_min = TMath::Abs(q1_zpos[i]);
   }
 }
 Double_t leff=IntBdl/central_field;
 cout << " Leff = " << leff << " Central field = " << central_field << " Bdl (kG-cm) = " << IntBdl << endl;
 //
TCanvas *can = new TCanvas("can","Q1 I=200A ang pos number 9 ",800,800);
 can->Divide(1,1);
 can->cd(1);
 TGraph *grq1 = new TGraph(q1_zpos.size(),&(q1_zpos[0]),&(q1_field[0]));
 grq1->Draw("AP");
 grq1->SetMarkerStyle(21);
 grq1->SetTitle("Q1 I=200A ang pos number 9");
 grq1->GetXaxis()->SetTitle("Z pos (cm)");
 grq1->GetYaxis()->SetTitle("Field (kG)");
 grq1->SetMarkerSize(.5);
 grq1->SetMaximum(0.);
 grq1->SetMinimum(-2.);

}
//
void Set_probe() {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //
vector<double> nominal_field_corfile_data;
vector<double> corr_field_corfile_data;
vector<double> resid;
 Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
 for (UInt_t i=0; i<corr_field_corfile_data.size(); i++) {
   corr_field_corfile_data[i]=corr_field_corfile_data[i]/1000.;
 }
 //
 Double_t fit_val;
TCanvas *cprobe = new TCanvas("cprobe","Lakeshore probe B versus B_corr",800,800);
 cprobe->Divide(1,2);
 cprobe->cd(1);
 TGraph *grneg = new TGraph(corr_field_corfile_data.size(),&(nominal_field_corfile_data[0]),&(corr_field_corfile_data[0]));
 TGraph *grpos = new TGraph(corr_field_corfile_data.size(),&(nominal_field_corfile_data[0]),&(corr_field_corfile_data[0]));
 grneg->Draw("AP");
 grneg->SetMarkerStyle(21);
 grneg->SetTitle("Lakeshore probe B correction versus B");
 grneg->GetXaxis()->SetTitle("Field (kG)");
 grneg->GetYaxis()->SetTitle("Field correction (kG)");
 grneg->SetMarkerSize(.5);
 grneg->SetMaximum(.3);
 grneg->SetMinimum(-.3);
 grneg->Fit("pol8","Q","",-30.0,0.0);
 neg_field_corr = grneg->GetFunction("pol8");
 grpos->Fit("pol6","Q","",0.0,30.0);
 pos_field_corr = grpos->GetFunction("pol6");
 for (UInt_t i=0;i<corr_field_corfile_data.size();i++) {
   if (nominal_field_corfile_data[i] <0) fit_val=neg_field_corr->Eval(nominal_field_corfile_data[i],0.,0.);
   if (nominal_field_corfile_data[i] >0) fit_val=pos_field_corr->Eval(nominal_field_corfile_data[i],0.,0.);
   //   cout << fit_val << " " << yvec[i] << endl;
   resid.push_back((fit_val-corr_field_corfile_data[i])/corr_field_corfile_data[i]);
 }
 cprobe->cd(2);
 TGraph *gr_resid = new TGraph(corr_field_corfile_data.size(),&(nominal_field_corfile_data[0]),&(resid[0]));
 gr_resid->Draw("AP");
 gr_resid->SetTitle("Field Corr FIT Residual versus Current");
 gr_resid->GetXaxis()->SetTitle("Field (kG)");
 gr_resid->GetYaxis()->SetTitle("(FIT-Meas)/Meas");
 gr_resid->SetMarkerStyle(21);
 gr_resid->SetMarkerSize(.5);
 gr_resid->SetMaximum(.1);
 gr_resid->SetMinimum(-.1);
 cprobe->Close();
 //
}

