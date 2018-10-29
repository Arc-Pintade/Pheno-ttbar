#define analyze_cxx
#include "../include/analyze.hpp"
#include "../include/Const.hpp"

#include <iostream>
#include <fstream>
#include <vector>

#include <TString.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;

// cout for TMatrixD

void showMatrix(TMatrixD m){
    cout<<"matrix is : "<<endl;
    for(int i=0; i<4; i++)
        cout<<m(i,0)<<" "<<m(i,1)<<" "<<m(i,2)<<" "<<m(i,3)<<endl;
    cout<<endl;
}

//**********************************************************************************//
//************************* class EventsAnalyze definition *************************//
//**********************************************************************************//

//__________________ Constructor __________________//

EventsAnalyze::EventsAnalyze(int nEvents_user, string VMAPqqbar_user, string VMAP2g_user, string VMAF_user){

// Without Resize matrices are uncompatible
    averagePqqbar.ResizeTo(4,4,-1);
    averageP2g.ResizeTo(4,4,-1);
    averageF.ResizeTo(4,4,-1);

    nEvents         = nEvents_user;

    nPqqbar         = calculateEventNumber(VMAPqqbar_user);
    nP2g            = calculateEventNumber(VMAP2g_user);
    nF              = calculateEventNumber(VMAF_user);

    averagePqqbar   = calculateAverageMatrix(VMAPqqbar_user);
    averageP2g      = calculateAverageMatrix(VMAP2g_user);
    averageF        = calculateAverageMatrix(VMAF_user);
}

//__________________ Read Values __________________//

int EventsAnalyze::getnF(){return nF;}
int EventsAnalyze::getnPqqbar(){return nPqqbar;}
int EventsAnalyze::getnP2g(){return nP2g;}
TMatrixD EventsAnalyze::getAverageF(){return averageF;}
TMatrixD EventsAnalyze::getAveragePqqbar(){return averagePqqbar;}
TMatrixD EventsAnalyze::getAverageP2g(){return averageP2g;}


//________________ Calculated Values ______________//

// Number of events for a given cut 

int EventsAnalyze::calculateEventNumber(string s){
    fstream f(s, ios::in);

    int foo = nEvents;
    TMatrixD null(4,4);
    vector <TMatrixD> VM(nEvents);
    for(int i=0; i<nEvents; i++){
        VM[i].ResizeTo(4,4,-1);
    }

    for(int i=0; i<nEvents; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++)
                f>>VM[i].operator[](j).operator[](k);
        }
        if(VM[i]==null)
            foo--;
    }

    f.close();
    return foo;
}

// Average Amunu Matrices

TMatrixD EventsAnalyze::calculateAverageMatrix(string s){
    fstream f(s, ios::in);
    double n = (double)calculateEventNumber(s);

    TMatrixD foo(4,4);

    TMatrixD null(4,4);
    vector <TMatrixD> VM(nEvents);
    for(int i=0; i<nEvents; i++){
        VM[i].ResizeTo(4,4,-1);
    }

    for(int i=0; i<nEvents; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++)
                f>>VM[i].operator[](j).operator[](k);
        }
        foo += VM[i];
    }
    foo = (1/n) * foo;

    f.close();
    return foo;
}

// text file
void EventsAnalyze::writeNumMatrix(string s, int n, TMatrixD m){
    std::ofstream f("data/matrix/"+s+".txt", std::ios::out);
    f<<s<<" "<<endl<<endl;
    f<<"Events n= "<<n<<endl<<endl;
    for(int j =0; j<4; j++)
        f<<m(j,0)<<" "<<m(j,1)<<" "<<m(j,2)<<" "<<m(j,3)<<std::endl;
    f.close();
}

//**********************************************************************************//
//**********************************************************************************//


//**********************************************************************************//
//************************** class FunctionAnalyze definition **********************//
//**********************************************************************************//

//__________________ Constructor __________________//

FunctionAnalyze::FunctionAnalyze(int nPqqbar_user, int nP2g_user, int nF_user, TMatrixD averagePqqbar_user, TMatrixD averageP2g_user, TMatrixD averageF_user){

// Without Resize matrices are uncompatible
    averagePqqbar.ResizeTo(4,4,-1);
    averageP2g.ResizeTo(4,4,-1);
    averageF.ResizeTo(4,4,-1);
    averageA.ResizeTo(4,4,-1);
    averageAD0.ResizeTo(4,4,-1);

    nPqqbar         = nPqqbar_user;
    nP2g            = nP2g_user;
    nF              = nF_user;

    averagePqqbar   = averagePqqbar_user;
    averageP2g      = averageP2g_user;
    averageF        = averageF_user;

    averageA        = calculateAverageMatrixA(averagePqqbar, averageP2g, averageF, nPqqbar, nP2g, nF);
    averageAD0      = calculateAverageMatrixAD0(averagePqqbar, averageF);

    a1              = calculateCoefficent(1, averageA);
    a2              = calculateCoefficent(2, averageA);
    a3              = calculateCoefficent(3, averageA);
    a4              = calculateCoefficent(4, averageA);
    a5              = calculateCoefficent(5, averageA);
    a1D0            = calculateCoefficentD0(1, averageAD0);
    a2D0            = calculateCoefficentD0(2, averageAD0);
    a3D0            = calculateCoefficentD0(3, averageAD0);
    a4D0            = calculateCoefficentD0(4, averageAD0);
    a5D0            = calculateCoefficentD0(5, averageAD0);

}

//__________________ Read Values __________________//

int FunctionAnalyze::getnF(){return nF;}
int FunctionAnalyze::getnPqqbar(){return nPqqbar;}
int FunctionAnalyze::getnP2g(){return nP2g;}
TMatrixD FunctionAnalyze::getAverageF(){return averageF;}
TMatrixD FunctionAnalyze::getAveragePqqbar(){return averagePqqbar;}
TMatrixD FunctionAnalyze::getAverageP2g(){return averageP2g;}
TMatrixD FunctionAnalyze::getAverageA(){return averageA;}
TMatrixD FunctionAnalyze::getAverageAD0(){return averageAD0;}

double FunctionAnalyze::geta1(){return a1;}
double FunctionAnalyze::geta2(){return a2;}
double FunctionAnalyze::geta3(){return a3;}
double FunctionAnalyze::geta4(){return a4;}
double FunctionAnalyze::geta5(){return a5;}
double FunctionAnalyze::geta1D0(){return a1D0;}
double FunctionAnalyze::geta2D0(){return a2D0;}
double FunctionAnalyze::geta3D0(){return a3D0;}
double FunctionAnalyze::geta4D0(){return a4D0;}
double FunctionAnalyze::geta5D0(){return a5D0;}

// Text input
TMatrixD FunctionAnalyze::readMatrix(string s){
    TMatrixD foo(4,4);
    string tmp;
    int tmp2;
    fstream f("data/matrix/"+s+".txt", ios::in);
    f>>tmp>>tmp>>tmp;
    f>>tmp2;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            f>>foo(i,j);
    return foo;
}

int FunctionAnalyze::readNumber(string s){
    int foo;
    string tmp;
    fstream f("data/matrix/"+s+".txt", ios::in);
    f>>tmp>>tmp>>tmp;
    f>>foo;
    return foo;
}

//________________ Calculated Values ______________//

// Matrix A
TMatrixD FunctionAnalyze::calculateAverageMatrixA(TMatrixD Pqqbar, TMatrixD P2g, TMatrixD F, int nP1, int nP2, int n){
    TMatrixD foo(4,4);
    foo = ((double)nP1 / (double)n) * 0.5*Pqqbar + ((double)nP2 / (double)n) * 0.5*P2g + F;
    return foo;
}

TMatrixD FunctionAnalyze::calculateAverageMatrixAD0(TMatrixD Pqqbar, TMatrixD F){
    TMatrixD foo(4,4);
    foo = 0.5*Pqqbar + F;
    return foo;
}


// a and fSME fonctions values
double FunctionAnalyze::calculateCoefficent(int i, TMatrixD m){
    double foo = 0;

    if(i==1)
        foo = (s1*s1*s2*s2 + c1*c1) * m(0,0) + (s1*s1*c2*c2) * m(2,2);
    else if(i==2)
        foo = (c2*c2) * m(0,0) + (s2*s2) * m(2,2);
    else if(i==3)
        foo = (s1*c2*s2) * (m(2,2) - m(0,0));
    else if(i==4)
        foo = (c1*s1*c2*c2) * (m(2,2) - m(0,0));
    else if(i==5)
        foo = (c2*c1*s2) * (m(2,2) - m(0,0));
    else{
        cout<<"error with coefficient"<<endl;
    }
    return foo;
}

double FunctionAnalyze::calculateCoefficentD0(int i, TMatrixD m){
    double foo = 0;
    if(i==1)
        foo = (c1D0*c1D0*c2D0*c2D0 + s1D0*s1D0) * m(0,0) + (c1D0*c1D0*s2D0*s2D0) * m(2,2);
    else if(i==2)
        foo = (s2D0*s2D0) * m(0,0) + (c2D0*c2D0) * m(2,2);
    else if(i==3)
        foo = (c1D0*c2D0*s2D0) * (m(2,2) - m(0,0));
    else if(i==4)
        foo = (c1D0*s1D0*s2D0*s2D0) * (m(0,0) + m(2,2));
    else if(i==5)
        foo = (c2D0*s1D0*s2D0) * (m(0,0) - m(2,2));
    else{
        cout<<"error with coefficient"<<endl;
    }
    return foo;
}

double FunctionAnalyze::fSME(string s, double cmunu, double t){
    double foo = 0;
    if(s=="XX" || s=="YY")
        foo = 2 * cmunu * (((a1-a2)/2) * cos(2*omega*t) + a3 * sin(2*omega*t));
    else if(s=="XY" || s=="YX")
        foo = 2 * cmunu * (((a1-a2)/2) * sin(2*omega*t) - a3 * cos(2*omega*t));
    else if(s=="XZ" || s=="ZX")
        foo = 2 * cmunu * (a4 * cos(omega*t) + a5 * sin(omega*t));
    else if(s=="YZ" || s=="ZY")
        foo = 2 * cmunu * (a4 * sin(omega*t) - a5 * cos(omega*t));
    else
        cout<<"error f fonction"<<endl;
    return foo;
}

double FunctionAnalyze::fSMED0(string s, double cmunu, double t){
    double foo = 0;
    if(s=="XX" || s=="YY")
        foo = 2 * cmunu * (((a1D0-a2D0)/2) * cos(2*omega*t + deph) + a3D0 * sin(2*omega*t + deph));
    else if(s=="XY" || s=="YX")
        foo = 2 * cmunu * (((a1D0-a2D0)/2) * sin(2*omega*t + deph) - a3D0 * cos(2*omega*t + deph));
    else if(s=="XZ" || s=="ZX")
        foo = 2 * cmunu * (a4D0 * cos(omega*t + deph) + a5D0 * sin(omega*t + deph));
    else if(s=="YZ" || s=="ZY")
        foo = 2 * cmunu * (a4D0 * sin(omega*t + deph) - a5D0 * cos(omega*t + deph));
    else
        cout<<"error f fonction"<<endl;
    return foo;
}

//________________ Graphicals analysis __________________//

// fSME fuction of time
void FunctionAnalyze::quadriHisto(TString s, double cmunu, int maxTime, int timeStep){
    TCanvas* w = new TCanvas("f functions histogram "+s,"",200,10,800,600);
    TH1F* h1 = new TH1F("c_{XX} = - c_{YY}"+s,"", timeStep, 0, maxTime);
    TH1F* h2 = new TH1F("c_{XY} = c_{YX}"+s,"", timeStep, 0, maxTime);
    TH1F* h3 = new TH1F("c_{XZ} = c_{ZX}"+s,"", timeStep, 0, maxTime);
    TH1F* h4 = new TH1F("c_{YZ} = c_{ZY}"+s,"", timeStep, 0, maxTime);
    for (int i=0; i<timeStep; i++){
        h1->SetBinContent(i+1, fSME("XX", cmunu, i));
        h2->SetBinContent(i+1, fSME("XY", cmunu, i));
        h3->SetBinContent(i+1, fSME("XZ", cmunu, i));
        h4->SetBinContent(i+1, fSME("YZ", cmunu, i));
    }

    h1->Draw();
    h2->Draw("SAME");
    h3->Draw("SAME");
    h4->Draw("SAME");

    h1->GetYaxis()->SetTitle("f_{SME}(t)");
    h1->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");

    h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);
    h2->Write();   h2->SetLineWidth(2);   h2->SetLineColor(kMagenta);
    h3->Write();   h3->SetLineWidth(2);   h3->SetLineColor(kBlue);
    h4->Write();   h4->SetLineWidth(2);   h4->SetLineColor(kGreen);

    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram "+s,"C"); // option "C" allows to center the header
    legend->AddEntry(h1,"c_{XX} = - c_{YY}","l"); // option "l" is for line (form of legend)
    legend->AddEntry(h2,"c_{XY} = c_{YX}","l");
    legend->AddEntry(h3,"c_{XZ} = c_{ZX}","l");
    legend->AddEntry(h4,"c_{YZ} = c_{ZY}","l");
    legend->Draw();

    w->Update();
    w->SaveAs("results/"+s+"fComparaison.png");
}


void FunctionAnalyze::signalView(TString s, TString XX, double cmunu, int nbin, double nEventBkgd, double nEventTTbarSM){
    TCanvas* c = new TCanvas("Input"+s+" f"+XX,"Input"+s, 10,10,1000,800);
    double norm = (nEventBkgd+nEventTTbarSM)/((double)nbin);

    TH1F* hTTbarSM = new TH1F("t#bar{t} events"+s,"",nbin,0,nbin);
    for (int it=0; it<nbin; it++) 
        hTTbarSM->SetBinContent(it+1,(nEventBkgd+nEventTTbarSM)/((double)nbin));
    hTTbarSM->Write();
    hTTbarSM->Draw("SAME");
    hTTbarSM->SetFillColor(kRed+1);
    hTTbarSM->GetYaxis()->SetTitle("Events");
    hTTbarSM->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");

    TH1F* hBkgd = new TH1F("Background events","",nbin,0,nbin);
    for (int it=0; it<nbin; it++) 
        hBkgd->SetBinContent(it+1,nEventBkgd/((double)nbin));
    hBkgd->Write();
    hBkgd->Draw("SAME");
    hBkgd->SetFillColor(kBlue+1);

    TH1F* hSME;
    if(XX=="XX" || XX=="YY"){
        hSME = new TH1F("c_{xx} = - c_{yy}"+s,"",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XX", cmunu, it*3600));
        }
    }
    else if(XX=="XY" || XX=="YX"){
        hSME = new TH1F("c_{xy} = c_{yx}"+s,"",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XY", cmunu, it*3600));
        }
    }
    else if(XX=="XZ" || XX=="ZX"){
        hSME = new TH1F("c_{xz} = c_{zx}"+s,"",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XZ", cmunu, it*3600));
        }
    }
    else if(XX=="YZ" || XX=="ZY"){
        hSME = new TH1F("c_{yz} = c_{zy}"+s,"",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("YZ", cmunu, it*3600));
        }
    }
    else
        cout<<"error with f fuction"<<endl;
    hSME->Write();
    hSME->Draw("SAME");
    hSME->SetLineColor(kGreen+1);
    hSME->SetLineWidth(2);

    hTTbarSM->SetStats(0);

    TH1F* hAzimov = new TH1F("Azimov"+s, "", nbin, 0, nbin);
    for (int it=0; it<nbin; it++){
        hAzimov->SetBinContent(it+1,(nEventBkgd+nEventTTbarSM)/((double)nbin));
        hAzimov->SetBinError(it+1, sqrt((nEventBkgd+nEventTTbarSM)/((double)nbin)));
    }
    hAzimov->Write();
    hAzimov->Draw("E SAME");

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) signal histogram "+s,"C"); // option "C" allows to center the header
    legend->AddEntry(hTTbarSM,"t#bar{t} signal","f"); // option "f" is for fill (form of legend)
    legend->AddEntry(hBkgd,"Background ignal","f");
    legend->AddEntry(hSME,"Theoritical SME signal","l");
    legend->AddEntry(hAzimov,"Statistical errors","lep");
    legend->Draw();

    c->SaveAs("results/"+s+XX+"signal.png");
}

void FunctionAnalyze::signalHisto(TString s, TString XX, double cmunu, int nbin, double nEventTTbarSM){
    TCanvas* c = new TCanvas("Input"+s+" f"+XX,"Input"+s, 10,10,800,600);
    double norm = (nEventTTbarSM)/((double)nbin);

    TH1F* hSME;
    if(XX=="XX" || XX=="YY"){
        hSME = new TH1F("c_{xx} = - c_{yy}","",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XX", cmunu, it*3600));
        }
    }
    else if(XX=="XY" || XX=="YX"){
        hSME = new TH1F("c_{xy} = c_{yx}","",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XY", cmunu, it*3600));
        }
    }
    else if(XX=="XZ" || XX=="ZX"){
        hSME = new TH1F("c_{xz} = c_{zx}","",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("XZ", cmunu, it*3600));
        }
    }
    else if(XX=="YZ" || XX=="ZY"){
        hSME = new TH1F("c_{yz} = c_{zy}","",nbin,0,nbin);
        for (int it=0; it<nbin; it++) {
            hSME->SetBinContent(it+1, norm + norm*fSME("YZ", cmunu, it*3600));
        }
    }
    else
        cout<<"error with f fuction"<<endl;
    hSME->Write();
    hSME->Draw();
    hSME->SetLineColor(kBlue);
    hSME->SetLineWidth(2);

    hSME->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) binned signal "+s,"C"); // option "C" allows to center the header
    legend->AddEntry(hSME,"SME signal","f"); // option "f" is for fill (form of legend)
    legend->Draw();
}

void FunctionAnalyze::earthSignal(TString s, TString XX, double cmunu){

    double pas = 1000;
    double tmp = M_PI/pas;
    double tmp2 = 2*(M_PI/pas);
    double b1,b2,b3,b4,b5,arg1,arg2,b1b2;
    TCanvas* c = new TCanvas("max f_{SME}(#lambda, #theta) "+s+" f"+XX,"Input"+s, 10,10,800,600);
    TH2F* h = new TH2F(s+XX, "", pas, 0, M_PI, pas, 0, 2*M_PI);
    for(int i = 0; i<pas; i++)
        for(int j=0; j<pas; j++){

            b1 = (sin(i*tmp)*sin(i*tmp)*sin(j*tmp2)*sin(j*tmp2) + cos(i*tmp)*cos(i*tmp)) * averageA(0,0) + (sin(i*tmp)*sin(i*tmp)*cos(j*tmp2)*cos(j*tmp2)) * averageA(2,2);
            b2 = cos(j*tmp2)*cos(j*tmp2) * averageA(0,0) + sin(j*tmp2)*sin(j*tmp2) * averageA(2,2);
            b1b2 = (b1-b2)/2;
            b3 = sin(i*tmp)*cos(j*tmp2)*sin(j*tmp2) * (averageA(2,2) - averageA(0,0));
            b4 = cos(i*tmp)*sin(i*tmp)*cos(j*tmp2)*cos(j*tmp2) * (averageA(2,2) - averageA(0,0));
            b5 = cos(j*tmp2)*cos(i*tmp2)*sin(j*tmp2) * (averageA(2,2) - averageA(0,0));
            arg1 = abs( atan(b3/b1b2)/(omega) );
            arg2 = abs( atan(b3/b1b2)/(2*omega) );

            if((b1-b2) != 0){
                if(XX == "XX" || XX == "YY")
                    h->SetBinContent(i+1, j+1, abs( 2*cmunu*(b1b2*cos(2*omega*arg2) + b3*sin(2*omega*arg2)) ) );
                else if(XX == "XY" || XX == "YX")
                    h->SetBinContent(i+1, j+1, abs( 2*cmunu*(b1b2*sin(2*omega*arg2) - b3*cos(2*omega*arg2)) ) );
                else if(XX == "XZ" || XX == "ZX")
                    h->SetBinContent(i+1, j+1, abs( 2*cmunu*(b4*cos(omega*arg1) + b5*sin(2*omega*arg1)) ) );
                else if(XX == "YZ" || XX == "ZY")
                    h->SetBinContent(i+1, j+1, abs( 2*cmunu*(b4*sin(omega*arg1) - b5*cos(2*omega*arg1)) ) );
                else 
                    cout<<"error TH2F"<<endl;
            }
        }

    h->SetTitle("max(f_{SME}(#lambda, #theta)) "+s+" f"+XX);
    h->GetYaxis()->SetTitle("Azimuth #theta (in rad)");
    h->GetXaxis()->SetTitle("Latitude #lambda (in rad)");

    h->Write();
    h->Draw("colz");
    h->SetStats(0);

    c->SaveAs("results/"+s+"maxf_{SME}(lambda,theta)"+XX+".png");
}

TH1F* FunctionAnalyze::useHisto(TString s, TString XX, double cmunu, int maxTime, int timeStep, Color_t color){
//    TCanvas* w = new TCanvas("f functions histogram"+XX+s,"",200,10,800,600);
    TH1F* foo;

    if(XX=="XX" || XX=="YY"){
        foo = new TH1F("c_{XX} = c_{YY} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSME("XX", cmunu, i));
        }
    }
    else if(XX=="XY" || XX=="YX"){
        foo = new TH1F("c_{XY} = c_{YX} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSME("XY", cmunu, i));
        }
    }
    else if(XX=="XZ" || XX=="ZX"){
        foo = new TH1F("c_{XZ} = c_{ZX} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSME("XZ", cmunu, i));
        }
    }
    else if(XX=="YZ" || XX=="ZY"){
        foo = new TH1F("c_{YZ} = c_{ZY} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSME("YZ", cmunu, i));
        }
    }
    else
        cout<<"error with f fuction"<<endl;

    foo->SetLineWidth(2);   foo->SetLineColor(color);
    foo->SetStats(0);

    return foo;
}

TH1F* FunctionAnalyze::useHistoD0(TString s, TString XX, double cmunu, int maxTime, int timeStep, Color_t color){
//    TCanvas* w = new TCanvas("f functions histogram"+XX+s,"",200,10,800,600);
    TH1F* foo;

    if(XX=="XX" || XX=="YY"){
        foo = new TH1F("c_{XX} = c_{YY} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSMED0("XX", cmunu, i));
        }
    }
    else if(XX=="XY" || XX=="YX"){
        foo = new TH1F("c_{XY} = c_{YX} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSMED0("XY", cmunu, i));
        }
    }
    else if(XX=="XZ" || XX=="ZX"){
        foo = new TH1F("c_{XZ} = c_{ZX} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSMED0("XZ", cmunu, i));
        }
    }
    else if(XX=="YZ" || XX=="ZY"){
        foo = new TH1F("c_{YZ} = c_{ZY} "+s,"", timeStep, 0, maxTime);
        for (int i=0; i<timeStep; i++) {
            foo->SetBinContent(i+1, fSMED0("YZ", cmunu, i));
        }
    }
    else
        cout<<"error with f fuction"<<endl;

    foo->SetLineWidth(2);   foo->SetLineColor(color);
    foo->SetStats(0);

    return foo;
}

//________________ statics functions _________________//

TH1F* FunctionAnalyze::useConstHisto(TString s, double signal, int maxTime, int nbin){
    TH1F* foo = new TH1F(s,s,nbin,0,maxTime);
    for (int i=0; i<nbin; i++) 
        foo->SetBinContent(i+1,signal/((double)nbin));
    return foo;
}

void FunctionAnalyze::amplEnergy(TString s, TMatrixD m1, TMatrixD m2, TMatrixD m3, TMatrixD mT){
    TCanvas* w = new TCanvas("Amplitude = f(Energy)"+s,"",200,10,800,600);
    TH1F* h1 = new TH1F("A^{ZZ}"+s,"", 15, 0, 15);
    TH1F* h2 = new TH1F("A^{XX}"+s,"", 15, 0, 15);

    h1->SetBinContent(1, -mT(2,2));
    h1->SetBinContent(2,-m1(2,2));
    h1->SetBinContent(7,-m2(2,2));
    h1->SetBinContent(13,-m3(2,2));
    h2->SetBinContent(1, -mT(0,0));
    h2->SetBinContent(2,-m1(0,0));
    h2->SetBinContent(7,-m2(0,0));
    h2->SetBinContent(13,-m3(0,0));

    h1->Draw("");
    h2->Draw("SAME");

    h1->GetYaxis()->SetTitle("A^{#mu #nu}");
    h1->GetXaxis()->SetTitle("Energy beam (in TeV)");

    h1->SetStats(0);

    h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);
    h2->Write();   h2->SetLineWidth(2);   h2->SetLineColor(kBlue);

    w->Update();
    w->BuildLegend();

    w->SaveAs("results/"+s+"AmatricesComparaison.png");
}

void FunctionAnalyze::fatHisto(TString s, TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4){
    TCanvas* w = new TCanvas("fatHisto"+s,"",200,10,800,600);
    h1->Draw();
    h2->Draw("SAME");
    h3->Draw("SAME");
    h4->Draw("SAME");

    h1->SetStats(0);

    w->Update();

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram "+s,"C"); // option "C" allows to center the header
    legend->AddEntry(h1,"LHC 13TeV ","l"); // option "l" is for line (form of legend)
    legend->AddEntry(h2,"LHC 7TeV ","l");
    legend->AddEntry(h3,"LHC 2TeV ","l");
    legend->AddEntry(h4,"Tevatron ","l");
    legend->Draw();

    w->SaveAs("results/"+s+"energyComparaison.png");
}

void FunctionAnalyze::fatHistoSwitch(TString s, TH1F* h1, TH1F* h2){
    TCanvas* w = new TCanvas("fatHisto"+s,"",200,10,800,600);
    h1->Draw();
    h2->Draw("SAME");

    h1->SetStats(0);
    h1->SetMaximum(0.02);
    h1->SetMinimum(-0.02);

    w->Update();

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram switch angles "+s,"C"); // option "C" allows to center the header
    legend->AddEntry(h1,"LHC 2TeV "+s,"l");
    legend->AddEntry(h2,"Tevatron "+s,"l");
    legend->Draw();

    w->SaveAs("results/"+s+"energyComparaison.png");
}

