
#include "../include/Const.hpp"
#include "../include/analyze.hpp"

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1.h>
#include <TTree.h>

#include <iostream>
#include <fstream>

using namespace std;

# ifndef __CINT__

int main(){

//_______________ To keep application open until quit ROOT _______________//

    TApplication app("app", nullptr, nullptr);

//_____________________ All needed parameters _____________________//

    int nEvents13TeV = 5000000;
    int nEvents7TeV  = 1000000; //same as 2Tev and Tevatron

    double lumino    = 36;
    double Bkgd      = (646/2.2) * lumino;
    double TTbarSM   = (9921/2.2) * lumino;

    double cmunu = 0.01;
    int time = 24;              //hours
    int bin = 24;               //bins
    int step = 84600;           // 1 events per second

//_______________________ Condition program ______________________//

    string activation;
    do{
        cout<<"Enter a commande ('creation' or 'analyze') : ";
        cin>>activation;
    }
    while(activation != "analyze" && activation != "creation");
    cout<<endl;

//_________________ Write matrices and number of Events in text files ________________//

    if(activation=="creation"){

        EventsAnalyze LHC13TeVCMS(nEvents13TeV, "data/LHC13TeV/VMAPqqbarCMS.txt", "data/LHC13TeV/VMAP2gCMS.txt", "data/LHC13TeV/VMAFCMS.txt");
        EventsAnalyze LHC7TeV(nEvents7TeV, "data/LHC2-7TeV/VMAPqqbar7TeVCMS.txt", "data/LHC2-7TeV/VMAP2g7TeVCMS.txt", "data/LHC2-7TeV/VMAF7TeVCMS.txt");
        EventsAnalyze LHC2TeV(nEvents7TeV, "data/LHC2-7TeV/VMAPqqbar2TeVCMS.txt", "data/LHC2-7TeV/VMAP2g2TeVCMS.txt", "data/LHC2-7TeV/VMAF2TeVCMS.txt");
        EventsAnalyze TEV(nEvents7TeV, "data/Tevatron/VMAPqqbarTevatron.txt", "data/Tevatron/VMAP2gTevatron.txt", "data/Tevatron/VMAFTevatron.txt");

        cout<<"LHC 13 TeV"<<endl;
        cout<<endl;

        EventsAnalyze::writeNumMatrix("13TeVCMSF", LHC13TeVCMS.getnF(), LHC13TeVCMS.getAverageF());
        EventsAnalyze::writeNumMatrix("13TeVCMSPqqbar", LHC13TeVCMS.getnPqqbar(), LHC13TeVCMS.getAveragePqqbar());
        EventsAnalyze::writeNumMatrix("13TeVCMSP2g", LHC13TeVCMS.getnP2g(), LHC13TeVCMS.getAverageP2g());

        cout<<"LHC 7 TeV"<<endl;
        cout<<endl;

        EventsAnalyze::writeNumMatrix("7TeVCMSF", LHC7TeV.getnF(), LHC7TeV.getAverageF());
        EventsAnalyze::writeNumMatrix("7TeVCMSPqqbar", LHC7TeV.getnPqqbar(), LHC7TeV.getAveragePqqbar());
        EventsAnalyze::writeNumMatrix("7TeVCMSP2g", LHC7TeV.getnP2g(), LHC7TeV.getAverageP2g());

        cout<<"LHC 2 TeV"<<endl;
        cout<<endl;

        EventsAnalyze::writeNumMatrix("2TeVCMSF", LHC2TeV.getnF(), LHC2TeV.getAverageF());
        EventsAnalyze::writeNumMatrix("2TeVCMSPqqbar", LHC2TeV.getnPqqbar(), LHC2TeV.getAveragePqqbar());
        EventsAnalyze::writeNumMatrix("2TeVCMSP2g", LHC2TeV.getnP2g(), LHC2TeV.getAverageP2g());

        cout<<"Tevatron 2 TeV"<<endl;
        cout<<endl;

        EventsAnalyze::writeNumMatrix("TEVF", TEV.getnF(), TEV.getAverageF());
        EventsAnalyze::writeNumMatrix("TEVPqqbar", TEV.getnPqqbar(), TEV.getAveragePqqbar());
        EventsAnalyze::writeNumMatrix("TEVP2g", TEV.getnP2g(), TEV.getAverageP2g());

        cout<<"End Bye Bro !"<<endl;
        return 0;
    }

//_________________ Analyze the Lorentz violation ________________//

    else if(activation=="analyze"){

        TFile* f = new TFile("AnalyzeTTbar.root", "RECREATE");

        //LHC 13 TeV

        int nF              = FunctionAnalyze::readNumber("13TeVCMSF");
        int nPqqbar         = FunctionAnalyze::readNumber("13TeVCMSPqqbar");
        int nP2g            = FunctionAnalyze::readNumber("13TeVCMSP2g");

        TMatrixD F          = FunctionAnalyze::readMatrix("13TeVCMSF");
        TMatrixD Pqqbar     = FunctionAnalyze::readMatrix("13TeVCMSPqqbar");
        TMatrixD P2g        = FunctionAnalyze::readMatrix("13TeVCMSP2g");

        FunctionAnalyze LHC13TeV(nPqqbar, nP2g, nF, Pqqbar, P2g, F);
        LHC13TeV.quadriHisto("LHC 13Tev", cmunu, time, step);
        LHC13TeV.signalView("LHC 13TeV", "XX", cmunu, bin, Bkgd, TTbarSM);
        LHC13TeV.signalView("LHC 13TeV", "XZ", cmunu, bin, Bkgd, TTbarSM);
        TH1F* h13XX = LHC13TeV.useHisto("LHC 13TeV", "XX", cmunu, time, step, kRed);
        TH1F* h13XY = LHC13TeV.useHisto("LHC 13TeV", "XY", cmunu, time, step, kRed);
        TH1F* h13XZ = LHC13TeV.useHisto("LHC 13TeV", "XZ", cmunu, time, step, kRed);
        TH1F* h13YZ = LHC13TeV.useHisto("LHC 13TeV", "YZ", cmunu, time, step, kRed);

        //LHC 7 TeV

        nF                  = FunctionAnalyze::readNumber("7TeVCMSF");
        nPqqbar             = FunctionAnalyze::readNumber("7TeVCMSPqqbar");
        nP2g                = FunctionAnalyze::readNumber("7TeVCMSP2g");

        F                   = FunctionAnalyze::readMatrix("7TeVCMSF");
        Pqqbar              = FunctionAnalyze::readMatrix("7TeVCMSPqqbar");
        P2g                 = FunctionAnalyze::readMatrix("7TeVCMSP2g");

        FunctionAnalyze LHC7TeV(nPqqbar, nP2g, nF, Pqqbar, P2g, F);
        LHC7TeV.quadriHisto("LHC 7Tev", cmunu, time, step);
        TH1F* h7XX = LHC7TeV.useHisto("LHC 7TeV", "XX", cmunu, time, step, kBlue);
        TH1F* h7XY = LHC7TeV.useHisto("LHC 7TeV", "XY", cmunu, time, step, kBlue);
        TH1F* h7XZ = LHC7TeV.useHisto("LHC 7TeV", "XZ", cmunu, time, step, kBlue);
        TH1F* h7YZ = LHC7TeV.useHisto("LHC 7TeV", "YZ", cmunu, time, step, kBlue);

        //LHC 2 TeV

        nF                  = FunctionAnalyze::readNumber("2TeVCMSF");
        nPqqbar             = FunctionAnalyze::readNumber("2TeVCMSPqqbar");
        nP2g                = FunctionAnalyze::readNumber("2TeVCMSP2g");

        F                   = FunctionAnalyze::readMatrix("2TeVCMSF");
        Pqqbar              = FunctionAnalyze::readMatrix("2TeVCMSPqqbar");
        P2g                 = FunctionAnalyze::readMatrix("2TeVCMSP2g");

        FunctionAnalyze LHC2TeV(nPqqbar, nP2g, nF, Pqqbar, P2g, F);
        LHC2TeV.signalHisto("LHC 2TeV", "XX", cmunu, bin, TTbarSM);
        TH1F* h2XX = LHC2TeV.useHisto("LHC 2TeV", "XX", cmunu, time, step, kMagenta);
        TH1F* h2XY = LHC2TeV.useHisto("LHC 2TeV", "XY", cmunu, time, step, kMagenta);
        TH1F* h2XZ = LHC2TeV.useHisto("LHC 2TeV", "XZ", cmunu, time, step, kMagenta);
        TH1F* h2YZ = LHC2TeV.useHisto("LHC 2TeV", "YZ", cmunu, time, step, kMagenta);

        //Tevatron 2 TeV

        nF                  = FunctionAnalyze::readNumber("2TeVCMSF");
        nPqqbar             = FunctionAnalyze::readNumber("2TeVCMSPqqbar");
        nP2g                = FunctionAnalyze::readNumber("2TeVCMSP2g");

        F                   = FunctionAnalyze::readMatrix("TEVF");
        Pqqbar              = FunctionAnalyze::readMatrix("TEVPqqbar");
        P2g                 = FunctionAnalyze::readMatrix("TEVP2g");

        FunctionAnalyze TEV(nPqqbar, nP2g, nF, Pqqbar, P2g, F);
        TEV.quadriHisto("LHC 7Tev", cmunu, time, step);
        TH1F* hTEVXX = TEV.useHistoD0("TEV", "XX", cmunu, time, step, kGreen);
        TH1F* hTEVXY = TEV.useHistoD0("TEV", "XY", cmunu, time, step, kGreen);
        TH1F* hTEVXZ = TEV.useHistoD0("TEV", "XZ", cmunu, time, step, kGreen);
        TH1F* hTEVYZ = TEV.useHistoD0("TEV", "YZ", cmunu, time, step, kGreen);

        FunctionAnalyze::fatHisto("f^{XX}", h13XX, h7XX, h2XX, hTEVXX);
        FunctionAnalyze::fatHisto("f^{XY}", h13XY, h7XY, h2XY, hTEVXY);
        FunctionAnalyze::fatHisto("f^{XZ}", h13XZ, h7XZ, h2XZ, hTEVXZ);
        FunctionAnalyze::fatHisto("f^{YZ}", h13YZ, h7YZ, h2YZ, hTEVYZ);

        FunctionAnalyze::amplEnergy(" decay", LHC2TeV.getAverageF(), LHC7TeV.getAverageF(), LHC13TeV.getAverageF(), TEV.getAverageF());
        FunctionAnalyze::amplEnergy(" q#bar{q}", LHC2TeV.getAveragePqqbar(), LHC7TeV.getAveragePqqbar(), LHC13TeV.getAveragePqqbar(), TEV.getAveragePqqbar());
        FunctionAnalyze::amplEnergy(" gluglu", LHC2TeV.getAverageP2g(), LHC7TeV.getAverageP2g(), LHC13TeV.getAverageP2g(), TEV.getAverageP2g());
        FunctionAnalyze::amplEnergy("", LHC2TeV.getAverageA(), LHC7TeV.getAverageA(), LHC13TeV.getAverageA(), TEV.getAverageAD0());

        cout<<endl;
        cout<<"LHC13Tev a1-a2 : "<<LHC13TeV.geta1()-LHC13TeV.geta2()<<endl;
        cout<<"LHC7Tev  a1-a2 : "<<LHC7TeV.geta1()-LHC7TeV.geta2()<<endl;
        cout<<"LHC2Tev  a1-a2 : "<<LHC2TeV.geta1()-LHC2TeV.geta2()<<endl;
        cout<<"TEV      a1-a2 : "<<TEV.geta1D0()-TEV.geta2D0()<<endl;
        cout<<endl;
        cout<<"LHC13Tev a3 : "<<LHC13TeV.geta3()<<endl;
        cout<<"LHC7Tev  a3 : "<<LHC7TeV.geta3()<<endl;
        cout<<"LHC2Tev  a3 : "<<LHC2TeV.geta3()<<endl;
        cout<<"TEV      a3 : "<<TEV.geta3D0()<<endl;
        cout<<endl;
        cout<<"LHC13Tev a4 : "<<LHC13TeV.geta4()<<endl;
        cout<<"LHC7Tev  a4 : "<<LHC7TeV.geta4()<<endl;
        cout<<"LHC2Tev  a4 : "<<LHC2TeV.geta4()<<endl;
        cout<<"TEV      a4 : "<<TEV.geta4D0()<<endl;
        cout<<endl;
        cout<<"LHC13Tev a5 : "<<LHC13TeV.geta5()<<endl;
        cout<<"LHC7Tev  a5 : "<<LHC7TeV.geta5()<<endl;
        cout<<"LHC2Tev  a5 : "<<LHC2TeV.geta5()<<endl;
        cout<<"TEV      a5 : "<<TEV.geta5D0()<<endl;
        cout<<endl;

        TH1F* h2XZTEVangle = LHC2TeV.useHistoD0("LHC 2TeV TEV angles", "XZ", cmunu, time, step, kMagenta);
        TH1F* hTEVXZLHCangle = TEV.useHisto("TEV LHC angles", "XZ", cmunu, time, step, kGreen);
        FunctionAnalyze::fatHistoSwitch("Switch", h2XZTEVangle, hTEVXZLHCangle);
        FunctionAnalyze::fatHistoSwitch("Unswitch", h2XZ, hTEVXZ);

        LHC13TeV.earthSignal("13TeV", "XX", cmunu);
        LHC13TeV.earthSignal("13TeV", "XY", cmunu);
        LHC13TeV.earthSignal("13TeV", "XZ", cmunu);
        LHC13TeV.earthSignal("13TeV", "YZ", cmunu);
        LHC2TeV.earthSignal("2TeV", "XX", cmunu);
        LHC2TeV.earthSignal("2TeV", "XY", cmunu);
        LHC2TeV.earthSignal("2TeV", "XZ", cmunu);
        LHC2TeV.earthSignal("2TeV", "YZ", cmunu);

        //Over
        cout<<"End Bye Bro !"<<endl;
        app.Run();
        f->Write();
        f->Close();
        return 0;
    }
    else{
        cout<<"Warning !!! Problem with the code"<<endl;
        return 0;
    }

//Bye.
}

#endif

