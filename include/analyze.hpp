#ifndef analyze_h
#define analyze_h
#define analyze_cxx

#include <TString.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TH1F.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <cmath>

#ifndef __CINT__

using namespace std;

void showMatrix(TMatrixD m);

//**********************************************************************************//
//****************************** class Analyze declaration *************************//
//**********************************************************************************//

class EventsAnalyze {

    private :

        int nEvents;

// number of events and average Amunu matrices
        int nF;
        int nPqqbar;
        int nP2g;
        TMatrixD averageF;
        TMatrixD averagePqqbar;
        TMatrixD averageP2g;

// calculations
        int calculateEventNumber(string VM);
        TMatrixD calculateAverageMatrix(string VM);

    public :

// constructor
        EventsAnalyze();
        EventsAnalyze(int nEvents_user, string VMAPqqbar_user, string VMAP2g_user, string VMAF_user);
        ~EventsAnalyze(){};

// publics values
        int getnF();
        int getnPqqbar();
        int getnP2g();
        TMatrixD getAverageF();
        TMatrixD getAveragePqqbar();
        TMatrixD getAverageP2g();

// text file output 
        static void writeNumMatrix(string s, int n, TMatrixD m);

};

class FunctionAnalyze {

    private :

// average Amunu matrices A
        int nF;
        int nPqqbar;
        int nP2g;
        TMatrixD averageF;
        TMatrixD averagePqqbar;
        TMatrixD averageP2g;
        TMatrixD averageA;
        TMatrixD averageAD0;

// a coefficients of the fSME function
        double a1;
        double a2;
        double a3;
        double a4;
        double a5;
        double a1D0;
        double a2D0;
        double a3D0;
        double a4D0;
        double a5D0;

// calculations
        TMatrixD calculateAverageMatrixA(TMatrixD Pqqbar, TMatrixD P2g, TMatrixD F, int nP1, int nP2, int n);
        TMatrixD calculateAverageMatrixAD0(TMatrixD Pqqbar, TMatrixD F);
        double calculateCoefficent(int i, TMatrixD m);
        double calculateCoefficentD0(int i, TMatrixD m);

    public :

// constructor
        FunctionAnalyze();
        FunctionAnalyze(int nPqqbar_user, int nP2g_user, int nF_user, TMatrixD averagePqqbar_user, TMatrixD averageP2g_user, TMatrixD averageF_user);
        ~FunctionAnalyze(){};

// Text file input
        static int readNumber(string s);
        static TMatrixD readMatrix(string s);

// publics values
        int getnF();
        int getnPqqbar();
        int getnP2g();
        TMatrixD getAverageF();
        TMatrixD getAveragePqqbar();
        TMatrixD getAverageP2g();
        TMatrixD getAverageA();
        TMatrixD getAverageAD0();

        double geta1();
        double geta2();
        double geta3();
        double geta4();
        double geta5();
        double geta1D0();
        double geta2D0();
        double geta3D0();
        double geta4D0();
        double geta5D0();

// fSME function
        double fSME(string s, double cmunu, double t);
        double fSMED0(string s, double cmunu, double t);

// histograms
        void quadriHisto(TString s, double cmunu, int maxTime, int timeStep);
        void signalView(TString s, TString XX, double cmunu, int nbin, double nEventBkgd, double nEventTTbarSM);
        void signalHisto(TString s, TString XX, double cmunu, int nbin, double nEventTTbarSM);
        void earthSignal(TString s, TString XX, double cmunu);

// To compare signal with energy and use for stat
        TH1F* useHisto(TString s, TString XX, double cmunu, int maxTime, int timeStep, Color_t color);
        TH1F* useHistoD0(TString s, TString XX, double cmunu, int maxTime, int timeStep, Color_t color);
        static TH1F* useConstHisto(TString s, double signal, int maxTime, int nbin);

// Energy comparaison
        static void amplEnergy(TString s, TMatrixD m1, TMatrixD m2, TMatrixD m3, TMatrixD mT);
        static void fatHisto(TString s, TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4);
        static void fatHistoSwitch(TString s, TH1F* h1, TH1F* h2);
};


#endif
#endif
