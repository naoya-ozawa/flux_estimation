#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

double T_n_l (double KE, int A, double l){
    // KE of emitted n (MeV), A = mass # of compound nucleus
    // l: orbital angular momentum of compound nucleus
    double u = 931.5; // MeV/c^2
    double m_n = 1.0 * u; // MeV/c^2
    double M = double(A) * u; // MeV/c^2; mass of CN
    double RM = (M-m_n)*m_n / M; // MeV/c^2
    double hbarc = 197.0; // MeV fm
    double c = 3.0 * TMath::Power(10.,8); // fm/fs
    double Pf = TMath::Sqrt(KE / (KE + 40)); 
    double P = (4.0 * Pf) / ((1.0 + Pf) * (1.0 + Pf));
    double rho = TMath::Sqrt(l*(l+1) / (2.0 * RM * KE)) * hbarc; // fm
    double R = 1.2 * TMath::Power(double(A),1./3.) + 2.0; // fm
    double sqdiff = TMath::Abs(rho*rho - R*R);
    double coef = -2.0 * TMath::Sqrt(2.0 * RM * KE); // MeV/c
    double inlog = (rho + sqdiff) / R;
    return P * TMath::Exp( coef * ( sqdiff*sqdiff + rho*TMath::Log(inlog) ) );
}

double T_n (double *x, double *par){
    double e = x[0];
    int A = int(par[0]);
    int N = int(par[1]); // truncate summation at N-1
    double prob = 0.0;
    for (int l=0; l<N; ++l){
        prob += T_n_l(e,A,double(l));
    }
    return 100. * prob;
}

int main (int argc, char** argv){
    TRint rootapp("app",&argc,argv);
    TCanvas *c1 = new TCanvas();

    TF1 *neutron = new TF1("neutron",T_n,TMath::Power(10.,-9),TMath::Power(10.,2),2);
    neutron->SetParameters(215.,100.);
    neutron->SetTitle("Energy distribution of evaporated neutrons;Kinetic energy of neutron (MeV);Probability (%)");
    neutron->Draw();

    c1->Update();
    c1->Modified();

    rootapp.Run();

    return 0;
}
