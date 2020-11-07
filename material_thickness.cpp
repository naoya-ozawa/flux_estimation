#include <iostream>
#include <TVectorD.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TSpline.h>
#include <TLatex.h>
#include <TLine.h>
#include <TRint.h>
using namespace std;

int main (int argc, char** argv){
    TRint rootapp("app",&argc,argv);

    TCanvas *c1 = new TCanvas("c1");
    c1->Divide(2,1);

    c1->cd(1);

    const char* filename = "./X4sShowData_B10ABS.txt";

    TTree *tree = new TTree("tree","neutron absorption cross section data");
    double e_n, e_n_err, cs, cs_err;
    tree->Branch("e_n",&e_n,"e_n/D"); // MeV
    tree->Branch("e_n_err",&e_n_err,"e_n_err/D"); // MeV
    tree->Branch("cs",&cs,"cs/D"); // b
    tree->Branch("cs_err",&cs_err,"cs_err/D"); // b
    tree->ReadFile(filename,"e_n/D:e_n_err/D:cs/D:cs_err/D");

    tree->Project("","e_n:e_n_err:cs:cs_err");
    TVectorD *e_n_vec = new TVectorD(tree->GetSelectedRows(),tree->GetV1());
    TVectorD *e_n_err_vec = new TVectorD(tree->GetSelectedRows(),tree->GetV2());
    TVectorD *cs_vec = new TVectorD(tree->GetSelectedRows(),tree->GetV3());
    TVectorD *cs_err_vec = new TVectorD(tree->GetSelectedRows(),tree->GetV4());
    TGraphErrors *ge_cs = new TGraphErrors(*e_n_vec,*cs_vec,*e_n_err_vec,*cs_err_vec);
    ge_cs->SetTitle("Neutron absorption cross section data from EXFOR DB;Neutron energy E_{n} (MeV);Cross section #sigma (b)");
    ge_cs->Draw("AP*");
    TF1 *f_cs = new TF1("f_cs","[0]*TMath::Exp(-[1]*x**[2])/TMath::Sqrt(x)");
    f_cs->SetParameters(0.8,15.,5.0);
    ge_cs->Fit("f_cs","ME","",0.,1.);
//    f_cs->Draw("SAME");

    c1->cd(1)->SetLogx();
    c1->cd(1)->SetLogy();


    c1->cd(2);

    double P_ing = 0.5*(4./5.); // 50% B4C
    double rho_ing = 2.08; // g/cm3
    double E_isot = 0.2; // B-10 abundance
    double NA = 6.*TMath::Power(10.,23); // /mol
    double A = 10.; // g/mol
    double N_isot = P_ing*rho_ing*E_isot*NA/A; // /cm3

    double E_neutron = 0.025*TMath::Power(10.,-6); // MeV, thermal
    double cs_isot = f_cs->Eval(E_neutron);

    TGraph *g_shield = new TGraph(1000);
    for (int i=0; i<1000; ++i){
        double factor = 100.*(1. - TMath::Exp(-cs_isot*TMath::Power(10.,-24)*N_isot*double(i+1)/1000.));
        g_shield->SetPoint(i,(i+1)/1000.,factor);
    }
    g_shield->SetTitle("Shielding factor for thermal neutrons;Shield thickness (cm);Shielding factor (%)");
    g_shield->Draw("AL");

    c1->cd(2)->SetLogy();


    c1->Update();
    c1->Modified();

    rootapp.Run();
    return 0;
}