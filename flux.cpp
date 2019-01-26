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

int main (int argc, char** argv){

    // Check input number
    if (argc != 4){
        cout << "usage: ./flux <O-18 intensity [puA]> <neutron energy [eV]> <distance [m]>" << endl;
        exit(1);
    }
    double I = strtod(argv[1],NULL); // puA
    double E_n = strtod(argv[2],NULL); // eV
    double distance = strtod(argv[3],NULL); // m

    cout << "======================================================" << endl;
    cout << "O-18(6+) Intensity = " << I << " puA" << endl;
    cout << "Assume neutrons of energy E_n = " << E_n << " eV" << endl;
    cout << "Calculate at distance = " << distance << " m" << endl;
    cout << "======================================================" << endl;

    // Constants
    double half_life = 2.7*24.; // hours
    double irradiation_time = 12.; // hours
    double time_elapsed = 3.; // hours

    double Au_ram = 197.; // g/mol
    double Na = 6.*TMath::Power(10.,23); // #/mol
    double Au_density = 19.32*TMath::Power(10.,6); // g/m^3
    double Au_area = 1.*1.*TMath::Power(10.,-4); // m^2

    // Constants used for calculation
    double T = 0.2*TMath::Power(10.,-3); // m
    double N_197Au = Au_density*Na/Au_ram; // #/m^3
    double N_0 = N_197Au*Au_area*T; // #
    double tau_l = half_life/TMath::Log(2.0); // hours <lifetime>
    double t_0 = irradiation_time; // hours
    double t_1 = t_0 + time_elapsed; // hours

    cout << "Gamma radiation data by CYRIC:" << endl;
    cout << "O-18(5+) Intensity = 0.3 puA" << endl;
    cout << "Irradiation time = 12 hours" << endl;
    cout << "Elapsed time after irradiation = 3 hours" << endl;
    cout << "======================================================" << endl;


    TRint rootapp("app",&argc,argv);
    TCanvas *c1 = new TCanvas();
    c1->Divide(2,2);


    c1->cd(1);
    // 1. Plot the data from CYRIC and estimate the values at the set intensity
    const char* data = "flux.dat";

    TGraphErrors *cyric = new TGraphErrors(data,"%lg %lg %lg");
    cyric->SetTitle("CYRIC (0.3 p#muA)");
    cyric->SetMarkerColor(6);
    cyric->SetMarkerStyle(11);

    TF1 *cyric_fit = new TF1("cyric_fit","[0]/(x*x) + [1]",0.5,8.0);
    cyric_fit->SetLineColor(6);
    cyric_fit->SetParameters(10.,0.);
    cyric->Fit(cyric_fit,"M","",0.5,8.0);

    double conversion = I/0.3; 

    TGraphErrors *cns = new TGraphErrors();
    for(int i=0; i<cyric->GetN(); ++i){
        cns->SetPoint(i,cyric->GetX()[i],conversion*cyric->GetY()[i]);
        cns->SetPointError(i,0.,conversion*cyric->GetEY()[i]);
    }
    cns->SetTitle(Form("CNS (%f p#muA)",I));
    cns->SetMarkerColor(4);
    cns->SetMarkerStyle(22);

    TF1 *cns_estimate = new TF1("cns_estimate","[0]/(x*x) + [1]",0.5,8.0);
    double a0 = conversion*cyric_fit->GetParameter(0);
    double a0E = conversion*cyric_fit->GetParError(0);
    double a1 = conversion*cyric_fit->GetParameter(1);
    cns_estimate->SetLineColor(4);
    cns_estimate->SetParameters(a0,a1);

    TMultiGraph *comp = new TMultiGraph();
    comp->Add(cyric);
    comp->Add(cns);
    comp->SetTitle("Radiation from ^{197}Au(n,#gamma)^{198}Au; Distance (m); Radiation (Bq)");

    comp->Draw("AP*");
    comp->GetYaxis()->SetRange(0.,300.);
    gPad->BuildLegend();
    cns_estimate->Draw("SAME");




    c1->cd(2);
    // 2. Retrieve the neutron cross section data

    double E[10] = { 1.338e-8, 1.213e-7, 1.001e-6, 8.039e-6, 9.576e-4, 5.0795e-3, 1.0042e-2, 1.0042e-1, 1.12, 7.6 };
    double Ee[10] = { 0., 0., 0., 0., 0., 2.35e-5, 4.65e-5, 0.0004625, 0.12, 0.05 };
    double C[10] = { 133.46, 47.139, 23.999, 15.323, 5.53, 0.658, 1.471, 0.301, 0.060441, 0.0008 };
    double Ce[10] = { 3.3168, 1.1516, 0.61878, 0.4052, 1.1613, 0.1066, 0.13974, 0.049063, 0.003403, 0.0005 };

    TGraphErrors *csdata = new TGraphErrors(10,E,C,Ee,Ce);
    csdata->SetTitle("Cross section #sigma (E_n) of ^{197}Au(n,#gamma)^{198}Au; Neutron Energy E_n (MeV); Cross Section (b)");
    csdata->Draw("AP*");    

    TF1 *csfit = new TF1("csfit","[0]*TMath::Exp([1]+[2]*x)*x**(-[3])");
    csdata->Fit(csfit,"M","",0.,20.);
    
    c1->cd(2)->SetLogx();
    c1->cd(2)->SetLogy();

    double cross_section_b = csfit->GetParameter(0)*TMath::Exp(csfit->GetParameter(1)+csfit->GetParameter(2)*E_n*TMath::Power(10.,-6)) / TMath::Power(E_n*TMath::Power(10.,-6), csfit->GetParameter(3));

    cout << "======================================================" << endl;
    cout << "Extrapolated cross section at energy E_n = " << E_n << " eV: " << cross_section_b << " b" << endl;


    // 3. Calculate the estimated neutron flux
    // 4. Simulate the gamma radiation that will be obtained

    c1->Update();
    c1->Modified();

    rootapp.Run();

    return 0;
}
