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



double N_ex(double *x, double *par){
    double par0 = par[0];
    double par1 = par[1];
    double par2 = par[2];
    double par3 = par[3];
    double par4 = par[4];
    double par5 = par[5];

    double num = 0.;
    if (x[0] <= par0){
        num = par1*(1.0 - TMath::Exp(-par2*x[0]));
    }else{
        num = par3*TMath::Exp(-(x[0]-par4)/par5);
    }
    return num;
}




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
    cns->SetTitle(Form("CNS (%g p#muA)",I));
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
    csdata->SetTitle("Cross section #sigma (E_{n}) of ^{197}Au(n,#gamma)^{198}Au; Neutron Energy E_{n} (MeV); Cross Section (b)");
    csdata->Draw("AP*");    

    TF1 *csfit = new TF1("csfit","[0]*TMath::Exp([1]+[2]*x)*x**(-[3])");
    csdata->Fit(csfit,"M","",0.,20.);
    
    c1->cd(2)->SetLogx();
    c1->cd(2)->SetLogy();

    double cross_section_b = csfit->GetParameter(0)*TMath::Exp(csfit->GetParameter(1)+csfit->GetParameter(2)*E_n*TMath::Power(10.,-6)) / TMath::Power(E_n*TMath::Power(10.,-6), csfit->GetParameter(3));

    cout << "======================================================" << endl;
    cout << "Extrapolated cross section at energy E_n = " << E_n << " eV: " << cross_section_b << " b" << endl;
    cout << "Data samples (10 points) taken from EXFOR database" << endl;
    cout << "======================================================" << endl;
 




    c1->cd(3);
    // 3. Calculate the estimated neutron flux

    double radiation = a0/(distance*distance); // Bq (ignore background)
    double radiationE = a0E/(distance*distance); // Bq

    cout << "Expected radiation = " << radiation << " +- " << radiationE << " Bq" << endl;

    double cross_section = cross_section_b*TMath::Power(10.,-28); // m^2

    double damp_factor = 1.0 - TMath::Exp(-cross_section*N_197Au*T);
    double time_factor = (1.0 - TMath::Exp(-t_0/tau_l)) * TMath::Exp(-(t_1-t_0)/tau_l);

    double alpha = radiation * tau_l*60.*60. / (N_0*damp_factor*time_factor);
    double alphaE = radiationE * tau_l*60.*60. / (N_0*damp_factor*time_factor);

    cout << "Parameter alpha = " << alpha << " +- " << alphaE << endl;

    double n_flux = alpha*N_0/(tau_l*60.*60.*damp_factor);
    double n_fluxE = alphaE*N_0/(tau_l*60.*60.*damp_factor);

    cout << "Expected neutron flux: " << n_flux << " +- " << n_fluxE << " /cm^2/s" << endl;

    cout << "Expected neutron capture in gold foil: " << n_flux*damp_factor << " +- " << n_fluxE*damp_factor << " /cm^2/s" << endl;
    cout << "Use the appropriate cross section values to estimate for other materials" << endl;
    cout << "======================================================" << endl;

    TLatex *l = new TLatex();
    l->SetTextAlign(12);
    l->SetTextSize(0.05);
    l->DrawLatex(0.1,0.8,"Calculation condition:");
    l->DrawLatex(0.2,0.7,Form("Beam intensity: %g p#muA",I));
    l->DrawLatex(0.2,0.6,Form("Neutron energy: %g eV",E_n));
    l->DrawLatex(0.2,0.5,Form("Distance: %g m",distance));
    l->DrawLatex(0.1,0.3,"Estimated neutron flux:");
    l->DrawLatex(0.2,0.2,Form("Expected neutron flux: %g+-%g /cm^{2}/s",n_flux,n_fluxE));
    l->DrawLatex(0.2,0.1,Form("Expected neutron capture in gold foil: %g+-%g /cm^{2}/s",n_flux*damp_factor,n_fluxE*damp_factor));




    c1->cd(4);
    // 4. Simulate the gamma radiation that will be obtained

    TF1 *Nex = new TF1("Nex",N_ex,0.,120.,6);
    Nex->SetParameters(t_0,N_0*alpha,1.0/((1.0-alpha)*tau_l),N_0*alpha*(1.0-TMath::Exp(-t_0/((1.0-alpha)*tau_l))),t_0,tau_l);

    Nex->SetTitle("Expected number of excited ^{198}Au nuclei");
    Nex->Draw();
    Nex->GetXaxis()->SetTitle("Time (h)");
    Nex->GetYaxis()->SetTitle("N_{ex}(t)");



    c1->Update();
    c1->Modified();

    rootapp.Run();

    return 0;
}
