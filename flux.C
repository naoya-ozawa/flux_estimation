void flux(){
    double distance = 10.0; // m

    TCanvas *c1 = new TCanvas();

    const char* data = "flux.dat";

    TGraphErrors *cyric = new TGraphErrors(data,"%lg %lg %lg");
    cyric->SetTitle("CYRIC (0.25 p#muA)");
    cyric->SetMarkerColor(6);
    cyric->SetMarkerStyle(11);


    TF1 *cyric_fit = new TF1("cyric_fit","[0]/(x*x) + [1]",0.5,8.0);
    cyric_fit->SetLineColor(6);
    cyric_fit->SetParameters(10.0,0.0);
    cyric->Fit(cyric_fit,"M","",0.5,8.0);

    // CYRIC data: 1.5 uA = 0.3 puA (18-O 5+ ions)
    // CNS data: 3 puA = 18 uA (18-O 6+ ions)
    double conversion = 3.0/0.3;

    TGraphErrors *cns = new TGraphErrors();
    for (int i =0; i < cyric->GetN(); ++i){
        cns->SetPoint(i,cyric->GetX()[i],conversion*cyric->GetY()[i]);
        cns->SetPointError(i,0.0,conversion*cyric->GetEY()[i]);
    }
    cns->SetTitle("CNS (3.0 p#muA)");
    cns->SetMarkerColor(4);
    cns->SetMarkerStyle(22);

    TF1 *cns_estimate = new TF1("cns_estimate","[0]/(x*x) + [1]",0.5,8.0);
    cns_estimate->SetLineColor(4);
    cns_estimate->SetParameters(conversion*cyric_fit->GetParameter(0),conversion*cyric_fit->GetParameter(1));

    TMultiGraph *comp = new TMultiGraph();
    comp->Add(cyric);
    comp->Add(cns);
    comp->SetTitle("Neutron Radiation; Distance (m); Radiation (Bq)");

    comp->Draw("AP*");

    comp->GetYaxis()->SetRange(0.0,300.0);

    gPad->BuildLegend();

    cns_estimate->Draw("SAME");


    c1->Update();
    c1->Modified();
   
    conversion = 1.0;
 
    double a0 = conversion*cyric_fit->GetParameter(0);
    double a0E = conversion*cyric_fit->GetParError(0);
    double a1 = conversion*cyric_fit->GetParameter(1);
    double a1E = conversion*cyric_fit->GetParError(1);
    double half_life = 2.7*24.0; // hours
    double lifetime = half_life/TMath::Log(2.0); // hours
    double time_elapsed = 3.0; // hours

    // 197Au + n -> 198Au e.s. -> 198Au g.s. + gamma (411keV)

    cout << "===========================================================" << endl;
    cout << "At CYRIC, " << distance << " m from the target:" << endl;

    cout << "Lifetime of Au-198 = " << lifetime << " hours." << endl;

    double gamma_radiation = a0/(distance*distance) + a1; // Bq
    double gamma_radiationE = TMath::Sqrt( (a0E*a0E)/(distance*distance*distance*distance) + (a1E*a1E) ); // Bq

    cout << "Expected gamma radiation of Au-198: " << gamma_radiation << " +- " << gamma_radiationE << " Bq" << endl;

    double n_0 = lifetime * gamma_radiation * TMath::Exp(time_elapsed/lifetime);
    double n_0E = lifetime * gamma_radiationE * TMath::Exp(time_elapsed/lifetime);

    cout << "Estimated number N(0) of Au-198 produced on film by n-irradiation: " << n_0 << " +- " << n_0E << endl;

    double cross_section_barn = 100.0; // b at energy En
    double cross_section = cross_section_barn*TMath::Power(10.,-28); // m^2

    cout << "===========================================================" << endl;
    cout << "Assume monochromatic neutron at energy XX; cross section = " << cross_section << " m^2" << endl;

    double au_density = 19.32*TMath::Power(10.,6); // m^-3
    double au_thickness = 0.2*TMath::Power(10.,-3); // m
    double n_on_film = n_0/(cross_section*au_density*au_thickness);
    double n_on_filmE = n_0E/(cross_section*au_density*au_thickness);

    cout << "Estimated number of neutrons that penetrated the film in total: " << n_on_film << " +- " << n_on_filmE << endl;

    double irradiation_time = 12.0*60.0*60.0; // seconds

    cout << "Neutron flux at distance " << distance << " m = " << n_on_film/irradiation_time << " +- " << n_on_filmE/irradiation_time << " /s/cm^2" << endl;


}
