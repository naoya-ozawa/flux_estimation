void flux(){
    double distance = 5.0; // m

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

    double half_life = 2.7*24.0; // hours
    double lifetime = half_life/TMath::Log(2.0); // hours
    double irradiation_time = 12.0; // hours
    double time_elapsed = 3.0; // hours

    // 197Au + n -> 198Au e.s. -> 198Au g.s. + gamma (411keV)

    cout << "===========================================================" << endl;
    cout << "At CYRIC, " << distance << " m from the target:" << endl;

    cout << "Lifetime of Au-198 = " << lifetime << " hours." << endl;

    double gamma_radiation = a0/(distance*distance); // Bq
    double gamma_radiationE = a0E/(distance*distance); // Bq

    cout << "Expected gamma radiation of Au-198: " << gamma_radiation << " +- " << gamma_radiationE << " Bq" << endl;

    cout << "===========================================================" << endl;
    cout << "Assuming a monochromatic neutron of energy En = 1MeV" << endl;

    double cross_section_b = 0.1; // b
    double cross_section = cross_section_b*TMath::Power(10.,-28);

    double time_factor = ( 1.0 - TMath::Exp(-(irradiation_time/lifetime)) ) * TMath::Exp(-(time_elapsed/lifetime));

    double n_flux = gamma_radiation/(cross_section*time_factor); // /s/m^2
    double n_fluxE = gamma_radiationE/(cross_section*time_factor); // /s/m^2

    cout << "Neutron flux at distance " << distance << " m: " << n_flux << " +- " << n_fluxE << " /s/m^2" << endl;

    cout << "Number of neutrons hitting 197-Au: " << cross_section*n_flux << " +- " << cross_section*n_fluxE << endl;


}
