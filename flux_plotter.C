void flux_plotter(){

  double distance[3] = {2.,5.,10.}; // m
  double distance_err[3] = {0.,0.,0.};

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);

  // slow neutrons (0.025 eV)
  c1->cd(1);

//  // soft beam (0.02 puA)
//  double slow_soft[3] = {8.9,1.42,0.355};
//  double slow_soft_err[3] = {1.1,0.17,0.042};
//  TGraphErrors *ge_slow_soft = new TGraphErrors(3,distance,slow_soft,distance_err,slow_soft_err);
//  ge_slow_soft->SetTitle("0.02 p#muA");
//  ge_slow_soft->SetLineColor(2);

  // normal beam (0.7 puA)
  double slow_norm[3] = {3.11e+2,49.8,12.4};
  double slow_norm_err[3] = {3.7e+1,5.9,1.5};
  TGraphErrors *ge_slow_norm = new TGraphErrors(3,distance,slow_norm,distance_err,slow_norm_err);
  ge_slow_norm->SetTitle("0.7 p#muA");
  ge_slow_norm->SetLineColor(2);
  ge_slow_norm->SetLineWidth(2);

  // intense beam (3.0 puA)
  double slow_inte[3] = {1.33e+3,2.13e+2,53.3};
  double slow_inte_err[3] = {1.6e+2,2.5e+1,6.3};
  TGraphErrors *ge_slow_inte = new TGraphErrors(3,distance,slow_inte,distance_err,slow_inte_err);
  ge_slow_inte->SetTitle("3.0 p#muA");
  ge_slow_inte->SetLineColor(4);
  ge_slow_inte->SetLineWidth(2);

  TMultiGraph *mg_slow = new TMultiGraph();
  mg_slow->SetTitle("Estimated Thermal Neutron (0.025 eV) Flux;Distance from Ion Source (m);Flux (/(cm^{2} s))");
  mg_slow->Add(ge_slow_norm);
  mg_slow->Add(ge_slow_inte);
  mg_slow->Draw("ALP");
  c1->cd(1)->BuildLegend();

  // fast neutrons (1 MeV)
  c1->cd(2);

  // normal beam (0.7 puA)
  double fast_norm[3] = {3.93e+5,6.29e+4,1.57e+4};
  double fast_norm_err[3] = {4.6e+4,7.4e+3,1.9e+3};
  TGraphErrors *ge_fast_norm = new TGraphErrors(3,distance,fast_norm,distance_err,fast_norm_err);
  ge_fast_norm->SetTitle("0.7 p#muA");
  ge_fast_norm->SetLineColor(2);
  ge_fast_norm->SetLineWidth(2);

  // intense beam (3.0 puA)
  double fast_inte[3] = {1.69e+6,2.70e+5,6.74e+4};
  double fast_inte_err[3] = {2.0e+5,3.2e+4,8.0e+3};
  TGraphErrors *ge_fast_inte = new TGraphErrors(3,distance,fast_inte,distance_err,fast_inte_err);
  ge_fast_inte->SetTitle("3.0 p#muA");
  ge_fast_inte->SetLineColor(4);
  ge_fast_inte->SetLineWidth(2);

  TMultiGraph *mg_fast = new TMultiGraph();
  mg_fast->SetTitle("Estimated Fast Neutron (1 MeV) Flux;Distance from Ion Source (m);Flux (/(cm^{2} s))");
  mg_fast->Add(ge_fast_norm);
  mg_fast->Add(ge_fast_inte);
  mg_fast->Draw("ALP");
  c1->cd(2)->BuildLegend();


  c1->Update();
  c1->Modified();


}
