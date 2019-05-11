// A code to simulate the percentage of neutrons that would penetrate the Au foil
// if it would be emitted in all directions uniformly
// 1. Place a rectangular "foil" at distance d from the neutron source
// 2. Fly neutrons in all directions: uniform 3D spherical distribution
// 3. Count how many of them passed through the foil
// 4. Scan d and plot the d-dependence of "detection efficiency"


#include <iostream>
#include <random>
#include <TMath.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TCanvas.h>

#include <TPolyLine3D.h>
#include <TLatex.h>

#include <TRint.h>

using namespace std;


int main (int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);

	// Number of neutrons to be flown
	int N_n = 10000;
	// Number of samples to be averaged with
	int N_Average = 100;


	random_device rnd;
	default_random_engine engine(rnd());

	normal_distribution<> randnorm(0.,1.);

	TGraph2D *traj = new TGraph2D();
	traj->SetName("traj");
	traj->SetTitle("Simulated Particle Trajectories; x (m); y (m); z (m)");

	TGraphErrors *rdep = new TGraphErrors();	
	rdep->SetTitle("Distance Dependence of Neutron Penetration; Distance r (m); Penetration #epsilon (r) (%)");
	

	double r_start = 0.01; // m
	double r_end = 0.6; // m
	double r_step = 0.01; // m
	int N_steps = round((r_end-r_start)/r_step);

	// Dimensions of the foil
	double F_H = 100.0; // mm
	double F_W = 100.0; // mm

	double x_ll = -F_W/(1000.*2.0); 
	double x_ul = F_W/(1000.*2.0);
	double y_ll = -F_H/(1000.*2.0);
	double y_ul = F_H/(1000.*2.0);

	// Neutron starting point
	double x_0 = 0.0;
	double y_0 = 0.0;
	double z_0 = 0.0;



	double r = r_start;
	// r-scan starts here
	for (int k = 0; k < N_steps; ++k){

		cout << "At r = " << r << " m: ";

		int hits = 0;
		int hitsq = 0;
		double stdev = 0.0;
		// averaging starts here
		for (int j = 0; j < N_Average; ++j){
			int counter = 0;
//			if ((j==0)&&(k==0)){
//				// Draw foil
//				c1->cd(1);
//				TPolyLine3D *foil = new TPolyLine3D(5);
//				foil->SetPoint(0, x_ll, y_ll, r);
//				foil->SetPoint(1, x_ul, y_ll, r);
//				foil->SetPoint(2, x_ul, y_ul, r);
//				foil->SetPoint(3, x_ll, y_ul, r);
//				foil->SetPoint(4, x_ll, y_ll, r);
//				foil->SetLineWidth(4);
//				foil->SetLineColor(2);
//				foil->Draw();
//			}

			for (int i = 0; i < N_n; ++i){
				double a = 0.0;
				double a_x, a_y, a_z;
				while (a == 0.0){
					a_x = randnorm(engine);
					a_y = randnorm(engine);
					a_z = randnorm(engine);
					a = TMath::Sqrt(a_x*a_x + a_y*a_y + a_z*a_z);
				}
	
				TPolyLine3D *trajectory = new TPolyLine3D(-1);
//				if ((j==0)&&(k==0)){
//					trajectory->SetPoint(0,x_0,y_0,z_0);
//					trajectory->SetPoint(1,x_0+r*a_x,y_0+r*a_y,z_0+r*a_z);
//					trajectory->SetLineWidth(1);
//					trajectory->SetLineStyle(2);
//					trajectory->Draw("SAME");
//				}

				if (a_z > 0.0){
					double t_foil = (r-z_0)/a_z;
					double x_f = x_0 + t_foil*a_x;
					double y_f = y_0 + t_foil*a_y;
					bool ul = (x_f <= x_ul)&&(y_f <= y_ul);
					bool ll = (x_f >= x_ll)&&(y_f >= y_ll);
					if ( ul&&ll ){
//						if ((j==0)&&(k==0)){
//							traj->SetPoint(2*counter,x_0,y_0,z_0);
//							traj->SetPoint(2*counter+1,x_f,y_f,r);
//							trajectory->SetLineColor(2);
//						}
						++counter;
					}
				}
			}
	
//			if ((j==0)&&(k==0)){
//				traj->Draw("SAME,P0,ah,fb,bb");
//				traj->GetXaxis()->SetLimits(-0.01,0.01);
//				traj->GetYaxis()->SetLimits(-0.01,0.01);
//				traj->GetZaxis()->SetLimits(r-0.1, r+0.1);
//			}
			hits += counter;
			hitsq += counter*counter;
		}	
		// averaging ends here
		double h_av = double(hits)/double(N_Average);
		double hsq_av = double(hitsq)/double(N_Average);
//		cout << "<hits> = " << h_av << ", <hits^2> = " << hsq_av << endl;
		stdev = TMath::Sqrt(hsq_av - (h_av*h_av));

		double penetration = 100.*double(h_av)/double(N_n);
		double penetration_err = 100.*double(stdev)/double(N_n);
		cout << penetration << " +- " << penetration_err << "% of neutrons penetrated the foil" << endl;
		rdep->SetPoint(k,r,penetration);
		rdep->SetPointError(k,0.,penetration_err);
		r += r_step;
	// r-scan ends here
	}

//	c1->cd(2);
	c1->cd(1);
	rdep->Draw("AL*");
//	c1->cd(2)->SetLogx();
//	c1->cd(2)->SetLogy();
	c1->cd(1)->SetLogx();
	c1->cd(1)->SetLogy();

	TF1* rdep_fit = new TF1("rdep_fit","[0]/(x*x) + [1]");
	rdep_fit->SetParameters(1.0,1.0);
	rdep->Fit("rdep_fit","M","",0.5,0.6);


	TGraphErrors *rdep_residual = new TGraphErrors();
	rdep_residual->SetTitle("Residual of 1/r^{2} Fit; Distance r (m); #epsilon(r) - f(r) (%)");
	for (int i = 0; i < rdep->GetN(); ++i){
		double r = r_start + double(i)*r_step;
		double res = rdep->GetY()[i] - rdep_fit->Eval(r);
		double res_err = rdep->GetEY()[i];
		rdep_residual->SetPoint(i,r,res);
		rdep_residual->SetPointError(i,0.0,res_err);
	}
	c1->cd(2);
	rdep_residual->Draw("AL*");
	
	TF1 *residual_fit = new TF1("residual_fit","1.0 - TMath::Exp([0]/(TMath::Power(x,[1])))");
	residual_fit->SetParameters(1.0,1.0);
	rdep_residual->Fit("residual_fit","M","",0.01,0.1);
//	rdep_residual->Draw("SAME");

	c1->cd(1);
	TF1 *improved_fit = new TF1("improved_fit","[0]/(x*x) + [1] + 1.0 - TMath::Exp([2]/(TMath::Power(x,[3])))");
	improved_fit->SetParameters(rdep_fit->GetParameter(0),rdep_fit->GetParameter(1),residual_fit->GetParameter(0),residual_fit->GetParameter(1));
	improved_fit->SetLineColor(4);
	improved_fit->Draw("SAME");


	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
