#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <sstream>
#include <time.h>
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TMath.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMinuit.h"
#include "TView3D.h"
#include "TView.h" 
#include "TVirtualPad.h" 

// This code does the following :
//  - creates a 2-D function of an annulus.
//  - creates a histogram and fills it randomly via the function.  
//  - draws the function and the histogram.
//  - splits the 2-D histogram in its angular and radial components from the given origin. 
//  - fills a new histogram with the reintegrated components (the precision is not yet the best... )
//  - and leaves the possibility to turn and 
// This example can be executed via the interpreter or ACLIC
//   root > .x rotateAnnulusFunc.C
//   root > .x rotateAnnulusFunc.C++

Double_t g3(Double_t *x, Double_t *par); 
Double_t fun3(Double_t *x, Double_t *par); 
void PixToMMzoom(TH2D * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax); 


bool ShowNicePalette = true; 

void rot_annulus(int i_run) {
    if (ShowNicePalette) {
        int NRGBs = 7, NCont = 99;
        gStyle->SetNumberContours(NCont);
        Double_t stops[7] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
        Double_t red[7]   = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
        Double_t green[7] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
        Double_t blue[7]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
    Double_t w = 700;   Double_t h = 700;
    TCanvas * can_fun = new TCanvas(" Annulus function ", "Annulus function", w, h );                       // name, title, width, height. 
    can_fun  -> Divide(2,2); 

    // Variables and function variables. 
    //Double_t pi = 3.14159265359; 
    //const double pi = 3.1415926535897;
    //Double_t pi = 4*atan(1); 
    const long double pi = 3.141592653589793238L;
    const Int_t npar = 9;
    Double_t theta = 0.478;             // theta = (0.44 ± 0.03) in radiants
    Double_t a = 1.25;                  // a     = (1.31 ± 0.07)mm
    Double_t b = 1.07;                  // b     = (1.15 ± 0.07)mm 
    Double_t xOrigin = 6.14;           // x0    = (6.126±0.016)mm 
    Double_t yOrigin = 7.31;           // z0    = (7.291±0.016)mm

 

    // Bin number hard-coded: 
    Int_t noBin = 1048*3; 
    Double_t z_upper = 1300.;  

   // Read histogram from root file. 
   char nameFile1[30];
   sprintf(nameFile1, "../root_files/fine_binned_gaussian_filtered/%d_fine_bin_gauss_filter.root", i_run);
   //sprintf(nameFile1, "../root_files/raw_data/%d_annulus.root", i_run);
   TFile* file = new TFile(nameFile1, "READ"); 
   TH2D * histo = (TH2D*)file->Get("histo_gauss"); 

    //********************** To be changed.  
    int sizeAxis = histo -> GetXaxis() -> GetXmax();
    cout << "sizeAxis : " << "\t" << sizeAxis << endl; 
    Double_t dimensions = ((1050-2)/43.5)*0.48;
    //Double_t dimensions = ((sizeAxis-2)/43.5)*0.48; 

    Double_t x0 = 0.; Double_t x1 = dimensions;               // good values for a zoom in: Double_t x0 = 2.8; Double_t x1 = 13.2; 
    Double_t y0 = 0.; Double_t y1 = dimensions;               // good values for a zoom in: Double_t y0 = 3.2; Double_t y1 = 13.6; 

    // Convert pixels to mm in the original image and zoom in : 
    PixToMMzoom(histo, x0, x1, y0, y1, 0., z_upper);  //void PixToMMzoom(TH2D * h1, xmin, xmax, ymin, ymax, zmin, zmax); 

   // Get the first, the last bin numbers of the x and y axis, respectively. 
   Int_t binXmin = histo -> GetXaxis() -> GetFirst();   Int_t binXmax = histo -> GetXaxis() -> GetLast();
   Int_t binYmin = histo -> GetYaxis() -> GetFirst();   Int_t binYmax = histo -> GetYaxis() -> GetLast();

    // Calculate the integral of the histoggram.  
   Int_t intHisto = -1; 
   intHisto = histo->Integral(binXmin, binXmax, binYmin, binYmax); 
   cout << "Integral of Histogram : " << intHisto << endl; 

   // Variable declaration for the unraveling of the histogram data in terms of circular coordinates. 
   Double_t tempX = -1.;        Double_t tempY = -1.;    Int_t tempZ = -1.; 
   vector<double> x;            vector<double> y;        vector<int> z; 
   vector<double> al;           vector<double> dx;       vector<double> dy; 
   Double_t r  = -1.;           Double_t alpha = -1.;    Double_t alpha1 = -2;   Double_t alphaDiff = -3; 
   Double_t r1 = -1.;           Double_t r2    = -1.; 
   Double_t tX = -1.;           Double_t tY    = -1.; 
   Double_t ba = sqrt(b/a);     Double_t ab = sqrt(a/b); 

   Double_t widthX = histo -> GetXaxis() -> GetBinWidth(1); 
   Double_t widthY = histo -> GetYaxis() -> GetBinWidth(1); 


   //Loop over all bins: incrementing a y bin, then looping over all corresponding x bins. 
   for(Int_t j = binYmin; j<binYmax+1; j++){
        for(Int_t i = binXmin; i<binXmax+1; i++){
            tempZ = histo -> GetBinContent(i,j); 
            // If a bin corresponding to coordinate (x,y) has content, the bin content is assigned to a new (x',y') coordinate
            if(tempZ!=0){
                tempX  = histo -> GetXaxis() -> GetBinCenter(i);             // x coordinate. 
                //tempX -= widthX; 
                tempY  = histo -> GetYaxis() -> GetBinCenter(j);             // y coordinate. 
                //tempY -= widthY; 
                r1 = tempX - xOrigin; 
                r2 = tempY - yOrigin; 
                r  = sqrt(r1 * r1 + r2 * r2);
                // Retrieving the angle alpha of the coordinate-origin to the x-axis. 
                if((r1 > 0) && (r2 > 0)){                                              // Quadrant I
                    alpha = acos(r1/r); 
                    alpha1 = asin(r2/r); 
                    alphaDiff = alpha-alpha1;
                    //cout << " alpha - alpha1 = " << alphaDiff << endl; 
                    alpha = alpha1; 
 
                    //cout << "alpha/pi [Quadrant I] : " << "\t" << alpha/pi << endl; 
                }
                else if ((r1 < 0) && (r2 > 0)){                                        // Quadrant II
                    alpha = acos(r1/r);
                    alpha1 = pi - asin(r2/r); 
                    alphaDiff = alpha-alpha1;
                    //cout << " alpha - alpha1 = " << alphaDiff << endl; 
                    alpha = alpha1; 
                    //cout << "alpha/pi [Quadrant II] : " << "\t" << alpha/pi << endl;
                }
                else if ((r1 < 0) && (r2 < 0)){                                        // Quadrant III
                    alpha  = pi + acos(abs(r1)/r); 
                    alpha1 = pi + asin(abs(r2)/r); 
                    alphaDiff = alpha-alpha1;
                    //cout << " alpha - alpha1 = " << alphaDiff << endl; 
                    alpha = alpha1; 
                    //cout << "alpha [Quadrant III] : " << "\t" << alpha << endl; 
                }
                else if ((r1 > 0) && (r2 < 0)){                                        // Quadrant IV
                    alpha = 2*pi - acos(abs(r1)/r);
                    alpha1 =2*pi - asin(abs(r2)/r); 
                    alphaDiff = alpha-alpha1;
                    //cout << " alpha - alpha1 = " << alphaDiff << endl; 
                    alpha = alpha1;  
                    //cout << "alpha [Quadrant IV] : " << "\t" << alpha << endl; 
                }
                else { 
                    cout << "Something went wrong. Please check for errors. " << endl; 
                }

                // If theta is substracted from alpha. 
                //alpha += pi/4.; 
                // The transformed coordinates (x',y') are calculated. 
                //tX = xOrigin + r*cos(alpha + theta);   
                //tY = yOrigin + r*sin(alpha + theta); 
                tX = xOrigin + ba*r*cos(alpha-theta);   
                tY = yOrigin + ab*r*sin(alpha-theta);  
                //cout << " x - x' : " << tX - (r1 + xOrigin) << "\t" << " y - y' : " << tY -(r2 + yOrigin) << endl; 
                //cout << " tX : " << tX << "\t" << "r1 + xOrigin : " << r1 + xOrigin << endl; 
                //cout << " tY : " << tY << "\t" << "r1 + yOrigin : " << r2 + yOrigin << endl; 
                // The transformed coordinates (x', y', z'=z ) are written into the vectors x, y, z. 
                x.push_back(tX);
                y.push_back(tY); 
                z.push_back(tempZ); 

                al.push_back(alpha); 
                dx.push_back(tX-r1+xOrigin); 
                dy.push_back(tY-r2+yOrigin);

            }
        }
   }
   

   // Create & Fill the new histogram htrans. 
   TH2D * htrans = new TH2D("htrans", "Transformed Data: ", noBin, x0, x1, noBin, y0, y1);            // name, title, ...

   // Loop over the length of the vectors x, y and z. Their content is filled into the histogram htrans. 
   Int_t N = z.size(); 
   for(Int_t i=0; i<N; i++) {
       htrans -> Fill(x[i], y[i], z[i]); 
   }

   // Set the histogram z-value to an adequate range. 
   htrans -> GetZaxis() -> SetRangeUser(0., z_upper);  

   Int_t intHtrans = -1; 
   intHtrans = htrans->Integral(binXmin, binXmax, binYmin, binYmax); 
   cout << "Integral of Histogram htrans : " << intHtrans << endl; 
   cout << "Integral of Histogram histo  : " << intHisto << endl; 
   if(intHtrans == intHisto){
       cout << " No counts got lost doing the transformation " << endl; 
   }
   else{ cout<< "We lost (or gained?) some counts doing the transformation. "<< endl; }

   // The initial and reconstructed and transformed histograms are being drawn as surf2 plot and surf2 plot from top view. 
   can_fun -> cd(1);    histo -> Draw("surf2"); 
   can_fun -> cd(2);    htrans -> Draw("surf2"); 

   // Draw the histogram from the top view.
   can_fun -> cd(3);    histo -> Draw("surf2");  
   gPad -> Modified(); gPad -> Update(); 
   gPad -> GetView() -> TopView(); 

    // Draw the histogram from the top view. 
   can_fun -> cd(4);    htrans -> Draw("surf2"); 
   gPad -> Modified(); gPad -> Update(); 
   gPad -> GetView() -> TopView(); 

   cout << " bin width Y : "  << "\t" << widthY << "\t" << " bin width X : "  << "\t" << widthX << endl; 

   //int iRun = 79768; 
   //char my_name[60];
   //sprintf(my_name, "../root_files/fine_binned_gaussian_filtered/%d_rotated_annulus_fine_binned.root", iRun);
   TFile* my_file = new TFile( "../root_files/rotated_coord_trans/79768_annulus_fine_binned_rot.root" , "RECREATE");  
   cout << "Did not break yet : Before write histo. " << endl; 
   histo -> Write(); 
   cout << "Did not break yet : Before write htrans. " << endl; 
   htrans -> Write(); 
   cout << "Did not break yet : Before delete my_file. " << endl; 
   delete my_file; 
}

Double_t g3(Double_t *x, Double_t *par) {                                                   // Function of the form of an elliptical annulus. 
    Double_t x1 = Double_t((x[0]-par[7])*cos(par[3])+(x[1]-par[8])*sin(par[3]));            //     It is created from to elliptically parameterised
    Double_t x2 = Double_t(-(x[0]-par[7])*sin(par[3])+(x[1]-par[8])*cos(par[3]));           //     Fermi functions substracted from each other.
    Double_t  a = Double_t(par[1]); 
    Double_t  b = Double_t(par[2]); 
    Double_t r1 = Double_t((x1)*(x1)/(a*a));
    Double_t r2 = Double_t((x2)*(x2)/(b*b));     
    Double_t A1 = Double_t((r1+r2-par[5])/par[4]);
    Double_t A2 = Double_t((r1+r2-par[6])/par[4]);
    return par[0]*(-(1/(exp(A2)+1))+1/(exp(A1)+1));  
  }

  Double_t fun3(Double_t *x, Double_t *par) {                     // Both x and par are vectors. 
    Double_t *p1 = &par[0];                                      // p1 is a pointer to the address of parameter 0. 
    Double_t result = g3(x,p1); 
    return result;
 }

 // Converts the histogram from pixels to mm and sets the range. 
  void PixToMMzoom(TH2D * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    int sizeAxis = h1 -> GetXaxis() -> GetXmax();
    Double_t dimensions = ((sizeAxis-2)/43.5)*0.48;                  
    h1-> GetXaxis() -> SetLimits(0.,dimensions); 	
    h1-> GetYaxis() -> SetLimits(0.,dimensions); 
    h1-> GetXaxis()->SetTitle("mm");		
    h1-> GetYaxis()->SetTitle("mm");
    h1 -> GetZaxis() -> SetRangeUser(zmin, zmax);  
    h1 -> GetXaxis() -> SetRangeUser(xmin, xmax);                                       // Set zAxis to some value you like!!! 
    h1 -> GetYaxis() -> SetRangeUser(ymin, ymax);
    cout << "xmin: " << xmin << "\t" << "xmax: " << xmax << "\t" << (xmax-xmin) << endl; 
    cout << "ymin: " << ymin << "\t" << "ymax: " << ymax << "\t" << (ymax-ymin) << endl; 
    
  }