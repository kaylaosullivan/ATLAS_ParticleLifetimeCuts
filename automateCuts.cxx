#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TDatabasePDG.h"  
#include "TParticlePDG.h" 



std::ofstream file;

TH1D *hist = new TH1D ();
TChain *t = new TChain();

//int numPlots = 7;
int numPlotted = 0;

//TCanvas *c1 = new TCanvas("c1","",1080,700);


// Areas and ratio
double Asig = 0;
double Abkg = 0;
double Adata = 0;
double BSratio = 0;

// Fit Parameters
double c0 = 0;
double cB = 0;
double rateB = 0;
double cS = 0;
double rateS = 0;
double lifetime = 0;
double chisq = 0;
double error = 0;

bool passed = false;
int status = 0;
bool noCut = true;

// Type of cut
//enum cutType {cosTheta, delR, pT, massK, massL};

struct cutType{
  float min;
  float max;// = 0;
  string name;// = "NoCut";
} cosTheta, delR, pT, massK, massL,none;

struct particle{
  int id;
  double mass;
  string symbol;
  double lifetime;
  double life_err;
} ptcl;



// Branch or variable to draw
//const char *branch = "(RecDelR/1000/(( sqrt(Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt)/sqrt((Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt) + 497.611*497.611))*299792458))>>hist1";

string branch = "";

//const char *branch;






// Fit Functions
double fit(const double *x, const double *par){  // Signal + Background Fit: const + 2 exponentials
    return  TMath::Exp(par[0] + par[1]*x[0]) + TMath::Exp(par[2] + par[3]*x[0]) + par[4];
}

double expon(const double *x, const double *par){ // Signal or Background: single exponential
    return par[0] + TMath::Exp(par[1] + par[2]*x[0]);
}

// Prototypes
void performCut(cutType,bool=false);
void doFit(cutType,int);
void recordResults(cutType,int);
string getTime();

void initializeCuts(){
  cosTheta.min = 0;
  cosTheta.max = 0;
  cosTheta.name = "RecCosTheta";
  
  delR.min = 0;
  delR.max = 0;
  delR.name = "RecDelR";
  
  pT.min = 0;
  pT.max = 0;
  pT.name = "RecPt";
  
  massK.min = 0;
  massK.max = 0;
  massK.name = "RecMass";
  
  massL.min = 0;
  massL.max = 0;
  massL.name = "RecMassLambda";
  
  none.min = 0;
  none.max = 0;
  none.name = "NoCut";
}

void initializePtcl(int id, string symb){


  TDatabasePDG *pdg = new TDatabasePDG();
  TParticlePDG *p = pdg->GetParticle(id);
  
  ptcl.id = id;
  ptcl.mass = p->Mass()*1000; //*1000 to be in MeV
  ptcl.symbol = symb;
  //ptcl.lifetime = p->Lifetime(); // Bug in TParticlePDG, lifetime gives 0. Added manually below.

  if (id == 310){
     ptcl.lifetime = 8.954e-11;
     ptcl.life_err = 0.004e-11;
  }
  if (id == 3122){
    ptcl.lifetime = 2.632e-10;
    ptcl.life_err = 0.02e-10;
  }
  
  file.open (ptcl.symbol+"_Results"+getTime()+".csv");
  initializeCuts();
 // performCut(none);
  
  //delete branch;
    branch = ("RecDelR/1000/(( sqrt(Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt)/sqrt((Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt) +"+  std::to_string(ptcl.mass)+"*"+ std::to_string(ptcl.mass)+"))*299792458)>>hist1").c_str();

  }


//void newParticle

void performCut(cutType cut, bool both = false){    
   //const char* branch = ("RecDelR/1000/(( sqrt(Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt)/sqrt((Pz*Pz*RecPt*RecPt/(1-Pz*Pz)+RecPt*RecPt) +"+  std::to_string(ptcl.mass)+"*"+ std::to_string(ptcl.mass)+"))*299792458)>>hist1").c_str();
    // 1 = only min, 2 = only max, 3 = both at once, 4 = none
    const char* branch2 = branch.c_str();
    //if (cut.name.compare("NoCut")){
    if (noCut){
      t->Draw(branch2);
      doFit(none,4);
      noCut = false;
    }
    if (!both && cut.min != 0){
      t->Draw(branch2,(cut.name+">"+std::to_string(cut.min) ).c_str());
      doFit(cut,1);
    }
    if (!both && cut.max != 0){
      t->Draw(branch2,(cut.name+"<"+std::to_string(cut.max) ).c_str());
      doFit(cut,2);
    }
    if (both){
      t->Draw(branch2,(cut.name+">"+std::to_string(cut.min)+" && " + cut.name+"<"+std::to_string(cut.max)).c_str());
      doFit(cut,3);  
    }
}

void performCuts(cutType cut, float* minCuts, float* maxCuts, int numMin, int numMax, bool both = false){   
 // float numMax = sizeof(maxCuts)/sizeof(maxCuts[0]);
  //float numMin = sizeof(minCuts)/sizeof(minCuts[0]);
  
  
  if(!both){
    cut.max = 0;
    for (int i = 0; i < numMin; i++){
   // cout << "min " << i<< endl;
      cut.min = minCuts[i];
      performCut(cut,false);  
    }
    cut.min = 0;
    for (int i=0; i<numMax; i++){
      cut.max = maxCuts[i];
      performCut(cut,false);
     // cout << "max " << i<< endl;
    }
  } else {
    for (int i =0; i < numMin; i++){
      cut.min = minCuts[i];
      cut.max = maxCuts[i];
      performCut(cut,true);
    }
  }
}
  
void errorCuts(cutType cut, float minCut, float minErr, float maxCut, float maxErr, bool both = false){
  float mins[3] = {(minCut-minErr), minCut, (minCut+minErr)};
  float maxs[3] = {(maxCut-maxErr), maxCut, (maxCut+maxErr)};
  performCuts(cut, mins, maxs,3,3,both);  
}


void doFit(cutType cut, int code){
numPlotted++;

if (code == 2) std::cout << "CUT: " << cut.name << ", MIN: - " << ", MAX: " << cut.max << std::endl;
else if (code == 1)  std::cout << "CUT: " << cut.name << ", MIN: " << cut.min<< ", MAX: " << "-" << std::endl;
else if (code ==3) std::cout << "CUT: " << cut.name << ", MIN: " << cut.min<< ", MAX: " << cut.max << std::endl;
else if (code ==4) std::cout<< "NO CUTS" <<std::endl;

TF1 *fn = new TF1();

int numBins = 189; //200;
       double lowLim = .2e-9; //.1e-9; 
       double highLim = 2e-9;
       
if (ptcl.id==310){

  numBins =178;// 200;//178;//178;
        lowLim = 2e-10;//1e-10;//2e-10;//0.2e-9;
        highLim = 1e-9;
  //double lowFit = 2.148e-10;
  //double highFit = 0.94e-9;
  
  // Fit Function
  fn = new TF1("fn", fit,2.148e-10,.94e-9,5);// 2.148e-10, .94e-9,5);

  // Set Parameter Names
  fn->SetParName(4, "C0");
  fn->SetParName(0, "Cb");
  fn->SetParName(1, "Rate_b");
  fn->SetParName(2, "Cs");
  fn->SetParName(3, "Rate_s");
/*
fn->SetParameter(0,9.67147);
        fn->SetParameter(1,-3.94e9);
        fn->SetParameter(2, 0.00001);//0.002);
        //fn->SetParLimits(2,0,1000);
	fn->SetParameter(3, -1.11e10);
 */
   fn->SetParameter(0,9.67147);
  fn->SetParameter(1,-3.8e9); //3.94
  fn->SetParError(1,1e-5);
  fn->SetParameter(2, 0.002);//001);//0.002);
  //fn->SetParLimits(2,0,1000);
	fn->SetParameter(3, -1.11e10);
   fn->SetParError(3,1e-3);
   //fn->SetParLimits(3,-6e10,-5e9);
 fn->SetParameter(4,-1000);
 fn->SetParLimits(4,-80000,0);

   //fn->FixParameter(3, -1.10630e10);  // Change to SetParameter to re-fit rate value
  fn->SetLineColor(kMagenta);
 }else{
  
       int numBins = 189; //200;
       double lowLim = .2e-9; //.1e-9; 
       double highLim = 2e-9;
  
  fn = new TF1("fn", fit,.175e-9 , 5.8e-9, 5);// .175 - 1.4

        // Set Parameter Names
        fn->SetParName(0, "Cb");
        fn->SetParName(1, "Rate_b");
        fn->SetParName(2, "Cs");
        fn->SetParName(3, "Rate_s");
        fn->SetParName(4, "C0");
        

        // Set Parameter Values
        fn->SetParameter(2, 0.001);
        fn->SetParLimits(2,0,1000);
        fn->SetParameter(3,-3.8e9);  // Rate can vary
         fn->SetParError(3,1e-5);
         fn->SetParameter(1,-6e7);
         fn->SetParError(1,1e-3);
       // fn->FixParameter(3,-3.871e9);  // Rate fixed
        //if (cut.name.compare("RecMass")) fn->SetParameter(0,-200000);
        //if (cut.name.compare("RecCosTheta")) fn->SetParameter(1,-4e7);
        fn->SetLineColor(kMagenta);
  /*      
  //ADD CODE TO PLOT ?
      c1->cd(numPlotted);
      hist->GetXaxis()->SetTitle("Time [s]");
      hist->GetYaxis()->SetTitle("Number of Events per Bin");
      
      hist->SetStats(0); // Remove statistics
      //gStyle->SetOptFit(0101);  // Display statistics
      hist->SetMarkerColor(1);
      hist->SetMarkerStyle(kFullDotLarge);
      hist->SetMinimum(0);     // Set minimum to 0 to see signal without background
      c1->SetRightMargin(1.2); // Extend right margin by 20% to avoid numbers being cut off
  
      // Draw and Fit histogram
      hist->Draw("P");
      hist->Draw("E1SAME");
  */
  }
    TFitResultPtr r = hist->Fit("fn","R");
    status = r;
  //  Int_t fitStatus = r;

    // Plot background function obtained from fit
    TF1 *timeBkgFunc = new TF1("timeBkgFunc", expon, lowLim, highLim, 3);
    timeBkgFunc->SetParameter(0, fn->GetParameter(4));
    timeBkgFunc->SetParameter(1, fn->GetParameter(0));
    timeBkgFunc->SetParameter(2, fn->GetParameter(1));
    timeBkgFunc->SetLineColor(kBlue);
    //timeBkgFunc->Draw("SAME");
            
    // Plot Signal function obtained from fit
    TF1 *timeSignalFunc = new TF1("timeSignalFunc", expon, lowLim, highLim, 3);
    timeSignalFunc->SetParameter(0, 0);
    timeSignalFunc->SetParameter(1, fn->GetParameter(2));
    timeSignalFunc->SetParameter(2, fn->GetParameter(3));
    timeSignalFunc->SetLineColor(kGreen);
    //timeSignalFunc->Draw("SAME");
    //c1->Update();  
    
    Asig = timeSignalFunc->Integral(lowLim,highLim);
    Abkg = timeBkgFunc->Integral(lowLim,highLim);
    BSratio = Abkg/Asig;
    Adata = hist->Integral(1,numBins);
    c0 = fn->GetParameter(4);
    cB = fn->GetParameter(0);
    rateB = fn->GetParameter(1);
    cS = fn->GetParameter(2);
    rateS = fn->GetParameter(3);
    lifetime = -1./rateS;
    chisq = fn->GetChisquare()/fn->GetNDF();
    error = (1/pow(rateS,2.0))*fn->GetParError(3);

    cout<< "Bkg/Sig: " << BSratio << endl;

  recordResults(cut,code);
}

void recordResults(cutType cut, int code){
  if (!passed){
    file << "Cut, Min, Max, Lifetime,Error,Chisq (red),bkg/sig,Total Area,C0,Cb,rate_b,Cs,rate_s,Status \n";
    file << "PDG Value,-,-,"<< (ptcl.lifetime)<<","<< (ptcl.life_err)<<",-,-,-,-,-,-,-,-,- \n";
    passed = true;
  }
  
  if (code == 1) // just min
    file << cut.name<<","<<cut.min<<",-,"<<lifetime<<","<<error<<","<<chisq<<","<<BSratio<<","<<Adata<<","<<c0<<","<<cB<<","<<rateB<<","<<cS<<","<<rateS<<","<<status<<"\n";
  else if (code == 2) // just max
    file << cut.name<<",-,"<<cut.max<<","<<lifetime<<","<<error<<","<<chisq<<","<<BSratio<<","<<Adata<<","<<c0<<","<<cB<<","<<rateB<<","<<cS<<","<<rateS<<","<<status<<"\n";
  else if (code == 3) // both
    file << cut.name<<","<<cut.min<<","<<cut.max<<","<<lifetime<<","<<error<<","<<chisq<<","<<BSratio<<","<<Adata<<","<<c0<<","<<cB<<","<<rateB<<","<<cS<<","<<rateS<<","<<status<<"\n";
   else if (code ==4) // none
      file << cut.name<<",-,-,"<<lifetime<<","<<error<<","<<chisq<<","<<BSratio<<","<<Adata<<","<<c0<<","<<cB<<","<<rateB<<","<<cS<<","<<rateS<<","<<status<<"\n";
}



string getTime(){
    string dateTime;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time ( &rawtime );
    
    timeinfo = localtime ( &rawtime );
    strftime (buffer,80,"%c.txt",timeinfo); 
    dateTime = buffer;
    
    return dateTime.c_str();
}

void automateCuts(){
  // Add chain and branch
  t = new TChain("KsValidation");
  t->Add("/home/steben/projects/rrg-steven/steben/MultiquarkSearch/newtril/run/dataoct26/user.kosulliv.00267358.physics_MinBias.recon.AOD.r10170_tid13005594_00.11_hist.210140925/*.root");
  
  
  //initializePtcl(310,"K"); //automatically initializes cuts as well
  initializePtcl(3122,"L");

  // Create histogram
  int numBins = 1;
  double lowLim = 0;
  double highLim = 1;
  if (ptcl.id==310){
    numBins =178;// 200;//178;//178;
    lowLim = 2e-10;//1e-10;//2e-10;//0.2e-9;
    highLim = 1e-9;
    
  }else if (ptcl.id ==3122){
   numBins = 189; //200;
   lowLim = .2e-9; //.1e-9; 
   highLim = 2e-9;
  }
        
  hist = new TH1D ("hist1","Lifetime Plot", numBins,lowLim,highLim);

 errorCuts(delR,25,1.325,475,13.225);
 errorCuts(pT,100,6.55,10500,254.275);
  errorCuts(massK, 295, 4.15, 490, 4.15);
  errorCuts(massL, 1110, 2.99, 1120, 2.99, true);
  
 
 
  float minR[] = {2,26,27,23};
  float maxR[] = {0};
  float minCos[] = {0.99995,0.99996,0.99997,0.99998,.99999};
  float minPT[] = {50,60,70,80,90};//,350};//,375,425,450};
  float maxPT[] = {0};//,9500};//,9750,10250,10500};//,2500,2750};
  float minKM[] = {290.85,299.15};
  float maxKM[] = {0};
 // float maxKM[] = {494.15,485.85};
  float minLM[] = {1090,1108};//,1109,1111};
  float maxLM[] = {1140,1122};//,1121,1119};
  
 performCuts(massL, minLM, maxLM, sizeof(minLM)/sizeof(minLM[0]), sizeof(maxLM)/sizeof(maxLM[0]),true);
 performCuts(cosTheta,minCos,{0},sizeof(minCos)/sizeof(minCos[0]),0,false);
  performCuts(delR,minR,maxR,sizeof(minR)/sizeof(minR[0]),sizeof(maxR)/sizeof(maxR[0]),false);
 perormCuts(pT,minPT,maxPT,sizeof(minPT)/sizeof(minPT[0]),sizeof(maxPT)/sizeof(maxPT[0]),false);
 performCuts(massK, minKM, maxKM, sizeof(minKM)/sizeof(minKM[0]),sizeof(maxKM)/sizeof(maxKM[0]),false);



  file.close();

  

}
