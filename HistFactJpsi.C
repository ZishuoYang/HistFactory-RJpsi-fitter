#include <iostream>
#include <stdio.h>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TH3.h"
#include "TLatex.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooParamHistFunc.h"
#include "RooHist.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"

#define UNBLIND


TDatime *date = new TDatime();

void HistFactJpsi() {
  using namespace RooFit;
  TLatex *t=new TLatex();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);
  t->SetTextFont(132);
  gROOT->ProcessLine("gStyle->SetLabelFont(132,\"xyz\");");
    gROOT->ProcessLine("gStyle->SetTitleFont(132,\"xyz\");");
      gROOT->ProcessLine("gStyle->SetTitleFont(132,\"t\");");
        gROOT->ProcessLine("gStyle->SetTitleSize(0.08,\"t\");");
          gROOT->ProcessLine("gStyle->SetTitleY(0.970);");
  char substr[128];
  RooRandom::randomGenerator()->SetSeed(date->Get()%100000);
    cout << date->Get()%100000 << endl; //For ToyMC, so I can run multiple copies 
                                        //with different seeds without recompiling

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ; //avoid accidental unblinding!

  // Below: Read histogram file to generate normalization constants required to make
  // each histo normalized to unity. Not totally necessary here, but convenient


  TFile q("Histos_20220901.root");
  TH1 *htemp;
  TString mchistos[8]={"sigmu","sigtau","Psi2SMu","JpsiDx","ChiMu","misID","fakeJpsi","Bd2JpsiX"
                        };
                 double mcN_sigmu,mcN_sigtau,mcN_Psi2SMu,mcN_JpsiDx,mcN_ChicMu,mcN_misID,mcN_fakeJpsi,mcN_Bd2JpsiX;
 double *mcnorms[8]={&mcN_sigmu,&mcN_sigtau,&mcN_Psi2SMu,&mcN_JpsiDx,&mcN_ChicMu,&mcN_misID,&mcN_fakeJpsi,&mcN_Bd2JpsiX};
  for(int i =0; i < 8; i++){
    q.GetObject("h_"+mchistos[i],htemp);
    assert(htemp!=NULL);
    *(mcnorms[i])=1./htemp->Integral();
    cout << "mcN_"+mchistos[i]+" = " << 1./ *(mcnorms[i]) << endl;
  }

  // Useful later to have the bin max and min for drawing
  TH3 *JUNK;
  q.GetObject("h_sigmu",JUNK);
  double Z_low=JUNK->GetZaxis()->GetXmin();
  double Z_high=JUNK->GetZaxis()->GetXmax();
  const int Z_bins = JUNK->GetZaxis()->GetNbins();
  JUNK->SetDirectory(0);
  q.Close();


  TStopwatch sw, sw2, sw3;

  

  TRandom *r3 = new TRandom3(date->Get());

  using namespace RooStats;
  using namespace HistFactory;

  // Many many flags for steering
  /* STEERING OPTIONS */ 
  const bool constrainDstst=true;
  const bool useMinos=true;
  const bool useMuShapeUncerts = true;
  const bool useTauShapeUncerts = false;
  const bool useDststShapeUncerts = false;
  const bool fixshapes = false;
  const bool fixshapesdstst = false;
  const bool dofit = true;
  const bool toyMC = false;
  const bool fitfirst = false;
  const bool slowplots = true;
  const bool BBon3d = true; //flag to enable Barlow-Beeston procedure for all histograms.
                           //Should allow easy comparison of fit errors with and 
			                     //without the technique. 3d or not is legacy from an old
                           //(3+1)d fit configuration
  const bool floatMisIDSmear = true; //turn on/off misID smearing variations
  const int numtoys = 1;
  const int toysize = 384236;
  // Set the prefix that will appear before
  // all output for this measurement
  RooStats::HistFactory::Measurement meas("my_measurement","my measurement");
  meas.SetOutputFilePrefix("results/jpsi");
  meas.SetExportOnly(kTRUE); //Tells histfactory to not run the fit and display
                             //results using its own 

  meas.SetPOI("RawRJpsi");

  // set the lumi for the measurement.
  // only matters for the data-driven
  // pdfs the way I've set it up. in invfb
  // variable rellumi gives the relative luminosity between the
  // data used to generate the pdfs and the sample
  // we are fitting

  // actually, now this is only used for the misID
  meas.SetLumi(1.000);
  meas.SetLumiRelErr(0.005);

  /******* Fit starting constants ***********/

  //ISOLATED FULL RANGE NONN
  //*
  const double expTau =0.02*0.8;//0.252*0.1742*0.781/0.85;
  const double e_iso = 0.314;
  double expMu = 14e3*0.8;
  //*/

  double RelLumi = 0.25;

  RooStats::HistFactory::Channel chan("Jpsimu_kinematic");
  chan.SetStatErrorConfig(1e-5,"Poisson");


  // tell histfactory what data to use
  chan.SetData("h_data_subset", "Histos_20220901.root");


  // Now that data is set up, start creating our samples
  // describing the processes to model the data

  /*********************** Bc->Jpsimunu (NORM) *******************************/

  RooStats::HistFactory::Sample sigmu("h_sigmu","h_sigmu", "Histos_20220901.root");
  if(useMuShapeUncerts)
  {
    sigmu.AddHistoSys("v1mu","h_sigmu_v1m","Histos_20220901.root","","h_sigmu_v1p","Histos_20220901.root","");
    sigmu.AddHistoSys("v2mu","h_sigmu_v2m","Histos_20220901.root","","h_sigmu_v2p","Histos_20220901.root","");
  }
  if(BBon3d) sigmu.ActivateStatError();
  sigmu.SetNormalizeByTheory(kFALSE);
  sigmu.AddNormFactor("Nmu", expMu, 1e-6, 21e3);
  sigmu.AddNormFactor("mcNorm_sigmu", mcN_sigmu, 1e-9, 1.);
  chan.AddSample(sigmu);
  
  /************************* Bc->Jpsitaunu (SIGNAL) *******************************/

  RooStats::HistFactory::Sample sigtau("h_sigtau","h_sigtau", "Histos_20220901.root");
  if(useTauShapeUncerts)
  {
    sigtau.AddHistoSys("v1mu","h_sigtau_v1m","Histos.root","","h_sigtau_v1p","Histos.root","");
  }
  if(BBon3d) sigtau.ActivateStatError();
  sigtau.SetNormalizeByTheory(kFALSE);
  sigtau.AddNormFactor("Nmu",expMu,1e-6,21e3);
  sigtau.AddNormFactor("RawRJpsi",expTau,1e-6,0.2);
  sigtau.AddNormFactor("mcNorm_sigtau", mcN_sigtau, 1e-9, 1.);
  chan.AddSample(sigtau);
  
  /************************* Bc->Psi2SMuNu *******************************/

  RooStats::HistFactory::Sample Psi2SMu("h_Psi2SMu","h_Psi2SMu", "Histos_20220901.root");
  if(BBon3d) Psi2SMu.ActivateStatError();
  Psi2SMu.SetNormalizeByTheory(kFALSE);
  Psi2SMu.AddNormFactor("mcNorm_Psi2SMu", mcN_Psi2SMu, 1e-9, 1.);
  Psi2SMu.AddNormFactor("NPsi2Mu",5e2*0.8,1e-6,10e2);
  chan.AddSample(Psi2SMu);
  
 /************************* Bc->JpsiDx *******************************/

  RooStats::HistFactory::Sample JpsiDx("h_JpsiDx","h_JpsiDx", "Histos_20220901.root");
  if(BBon3d) JpsiDx.ActivateStatError();
  JpsiDx.SetNormalizeByTheory(kFALSE);
  JpsiDx.AddNormFactor("mcNorm_JpsiDx", mcN_JpsiDx, 1e-9, 1.);
  JpsiDx.AddNormFactor("NJpsiDx",10e1*0.8,1e-6,5e2);
  chan.AddSample(JpsiDx);
 
 /************************* Bc->Chi_c{1,2}munu *******************************/

  RooStats::HistFactory::Sample ChicMu("h_ChiMu","h_ChiMu", "Histos_20220901.root");
  if(BBon3d) ChicMu.ActivateStatError();
  ChicMu.SetNormalizeByTheory(kFALSE);
  ChicMu.AddNormFactor("mcNorm_ChicMu", mcN_ChicMu, 1e-9, 1.);
  ChicMu.AddNormFactor("NChicMu",7e2*0.8,1e-6,10e2);
  chan.AddSample(ChicMu);
 
 /************************* Bu,d->JpsiX *******************************/

  RooStats::HistFactory::Sample Bd2JpsiX("h_Bd2JpsiX","h_Bd2JpsiX", "Histos_20220901.root");
  if(BBon3d) Bd2JpsiX.ActivateStatError();
  Bd2JpsiX.SetNormalizeByTheory(kFALSE);
  Bd2JpsiX.AddNormFactor("mcNorm_Bd2JpsiX", mcN_Bd2JpsiX, 1e-9, 1.);
  Bd2JpsiX.AddNormFactor("NBd2JpsiX",35e2*0.8,1e-6,56e2);
  chan.AddSample(Bd2JpsiX);

  /*********************** MisID BKG (FROM DATA)  *******************************/

  RooStats::HistFactory::Sample misID("h_misID","h_misID_smear", "Histos_20220901.root");
  if(BBon3d) misID.ActivateStatError();
  misID.SetNormalizeByTheory(kFALSE);
  //misID.AddNormFactor("NmisID",RelLumi,1e-6,1e4);
  if (floatMisIDSmear){ // MisID smearing
	misID.AddHistoSys("misID_smear","h_misID","Histos_20220901.root","","h_misID_doublesmear","Histos_20220901.root","");
  }
  misID.AddNormFactor("mcNorm_misID", mcN_misID, 1e-9, 1.);
  misID.AddNormFactor("NmisID",7e3*0.7,1e-6,21e3);
  chan.AddSample(misID);

  /*********************** Fake Jpsi BKG (FROM DATA)  *******************************/

  RooStats::HistFactory::Sample fakeJpsi("h_fakeJpsi","h_fakeJpsi", "Histos_20220901.root");
  if(BBon3d) fakeJpsi.ActivateStatError();
  fakeJpsi.SetNormalizeByTheory(kFALSE);
  fakeJpsi.AddNormFactor("mcNorm_fakeJpsi", mcN_fakeJpsi, 1e-9,1.);
  fakeJpsi.AddNormFactor("NfakeJpsi",1166.,1166.,1166.);
  chan.AddSample(fakeJpsi);
 
  /************************* B0->D1munu **************************************/
  /*
  RooStats::HistFactory::Sample d1mu("h_D1","h_D1", "DemoHistos.root");
  if(BBon3d) d1mu.ActivateStatError();
  if(useDststShapeUncerts)
  {
    d1mu.AddHistoSys("IW","h_D1IWp","DemoHistos.root","","h_D1IWm","DemoHistos.root","");
  }

  d1mu.SetNormalizeByTheory(kFALSE);
  d1mu.AddNormFactor("mcNorm_D1", mcN_D1, 1e-9, 1.);
  if (!constrainDstst)
  {
    d1mu.AddNormFactor("ND1",1e2,1e-6,1e5);
  }
  else
  {
    d1mu.AddNormFactor("NDstst0",0.102,1e-6,1e0);
    d1mu.AddNormFactor("Nmu",expMu,1e-6,1e6);
    d1mu.AddNormFactor("fD1",3.2,3.2,3.2);
    d1mu.AddOverallSys("BFD1",0.9,1.1);
  }
  chan.AddSample(d1mu);
  */

  /****** END SAMPLE CHANNELS *******/

  sw3.Reset();
  sw3.Start();
  meas.AddChannel(chan);
  ///RooStats::HistFactory::Channel chan2(chan);
  ///chan2.SetName("FakeD");
  //meas.AddChannel(chan2);

  meas.CollectHistograms();
  /*meas.AddConstantParam("mcNorm_sigmu");
  meas.AddConstantParam("mcNorm_sigtau");
  meas.AddConstantParam("mcNorm_D1");
  meas.AddConstantParam("fD1");
  meas.AddConstantParam("NDstst0");
  meas.AddConstantParam("NmisID");
  */


  RooWorkspace *w;
  w=RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
  
  ModelConfig *mc = (ModelConfig*) w->obj("ModelConfig"); // Get model manually
  RooSimultaneous *model = (RooSimultaneous*)mc->GetPdf();

  ///PiecewiseInterpolation *theIW = (PiecewiseInterpolation*) w->obj("h_D1_Jpsimu_kinematic_Hist_alpha");
  //theIW->disableCache(kTRUE);
  ///theIW->Print("V");
  

  RooRealVar* poi = (RooRealVar*) mc->GetParametersOfInterest()->createIterator()->Next();
  std::cout << "Param of Interest: " << poi->GetName() << std::endl;  


  // Lets tell roofit the right names for our histogram variables //
  RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  RooRealVar *x = (RooRealVar*) obs->find("obs_x_Jpsimu_kinematic");
  RooRealVar *y = (RooRealVar*) obs->find("obs_y_Jpsimu_kinematic");
  RooRealVar *z = (RooRealVar*) obs->find("obs_z_Jpsimu_kinematic");
  x->SetTitle("m^{2}_{miss}");
  x->setUnit("MeV^{2}");
  y->SetTitle("lifetime");
  y->setUnit("ps");
  z->SetTitle("Z(q^{2},E_{#mu})");
  z->setUnit("");

  // For simultaneous fits, this is the category histfactory uses to sort the channels

  RooCategory *idx = (RooCategory*) obs->find("channelCat");
  RooAbsData *data = (RooAbsData*) w->data("obsData");
  

/* FIX SOME MODEL PARAMS */
  for(int i =0; i < 8; i++){
    if (((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_"+mchistos[i])))!=NULL)
    {
      ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_"+mchistos[i])))->setConstant(kTRUE);
      cout << "mcNorm_"+mchistos[i]+" = " << ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_"+mchistos[i])))->getVal() << endl;
    }
  }

///  ((RooRealVar*)(mc->GetNuisanceParameters()->find("NDstst0")))->setVal(0.102);
///  ((RooRealVar*)(mc->GetNuisanceParameters()->find("NDstst0")))->setConstant(kTRUE);
///  ((RooRealVar*)(mc->GetNuisanceParameters()->find("fD1")))->setConstant(kTRUE);
///  ((RooRealVar*)(mc->GetNuisanceParameters()->find("NmisID")))->setConstant(kTRUE);
     ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_misID")))->setConstant(kTRUE);
     ((RooRealVar*)(mc->GetNuisanceParameters()->find("mcNorm_ChicMu")))->setConstant(kTRUE);
     ((RooRealVar*)(mc->GetNuisanceParameters()->find("NfakeJpsi")))->setConstant(kTRUE);
     ((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);

///  if(useDststShapeUncerts) ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_IW")))->setRange(-3.0,3.0);  
///  if(useMuShapeUncerts) ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v1mu")))->setRange(-8,8);
///  if(useMuShapeUncerts) ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v2mu")))->setRange(-8,8);
///  if(useMuShapeUncerts) ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v3mu")))->setRange(-8,8);
///  ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_BFD1")))->setRange(-3,3);


  if(fixshapes)
  {
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v1mu")))->setVal(1.06);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v1mu")))->setConstant(kTRUE);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v2mu")))->setVal(-0.159);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v2mu")))->setConstant(kTRUE);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v3mu")))->setVal(-1.75);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v3mu")))->setConstant(kTRUE);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v4tau")))->setVal(0.0002);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v4tau")))->setConstant(kTRUE);
  }
  if(fixshapesdstst)
  {
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_IW")))->setVal(-0.005);//-2.187);
    ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_IW")))->setConstant(kTRUE);
  }
  
  // This switches the model to a class written to handle analytic Barlow-Beeston lite.
  // Otherwise, every bin gets a minuit variable to minimize over!
  // This class, on the other hand, allows a likelihood where the bin parameters
  // are analyitically minimized at each step
  HistFactorySimultaneous* model_hf = new HistFactorySimultaneous( *model );
 
  RooFitResult *toyresult;
  RooAbsReal *nll_hf; 

  RooFitResult *result, *result2;

  cerr << "Saving PDF snapshot" << endl;
  RooArgSet *allpars;
  allpars=(RooArgSet*)((RooArgSet*) mc->GetNuisanceParameters())->Clone();
  allpars->add(*poi);
  RooArgSet *constraints;
  constraints = (RooArgSet*) mc->GetConstraintParameters();
  if(constraints != NULL) allpars->add(*constraints);
  w->saveSnapshot("TMCPARS",*allpars,kTRUE);
  RooRealVar poierror("poierror","poierror",0.00001,0.010);
  TIterator *iter = allpars->createIterator();
  RooAbsArg *tempvar;
  RooArgSet *theVars = (RooArgSet*) allpars->Clone();
  theVars->add(poierror);
  RooDataSet *toyresults = new RooDataSet("toyresults","toyresults",*theVars,StoreError(*theVars));
  RooDataSet *toyminos = new RooDataSet("toyminos","toyminos",*theVars,StoreError(*theVars));
  
  // The following code is very messy. Sorry.
  if (toyMC)
  {
    double checkvar;
    double checkvarmean=0;

    if(fitfirst)
    {
      nll_hf = model_hf->createNLL(*data);
      RooMinuit* minuit_temp = new RooMinuit(*nll_hf) ;
      minuit_temp->setPrintLevel(-1);
      minuit_temp->optimizeConst(1);
      minuit_temp->setErrorLevel(0.5); // 1-sigma = DLL of 0.5
      minuit_temp->setOffsetting(kTRUE);
      minuit_temp->fit("smh");
      poi->setVal(expTau*0.85);
      checkvar = ((RooRealVar*)allpars->find("Nmu"))->getVal();
      cout << checkvar << endl;
    }
    else
    {
      //lets set some params
      if(useMuShapeUncerts)
      {
        ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v1mu")))->setVal(1.06);
        ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v2mu")))->setVal(-0.159);
        ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v3mu")))->setVal(-1.75);
      }
      if(useTauShapeUncerts)
      {
        ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_v4tau")))->setVal(0.0002);
      }
      if(useDststShapeUncerts)
      {
        ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_IW")))->setVal(0.2);//-2.187);
      }
    }
  

      w->saveSnapshot("GENPARS",*allpars,kTRUE);
      
    RooDataSet* datanom;
    cerr << "Attempting to generate toyMC..." << endl;
    sw2.Reset();
    sw2.Start();
    double running_mean=0;
    double running_RMS=0;
    for (int runnum =0; runnum<numtoys ; runnum++)
    {
      w->loadSnapshot("GENPARS");
      cerr << "DEBUG CHECK: Ntau=" << poi->getVal() << endl; 
      cout << "PROGRESS: SAMPLE NUMBER " << runnum << " STARTING GENERATION... " << endl;
      RooDataSet* data2 = model->generate(RooArgSet(*x,*y,*z,model->indexCat()),Name("test"),AllBinned(),NumEvents(toysize),Extended());
      double wsum=0.;
      data2->Print();
      cout << "DONE" << endl;
      w->loadSnapshot("TMCPARS");
      nll_hf=model_hf->createNLL(*data2);
      RooMinuit *minuit_toy = new RooMinuit(*nll_hf);
      minuit_toy->optimizeConst(1);
      //minuit_toy->setPrintLevel(-1);
      minuit_toy->setOffsetting(kFALSE);
      minuit_toy->setOffsetting(kTRUE);
      minuit_toy->setErrorLevel(0.5);
      minuit_toy->setStrategy(2);
      //minuit_toy->setEps(5e-16);
      cout << "PROGRESS: FITTING SAMPLE " << runnum << " NOW...\t";
      minuit_toy->fit("smh");
       cout << " DONE. SAVING RESULT.\n" << endl;
      checkvarmean += ((RooRealVar*)allpars->find("Nmu"))->getVal();
      cout << "DEBUG CHECK: <Nmu>-Nmu_input = " << checkvarmean/(runnum+1) - checkvar << endl;
      poierror.setVal(poi->getError());
      toyresults->add(*theVars);
      cout << "PROGRESS: ATTEMPTING MINOS FOR SAMPLE " << runnum << " NOW...\t";
      minuit_toy->minos(RooArgSet(*poi));
      RooFitResult *testresult = minuit_toy->save("TOY","TOY");
      result=testresult;
      double edm=testresult->edm();
      if(edm > 0.1) 
      {
        cout << "BAD FIT. SKIPPING..." << endl;
        continue;
      }
      cout << " DONE. SAVING RESULT.\n" << endl;
      toyminos->add(*theVars);
      delete minuit_toy;
      delete nll_hf;
      if(runnum+1<numtoys) 
      {
        delete data2;
      delete testresult;

      }
      data=data2;
      double pulli=(poi->getVal()-expTau*0.85)/poi->getError();
      running_mean+=pulli;
      running_RMS+=pulli*pulli;
      cout << "RUNNING MEAN IS " << running_mean/(runnum+1) << "\tRUNNING RMS IS " << sqrt(running_RMS/(runnum+1)) << endl;
    }
    sw2.Stop();
    RooFormulaVar poi_pull("poi_pull","(@0-0.03428)/@1",RooArgList(*poi,poierror));
    RooRealVar *pulls = (RooRealVar*)toyminos->addColumn(poi_pull);
    pulls->setRange(-5,5);
    pulls->setBins(50);
    RooRealVar fitmu("fitmu","#mu",0.,-5,5);
    RooRealVar fitsig("fitsig","#sigma",1.,0.1,10);
    RooGaussian gaus("gaus","gaus",*pulls,fitmu,fitsig);
    gaus.fitTo(*toyminos,Range(-5,5));
    TCanvas *toytest=new TCanvas("toytest","toytest");
    toytest->Divide(2,1);
    toytest->cd(1);
    RooPlot *testframe = pulls->frame(Title("POI Pull"));
    toyminos->plotOn(testframe);
    gaus.plotOn(testframe);
    gaus.paramOn(testframe);
    testframe->Draw();
    //toytest->cd(2);
    //RooPlot *errframe = poierror.frame(Title("POI Error"));
    //toyresults->plotOn(errframe, Cut("poi_pull > -5 && poi_pull < 5"));
    //errframe->Draw();
    toytest->cd(2);
    RooPlot *valframe = poi->frame(Title("POI"));
    toyresults->plotOn(valframe);
    valframe->Draw();
  }


    if(dofit)
    {//return;
      nll_hf= model_hf->createNLL(*data,Offset(kTRUE), Verbose(kTRUE));

      RooMinuit* minuit_hf = new RooMinuit(*nll_hf) ;
      //minuit_hf->setVerbose(kTRUE);
      RooArgSet *temp = new RooArgSet();
      nll_hf->getParameters(temp)->Print("V");
      cout << "**********************************************************************" << endl;
      minuit_hf->setErrorLevel(0.5);
      #ifndef UNBLIND
        minuit_hf->setPrintLevel(-1);
      #endif


      std::cout << "Minimizing the Minuit (Migrad)" << std::endl;
      if(toyMC) {w->loadSnapshot("TMCPARS");}
      else {w->saveSnapshot("TMCPARS",*allpars,kTRUE);}
      sw3.Stop();
      sw.Reset();
      sw.Start();
      minuit_hf->setStrategy(2);
      minuit_hf->fit("smh");
      RooFitResult *tempResult=minuit_hf->save("TempResult","TempResult");

      cout << tempResult->edm() << endl;
      if (useMinos) minuit_hf->minos(RooArgSet(*poi));
      sw.Stop();
      result = minuit_hf->save("Result","Result");

      //EXAMPLE LL SCAN
      /*
      RooAbsReal *pll=nll_hf->createProfile(*poi);
      RooPlot *testFrame = poi->frame(Bins(10),Title("LL"),Range(0.02,0.06));
      // alternately RooPlot *testFrame = poi->frame(Bins(20),Title("LL"), Range(range_low,range_hi)); for restricted range
      pll->plotOn(testFrame,ShiftToZero());
      RooFormulaVar derp("derp","0.5+(@0-(0.03778+0.004301))*(@0-(0.03778-0.004219))*0.5/(0.004219*0.004301)",RooArgList(*poi));
      derp.plotOn(testFrame,ShiftToZero(),LineColor(kRed),LineStyle(kDashed));
      TCanvas *cLL=new TCanvas("cLL","cLL");
      testFrame->Draw();

      gROOT->ProcessLine(".q");
*/
        
    }

    RooPlot *mm2_frame = x->frame(Title("m^{2}_{miss}"));
    RooPlot *Dt_frame = y->frame(Title("decaytime"));
    RooPlot *Z_frame = z->frame(Title("Z"));
    RooPlot *mm2Z_frame[Z_bins];
    RooPlot *DtZ_frame[Z_bins];

    const int nframes = 3;
    RooPlot *drawframes[nframes] = {mm2_frame, Dt_frame, Z_frame};
    RooPlot *Zframes[2*Z_bins];
    RooPlot *Zbframes[2*Z_bins];

    for (int i=0; i < Z_bins; i++)
    {
      mm2Z_frame[i] = x->frame();
      DtZ_frame[i] = y->frame();
      Zframes[i] = mm2Z_frame[i];
      Zframes[i+Z_bins] = DtZ_frame[i];
      Zbframes[i] = x->frame();
      Zbframes[i+Z_bins] = y->frame();
    }

    const int ncomps = 10;


if(result != NULL)
    {
      printf("Fit ran with status %d\n",result->status());

      printf("Stat error on R(Jpsi) is %f\n",poi->getError());

      printf("EDM at end was %f\n",result->edm());
      
      result->floatParsInit().Print();
      
      cout << "CURRENT NUISANCE PARAMETERS:" << endl;
      //TIterator *paramiter = mc->GetNuisanceParameters()->createIterator();
      TIterator *paramiter = result->floatParsFinal().createIterator();
      RooRealVar *__temp= (RooRealVar *)paramiter->Next();
      int final_par_counter=0;
      while (__temp!=NULL)
      {
        if(!__temp->isConstant())
        {
          if(!(TString(__temp->GetName()).EqualTo(poi->GetName())))
          {
            cout << final_par_counter << ": "
            << __temp->GetName() << "\t\t\t = "
            << ((RooRealVar*)result->floatParsFinal().find(__temp->GetName()))->getVal()
            << " +/- "
            << ((RooRealVar*)result->floatParsFinal().find(__temp->GetName()))->getError() << endl;
          }
        }
        final_par_counter++;
        __temp=(RooRealVar *)paramiter->Next();
      }
      
      
      result->correlationMatrix().Print();
    
    if (dofit) printf("Stopwatch: fit ran in %f seconds with %f seconds in prep\n",sw.RealTime(), sw3.RealTime());
    //theIW->_cacheMgr.Print("V");
    //w->Print();
    //return;
    //gROOT->ProcessLine(".q");
  }
    if (toyMC) 
    {
      printf("Stopwatch: Generated test data in %f seconds\n",sw2.RealTime());

    }
      int colors[ncomps]={kRed,kBlue+1,kViolet,kViolet+1,kViolet+2,kGreen,kGreen+2,kOrange+1,kOrange+2,kOrange+3};
      const int ncomps2 = 8;
      TString names[ncomps2+1] = {"Data","Total Fit"
                                ,"B #rightarrow D*#mu#nu"
                                ,"B #rightarrow D**#mu#nu"
                                ,"B #rightarrow D**#tau#nu"
                                ,"B #rightarrow D*[D_{q} #rightarrow #mu#nuX]Y"
                                ,"Combinatoric (wrong-sign)"
                                ,"Misidentification BKG"
                                ,"Wrong-sign slow #pi"
      };
      
      
      RooHist* mm2resid;// = mm2_frame->pullHist() ;
      RooHist* Dtresid;// = Dt_frame->pullHist() ;
      RooHist* Zresid;// = Z_frame->pullHist() ;

      RooHist *resids[nframes];

      std::cout << "Starting plots" <<std::endl;
      for (int i = 0; i < nframes; i++){
        data->plotOn(drawframes[i],DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.4),DrawOption("ZP"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kRed));
        resids[i]=drawframes[i]->pullHist();
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kOrange-8),Components("*misID*,*sigmu*,*Psi2SMu*,*JpsiDx*,*fakeJpsi*,*Chi*,*Bd2JpsiX*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kGreen+2),Components("*misID*,*sigmu*,*Psi2SMu*,*JpsiDx*,*fakeJpsi*,*Chi*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kOrange+1),Components("*misID*,*sigmu*,*Psi2SMu*,*JpsiDx*,*fakeJpsi*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kGreen),Components("*misID*,*sigmu*,*Psi2SMu*,*JpsiDx*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kViolet),Components("*misID*,*sigmu*,*Psi2SMu*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kBlue+1),Components("*misID*,*sigmu*"));
        model_hf->plotOn(drawframes[i], Slice(*idx),ProjWData(*idx,*data),DrawOption("F"),FillColor(kOrange),Components("*misID*"));
        data->plotOn(drawframes[i],DataError(RooAbsData::Poisson),Cut("channelCat==0"),MarkerSize(0.4),DrawOption("ZP"));

      }

      mm2resid=resids[0];
      Dtresid=resids[1];
      Zresid=resids[2];
     
      char cutstrings[Z_bins][128];
      char rangenames[Z_bins][32];
      char rangelabels[Z_bins][128];
      RooHist *mm2Z_pulls[Z_bins];
      RooHist *DtZ_pulls[Z_bins];

      for (int i=0; i < Z_bins; i++)
      {
        double binlow = Z_low+i*(Z_high-Z_low)/Z_bins;
        double binhigh = Z_low+(i+1)*(Z_high-Z_low)/Z_bins;
        //sprintf(rangelabels[i],"%.2f < q^{2} < %.2f",binlow*1e-6,binhigh*1e-6);
        sprintf(rangelabels[i],"Z = %.0f",binlow);
        sprintf(cutstrings[i],"obs_z_Jpsimu_kinematic > %f && obs_z_Jpsimu_kinematic < %f && channelCat==0", Z_low+i*(Z_high-Z_low)/Z_bins, Z_low+(i+1)*(Z_high-Z_low)/Z_bins);
        sprintf(rangenames[i],"Zbin_%d",i);
        z->setRange(rangenames[i],binlow,binhigh);
      }

      if(slowplots == true)
      {
        cout << "Drawing Slow Plots" << endl;
        for (int i = 0; i < Z_bins; i++)
        {
          data->plotOn(mm2Z_frame[i],Cut(cutstrings[i]),DataError(RooAbsData::Poisson),MarkerSize(0.4),DrawOption("ZP"));
          data->plotOn(DtZ_frame[i],Cut(cutstrings[i]),DataError(RooAbsData::Poisson),MarkerSize(0.4),DrawOption("ZP"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kRed));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kRed));
          
          //Grab pulls
          mm2Z_pulls[i]=mm2Z_frame[i]->pullHist();
          DtZ_pulls[i]=DtZ_frame[i]->pullHist();
          
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kOrange-8),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*,*ChicMu*,*Bd2JpsiX*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kOrange-8),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*,*ChicMu*,*Bd2JpsiX*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kGreen+2),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*,*ChicMu*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kGreen+2),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*,*ChicMu*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kOrange+1),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kOrange+1),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*,*fakeJpsi*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kGreen),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kGreen),Components("*sigmu*,*Psi2SMu*,*JpsiDx*,*misID*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kViolet),Components("*sigmu*,*Psi2SMu*,*misID*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kViolet),Components("*sigmu*,*Psi2SMu*,*misID*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kBlue+1),Components("*sigmu*,*misID*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kBlue+1),Components("*sigmu*,*misID*"));
          model_hf->plotOn(mm2Z_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]),DrawOption("F"),FillColor(kOrange),Components("*misID*"));
          model_hf->plotOn(DtZ_frame[i], Slice(*idx),ProjWData(*idx,*data),ProjectionRange(rangenames[i]), DrawOption("F"),FillColor(kOrange),Components("*misID*"));
          data->plotOn(mm2Z_frame[i],Cut(cutstrings[i]),DataError(RooAbsData::Poisson),MarkerSize(0.4),DrawOption("ZP"));
          data->plotOn(DtZ_frame[i],Cut(cutstrings[i]),DataError(RooAbsData::Poisson),MarkerSize(0.4),DrawOption("ZP"));
        }
      }
  TCanvas *c1 = new TCanvas("c1","c1",1000,300);
  c1->SetTickx();
  c1->SetTicky();
  c1->Divide(3,1);
  TVirtualPad *curpad;
  curpad=c1->cd(1);
  curpad->SetTickx();
  curpad->SetTicky();
  curpad->SetRightMargin(0.02);
  curpad->SetLeftMargin(0.20);
  curpad->SetTopMargin(0.02);
  curpad->SetBottomMargin(0.13);
  mm2_frame->SetTitle("");
  mm2_frame->GetXaxis()->SetLabelSize(0.06);
  mm2_frame->GetXaxis()->SetTitleSize(0.06);
  mm2_frame->GetYaxis()->SetLabelSize(0.06);
  mm2_frame->GetYaxis()->SetTitleSize(0.06);
  mm2_frame->GetYaxis()->SetTitleOffset(1.75);
  mm2_frame->GetXaxis()->SetTitleOffset(0.9);
  TString thetitle=mm2_frame->GetYaxis()->GetTitle();
  thetitle.Replace(0,6,"Candidates");
  mm2_frame->GetYaxis()->SetTitle(thetitle);
  mm2_frame->Draw();
  t->DrawLatex(8.7,mm2_frame->GetMaximum()*0.95,"");
  curpad=c1->cd(2);
  curpad->SetTickx();
  curpad->SetTicky();
  curpad->SetRightMargin(0.02);
  curpad->SetLeftMargin(0.20);
  curpad->SetTopMargin(0.02);
  curpad->SetBottomMargin(0.13);
  Dt_frame->SetTitle("");
  Dt_frame->GetXaxis()->SetLabelSize(0.06);
  Dt_frame->GetXaxis()->SetTitleSize(0.06);
  Dt_frame->GetYaxis()->SetLabelSize(0.06);
  Dt_frame->GetYaxis()->SetTitleSize(0.06);
  Dt_frame->GetYaxis()->SetTitleOffset(1.75);
  Dt_frame->GetXaxis()->SetTitleOffset(0.9);
  thetitle=Dt_frame->GetYaxis()->GetTitle();
  thetitle.Replace(0,6,"Candidates");
  Dt_frame->GetYaxis()->SetTitle(thetitle);
  Dt_frame->Draw();
  t->DrawLatex(2250,Dt_frame->GetMaximum()*0.95,"");
  curpad=c1->cd(3);
  curpad->SetTickx();
  curpad->SetTicky();
  curpad->SetRightMargin(0.02);
  curpad->SetLeftMargin(0.20);
  curpad->SetTopMargin(0.02);
  curpad->SetBottomMargin(0.13);
  Z_frame->SetTitle("");
  Z_frame->GetXaxis()->SetLabelSize(0.06);
  Z_frame->GetXaxis()->SetTitleSize(0.06);
  Z_frame->GetYaxis()->SetLabelSize(0.06);
  Z_frame->GetYaxis()->SetTitleSize(0.06);
  Z_frame->GetYaxis()->SetTitleOffset(1.75);
  Z_frame->GetXaxis()->SetTitleOffset(0.9);
  thetitle=Z_frame->GetYaxis()->GetTitle();
  thetitle.Replace(0,6,"Candidates");
  Z_frame->GetYaxis()->SetTitle(thetitle);
  Z_frame->Draw();
  t->DrawLatex(11.1e6,Z_frame->GetMaximum()*0.95,"");


      RooPlot *mm2_resid_frame=x->frame(Title(" "));
      RooPlot *Dt_resid_frame=y->frame(Title(" "));
      RooPlot *Z_resid_frame=z->frame(Title(" "));
      RooPlot *DOCA_resid_frame;

      cerr << __LINE__ << endl;
      mm2_resid_frame->addPlotable(mm2resid,"P");
      cerr << __LINE__ << endl;
      Dt_resid_frame->addPlotable(Dtresid,"P");
      cerr << __LINE__ << endl;
      Z_resid_frame->addPlotable(Zresid,"P");
      cerr << __LINE__ << endl;

      TCanvas *c3 = new TCanvas("c3","c3",1000,200);
      c3->Divide(3,1);
      c3->cd(1);
      c3->cd(1)->SetGridy(1);
      mm2_resid_frame->SetMinimum(-3);
      mm2_resid_frame->SetMaximum(3);
      mm2_resid_frame->SetLabelSize(0.09,"Y");
      mm2_resid_frame->Draw();
      c3->cd(2);
      c3->cd(2)->SetGridy(1);
      Dt_resid_frame->SetMinimum(-3);
      Dt_resid_frame->SetMaximum(3);
      Dt_resid_frame->SetLabelSize(0.09,"Y");
      Dt_resid_frame->Draw();
      c3->cd(3);
      c3->cd(3)->SetGridy(1);
      Z_resid_frame->SetMinimum(-3);
      Z_resid_frame->SetMaximum(3);
      Z_resid_frame->SetLabelSize(0.09,"Y");
      Z_resid_frame->Draw();
      
      TCanvas *c2;
      if(slowplots == true)
      {

        c2 = new TCanvas("c2","c2",1200,600);
        c2->Divide(Z_bins,2);
        double max_scale=1.05;
        double max_scale2=1.05;
        char thename[32];
        for (int k = 0 ; k < Z_bins*2; k++)
        {
          c2->cd(k+1);
          /*
          Zframes[k]->SetTitle(rangelabels[(k % Z_bins)]);
          Zframes[k]->Draw();*/
          sprintf(thename,"bottompad_%d",k);
          //c2->cd((k<Z_bins)*(2*k+1)+(k>=Z_bins)*(2*(k+1-Z_bins)));
          TPad *padbottom = new TPad(thename,thename,0.,0.,1.,0.3);
          
          padbottom->SetFillColor(0);
          padbottom->SetGridy();
          padbottom->SetTickx();
          padbottom->SetTicky();
          padbottom->SetFillStyle(0);
          padbottom->Draw();
          padbottom->cd();
          padbottom->SetLeftMargin(padbottom->GetLeftMargin()+0.08);
          padbottom->SetTopMargin(0);//0.01);
          padbottom->SetRightMargin(0.04);
          //padbottom->SetBottomMargin(padbottom->GetBottomMargin()+0.23);
          padbottom->SetBottomMargin(0.5);
          
          //c2b->cd(k+1);
          TH1 *temphist2lo, *temphist2, *tempdathist;
          RooHist * temphist;
          if(k<Z_bins) {temphist=mm2Z_pulls[k];}
          else {temphist=DtZ_pulls[k-Z_bins];}
          temphist->SetFillColor(kBlue);
          temphist->SetLineColor(kWhite);
          Zbframes[k]->SetTitle(Zframes[k]->GetTitle());
          Zbframes[k]->addPlotable(temphist,"B");
          Zbframes[k]->GetXaxis()->SetLabelSize(0.33*0.22/0.3);
          Zbframes[k]->GetXaxis()->SetTitleSize(0.36*0.22/0.3);
          //Zbframes[k]->GetXaxis()->SetTitle("");
          Zbframes[k]->GetXaxis()->SetTickLength(0.10);
          Zbframes[k]->GetYaxis()->SetTickLength(0.05);
          Zbframes[k]->SetTitle("");
          Zbframes[k]->GetYaxis()->SetTitleSize(0.33*0.22/0.3);
          Zbframes[k]->GetYaxis()->SetTitle("Pulls");
          Zbframes[k]->GetYaxis()->SetTitleOffset(0.2);
          Zbframes[k]->GetXaxis()->SetTitleOffset(0.78);
          Zbframes[k]->GetYaxis()->SetLabelSize(0.33*0.22/0.3);
          Zbframes[k]->GetYaxis()->SetLabelOffset(99);
          Zbframes[k]->GetYaxis()->SetNdivisions(205);
          Zbframes[k]->Draw();
          Zbframes[k]->Draw();
          double xloc=-2.25;
          if(k>=Z_bins) xloc=50;
          t->SetTextSize(0.33*0.22/0.3);
          t->DrawLatex(xloc,-2,"-2");
          //t->DrawLatex(xloc,0," 0");
          t->DrawLatex(xloc*0.99,2," 2");
       
          
          
          c2->cd(k+1);
          //c2->cd((k<Z_bins)*(2*k+1)+(k>=Z_bins)*(2*(k+1-Z_bins)));
          sprintf(thename,"toppad_%d",k);
          TPad *padtop = new TPad(thename,thename,0.,0.3,1.,1.);
          padtop->SetLeftMargin(padtop->GetLeftMargin()+0.08);
          padtop->SetBottomMargin(0);//padtop->GetBottomMargin()+0.08);
          padtop->SetTopMargin(0.02);//padtop->GetBottomMargin()+0.08);
          padtop->SetRightMargin(0.04);
          padtop->SetFillColor(0);
          padtop->SetFillStyle(0);
          padtop->SetTickx();
          padtop->SetTicky();
          padtop->Draw();
          padtop->cd();
          Zframes[k]->SetMinimum(1e-4);
          if(k<Z_bins) Zframes[k]->SetMaximum(Zframes[k]->GetMaximum()*max_scale);
          if(k>=Z_bins) Zframes[k]->SetMaximum(Zframes[k]->GetMaximum()*max_scale2);
          //Zframes[k]->SetMaximum(1.05*Zframes[k]->GetMaximum());
          Zframes[k]->SetTitle(rangelabels[(k % Z_bins)]);
          Zframes[k]->SetTitleFont(132,"t");
          Zframes[k]->GetXaxis()->SetLabelSize(0.09*0.78/0.7);
          Zframes[k]->GetXaxis()->SetTitleSize(0.09*0.78/0.7);
          Zframes[k]->GetYaxis()->SetTitleSize(0.09*0.78/0.7);
          TString thetitle=Zframes[k]->GetYaxis()->GetTitle();
          /*thetitle.Replace(10,1,"");
          if(k < Z_bins)thetitle.Replace(27,1,"");
          if(k >= Z_bins)thetitle.Replace(16,1,"");
          thitle.Replace(0,6,"Candidates");*/
          //Zframes[k]->GetYaxis()->SetTitle("");
          Zframes[k]->GetYaxis()->SetLabelSize(0.09*0.78/0.7);
          Zframes[k]->GetXaxis()->SetTitleOffset(0.95);
          Zframes[k]->GetYaxis()->SetTitleOffset(0.95);
          Zframes[k]->GetYaxis()->SetNdivisions(506);
          Zframes[k]->Draw();
          t->SetTextSize(0.07);
          t->SetTextAlign(33);
          t->SetTextAngle(90);
          c2->cd((k<Z_bins)*(2*k+1)+(k>=Z_bins)*(2*(k+1-Z_bins)));
          if(k >= Z_bins) {t->DrawLatex(0.01,0.99,thetitle);}
          if(k < Z_bins) {t->DrawLatex(0.01,0.99,thetitle);}
          t->SetTextAlign(22);
          t->SetTextAngle(0);
          padtop->cd();
          t->SetTextSize(0.09*0.78/0.7);
          if(k >= Z_bins) {t->DrawLatex(2250,Zframes[k]->GetMaximum()*0.92,"");}
          if(k < Z_bins) {t->DrawLatex(8.7,Zframes[k]->GetMaximum()*0.92,"");}


        }
      
    }
    cerr << data->sumEntries() << '\t' << model_hf->expectedEvents(RooArgSet(*x,*y,*z,*idx)) << endl;

 }
