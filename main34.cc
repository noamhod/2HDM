// For this example to run, MadGraph 5 must be installed and the
// command "exe" (set by default as "mg5_aMC") must be available via
// the command line. Additionally, GZIP support must be enabled via
// the "--with-gzip" configuration option(s). Note that this example has
// only been tested with MadGraph 5 version 2.3.3; due to rapid
// MadGraph development, this example may not work with other
// versions. For more details on the LHAMadgraph class see the
// comments of Pythia8Plugins/LHAMadgraph.h.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/LHAMadgraph.h"

#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <sstream>

#define _FAT(x)                    { cout << "FATAL: " << __FILE__ << " " << __LINE__ << ": " << x << endl; exit(-1); }
#define _ERR(enable,x) if(enable==1) cout << "ERROR: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _WRN(enable,x) if(enable==1) cout << "WARNING: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _DBG(enable,x) if(enable==1) cout << "DEBUG: " << __FILE__ << " " << __LINE__ << ": " << x << endl;
#define _INF(enable,x) if(enable==1) cout << "INFO: "  << __FILE__ << " " << __LINE__ << ": " << x << endl;

using namespace Pythia8;

enum proc
{
	SM,SM_8TeV,MG_SM,SMNLO,SM_nobmass, SM_nobmass_8TeV,SMH,SMIA,SMIA_8TeV,SMIH,SMIHjj,IH,I,H,IA,A,A_8TeV,SMIA_valid,SMIH_valid
};

//==========================================================================
TFile* file;
TTree* tree;
vector<TLorentzVector>* p4;
vector<int>*            id;
float                   Q;
float                   Q2;
float                   aS;

void initFile(TString name)
{
	file = new TFile("tops."+name+".root","RECREATE");
	tree = new TTree(name,name+" gg->tt");
}
void setTree(int name)
{	
	gROOT->ProcessLine(".L Loader.C+"); // for the vector branches...
	
	p4 = new vector<TLorentzVector>;
	id = new vector<int>;
	
	switch(name)
	{
		case SM  :       initFile("SM");         break;
		case SM_8TeV:    initFile("SM_8TeV");    break;
		case MG_SM:      initFile("MG_SM");      break;
		case SMNLO:      initFile("SMNLO");      break;
		case SM_nobmass: initFile("SM_nobmass"); break;
		case SM_nobmass_8TeV: initFile("SM_nobmass_8TeV"); break;
		case SMH :       initFile("SMH");        break;
		case SMIA:       initFile("SMIA");       break;
		case SMIA_8TeV:  initFile("SMIA_8TeV");  break;
		case SMIH:       initFile("SMIH");       break;
		case SMIHjj:     initFile("SMIHjj");     break;
		case IH:         initFile("IH");         break;
		case I:          initFile("I");          break;
		case H:          initFile("H");          break;
		case IA:         initFile("IA");         break;
		case A:          initFile("A");          break;
		case A_8TeV:     initFile("A_8TeV");     break;
		case SMIA_valid: initFile("SMIA_valid"); break;
		case SMIH_valid: initFile("SMIH_valid"); break;
		default:   cout << "ENUM unkown: " << name << endl;  exit(-1);
	}	
	tree->Branch("p4",&p4);
	tree->Branch("id",&id);
	tree->Branch("Q",&Q);
	tree->Branch("Q2",&Q2);
	tree->Branch("aS",&aS);
}
void clearTree()
{
	p4->clear();
	id->clear();
	Q  = -1;
	Q2 = -1;
	aS = -1;
}
void fillTree()
{
	file->cd();
	tree->Fill();
}
void writeTree()
{
	file->cd();
	tree->Write();
}


//==========================================================================

// A simple method to run Pythia, analyze the events

void run(Pythia* pythia, int name, int nEvent) {
  pythia->readString("Random:setSeed = on");
  pythia->readString("Random:seed = 1");
  pythia->init();

  setTree(name);

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia->next()) continue;
    /////////////////  
	clearTree(); ////
	/////////////////
	
	Q  = pythia->info.QRen();
	Q2 = pythia->info.Q2Ren();
	aS = pythia->info.alphaS();
	
    int iGlu1(0), iGlu2(0);
	int iTop1(0), iTop2(0);
	// cout << "-----------------" << endl;
    for (int i = 0; i < pythia->process.size(); ++i) {
		// cout << i << " --> id=" << pythia->process[i].id() << " (" << pythia->process[i].e() << "," << pythia->process[i].px() << "," << pythia->process[i].py() << "," << pythia->process[i].pz() << ")" << endl;
      if (!iTop1 && pythia->process[i].id() ==  6) iTop1 = i;
      if (!iTop2 && pythia->process[i].id() == -6) iTop2 = i;
	  if (!iGlu1 && pythia->process[i].id() == 21) iGlu1 = i;
      if (!iGlu2 && pythia->process[i].id() == 21 && i!=iGlu1) iGlu2 = i;
      if (iTop1 && iTop2 && iGlu1 && iGlu2) {
        iTop1 = pythia->process[iTop1].iBotCopyId();
        iTop2 = pythia->process[iTop2].iBotCopyId();
		iGlu1 = pythia->process[iGlu1].iTopCopyId();
        iGlu2 = pythia->process[iGlu2].iTopCopyId();
		TLorentzVector pp1, pp2, pp3, pp4;
		pp1.SetPxPyPzE(pythia->process[iGlu1].px(),pythia->process[iGlu1].py(),pythia->process[iGlu1].pz(),pythia->process[iGlu1].e());
		pp2.SetPxPyPzE(pythia->process[iGlu2].px(),pythia->process[iGlu2].py(),pythia->process[iGlu2].pz(),pythia->process[iGlu2].e());
		pp3.SetPxPyPzE(pythia->process[iTop1].px(),pythia->process[iTop1].py(),pythia->process[iTop1].pz(),pythia->process[iTop1].e());
		pp4.SetPxPyPzE(pythia->process[iTop2].px(),pythia->process[iTop2].py(),pythia->process[iTop2].pz(),pythia->process[iTop2].e());
		p4->push_back(pp1);
		p4->push_back(pp2);
		p4->push_back(pp3);
		p4->push_back(pp4);
		id->push_back(pythia->process[iGlu1].id());
		id->push_back(pythia->process[iGlu2].id());
		id->push_back(pythia->process[iTop1].id());
		id->push_back(pythia->process[iTop2].id());
		
		// cout << "tops(" << iTop1 << "," << iTop2 << ")  gluons(" << iGlu1 << "," << iGlu2 << ")" << endl;
		// cout << "m(tt)=" << (pp3+pp4).M() << "  m(gg)=" << (pp1+pp2).M() << endl;
		
        break;
      }
    }
    fillTree();
  }
  writeTree();
  pythia->stat();
}

//==========================================================================

int main(int argc, char* argv[])
{
	if(argc!=2) _FAT("need to specify one arg");
	string name = argv[1];
	_INF(1,"arg="<<name);
	
	// The name of the MadGraph5_aMC@NLO executable.
	// You must prepend this string with the path to the executable
	// on your local installation, or otherwise make it available.
	string mgpath  = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/";
	string exe     = mgpath+"bin/mg5_aMC";
	string model   = mgpath+"models/Higgs_Effective_Couplings_FormFact/";
	string command = "/bin/cp -f "+model+"parameters.py_ORIG  "+model+"parameters.py";
	system(command.c_str());
	
	// Pythia object to be used a couple of times
	Pythia* pythia;
	int nEvents = 405000;
	string sEvents = "";
	stringstream strm;
	strm << nEvents;
	strm >> sEvents;
	
	bool dynamical_scale_choice3 = false;

	if(name=="SM")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the 2HDM model.
		pythia = new Pythia();
		system("rm -rf madgraphrun1");
		LHAupMadgraph madgraph1(pythia, true, "madgraphrun1", exe);
		// madgraph1.setEvents(nEvents);
		madgraph1.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph1.readString("generate g g > t t~ / h h1 HIW=0 HIG=0 QED=99 QCD=99");
		// Note the need for a blank character before "set".
		madgraph1.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph1.readString(" set dynamical_scale_choice 3");
		madgraph1.readString(" set MT 172.5");
		madgraph1.readString(" set nevents "+sEvents);
		madgraph1.readString(" set ebeam1 6500");
		madgraph1.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph1);
		run(pythia, SMNLO, nEvents);
		delete pythia;
	}
	if(name=="SM_8TeV")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the 2HDM model.
		pythia = new Pythia();
		system("rm -rf madgraphrun_SM_8TeV");
		LHAupMadgraph madgraph_SM_8TeV(pythia, true, "madgraphrun_SM_8TeV", exe);
		// madgraph_SM_8TeV.setEvents(nEvents);
		madgraph_SM_8TeV.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph_SM_8TeV.readString("generate g g > t t~ / h h1 HIW=0 HIG=0 QED=99 QCD=99");
		// Note the need for a blank character before "set".
		madgraph_SM_8TeV.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_SM_8TeV.readString(" set dynamical_scale_choice 3");
		madgraph_SM_8TeV.readString(" set MT 172.5");
		madgraph_SM_8TeV.readString(" set nevents "+sEvents);
		madgraph_SM_8TeV.readString(" set ebeam1 4000");
		madgraph_SM_8TeV.readString(" set ebeam2 4000");
		pythia->setLHAupPtr(&madgraph_SM_8TeV);
		run(pythia, SM_8TeV, nEvents);
		delete pythia;
	}
	else if(name=="MG_SM")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the MG shipped SM model.
		pythia = new Pythia();
		system("rm -rf madgraphrun_mgsm");
		LHAupMadgraph madgraph_mgsm(pythia, true, "madgraphrun_mgsm", exe);
		// madgraph_mgsm.setEvents(nEvents);
		madgraph_mgsm.readString("import model sm");
		madgraph_mgsm.readString("generate g g > t t~");
		// Note the need for a blank character before "set".
		madgraph_mgsm.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_mgsm.readString(" set dynamical_scale_choice 3");
		madgraph_mgsm.readString(" set MT 172.5");
		madgraph_mgsm.readString(" set nevents "+sEvents);
		madgraph_mgsm.readString(" set ebeam1 6500");
		madgraph_mgsm.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph_mgsm);
		run(pythia, MG_SM, nEvents);
		delete pythia;
	}
	else if(name=="SMNLO")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the MG shipped SM-loop model.
		pythia = new Pythia();
		system("rm -rf madgraphrun_nlo");
		LHAupMadgraph madgraph_nlo(pythia, true, "madgraphrun_nlo", exe);
		// madgraph_nlo.setEvents(nEvents);
		madgraph_nlo.readString("import model loop_sm");
		madgraph_nlo.readString("generate g g > t t~ [QCD]");
		// Note the need for a blank character before "set".
		madgraph_nlo.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_nlo.readString(" set dynamical_scale_choice 3");
		madgraph_nlo.readString(" set MT 172.5");
		madgraph_nlo.readString(" set nevents "+sEvents);
		madgraph_nlo.readString(" set ebeam1 6500");
		madgraph_nlo.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph_nlo);
		run(pythia, SMNLO, nEvents);
		delete pythia;
	}
	else if(name=="SM_nobmass")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the MG shipped SM-loop model.
		pythia = new Pythia();
		system("rm -rf madgraphrun_smnobmass");
		LHAupMadgraph madgraph_smnobmass(pythia, true, "madgraphrun_smnobmass", exe);
		// madgraph_smnobmass.setEvents(nEvents);
		madgraph_smnobmass.readString("import model sm-no_b_mass");
		madgraph_smnobmass.readString("generate g g > t t~");
		// madgraph_smnobmass.readString("define j = g u c d s b u~ c~ d~ s~ b~");
		// madgraph_smnobmass.readString("generate g g > t t~ @0");
		// madgraph_smnobmass.readString("add process g g > t t~ j @1");
		// madgraph_smnobmass.readString("add process g g > t t~ j j @2");
		// Note the need for a blank character before "set".
		madgraph_smnobmass.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_smnobmass.readString(" set dynamical_scale_choice 3");
		madgraph_smnobmass.readString(" set MT 172.5");
		madgraph_smnobmass.readString(" set nevents "+sEvents);
		madgraph_smnobmass.readString(" set ebeam1 6500");
		madgraph_smnobmass.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph_smnobmass);
		run(pythia, SM_nobmass, nEvents);
		delete pythia;
	}
	else if(name=="SM_nobmass_8TeV")
	{
		// Produce leading-order gg->tt SM events with MadGraph 5 using the MG shipped SM-loop model.
		pythia = new Pythia();
		system("rm -rf madgraphrun_smnobmass_8TeV");
		LHAupMadgraph madgraph_smnobmass_8TeV(pythia, true, "madgraphrun_smnobmass_8TeV", exe);
		// madgraph_smnobmass.setEvents(nEvents);
		madgraph_smnobmass_8TeV.readString("import model sm-no_b_mass");
		madgraph_smnobmass_8TeV.readString("generate g g > t t~");
		// madgraph_smnobmass_8TeV.readString("define j = g u c d s b u~ c~ d~ s~ b~");
		// madgraph_smnobmass_8TeV.readString("generate g g > t t~ @0");
		// madgraph_smnobmass_8TeV.readString("add process g g > t t~ j @1");
		// madgraph_smnobmass_8TeV.readString("add process g g > t t~ j j @2");
		// Note the need for a blank character before "set".
		madgraph_smnobmass_8TeV.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_smnobmass_8TeV.readString(" set dynamical_scale_choice 3");
		madgraph_smnobmass_8TeV.readString(" set MT 172.5");
		madgraph_smnobmass_8TeV.readString(" set nevents "+sEvents);
		madgraph_smnobmass_8TeV.readString(" set ebeam1 4000");
		madgraph_smnobmass_8TeV.readString(" set ebeam2 4000");
		pythia->setLHAupPtr(&madgraph_smnobmass_8TeV);
		run(pythia, SM_nobmass_8TeV, nEvents);
		delete pythia;
	}
	else if(name=="SMIA")
	{
		// Produce leading-order gg->tt SM+A events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun2");
		LHAupMadgraph madgraph2(pythia, true, "madgraphrun2", exe);
		// madgraph2.setEvents(nEvents);
		madgraph2.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph2.readString("generate g g > t t~ / h HIW=1 HIG=1 QED=99 QCD=99");
		// Note the need for a blank character before "set".
		madgraph2.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph2.readString(" set dynamical_scale_choice 3");
		madgraph2.readString(" set MT 172.5");
		madgraph2.readString(" set MP 500");
		madgraph2.readString(" set WH1 143.907354");
		madgraph2.readString(" set YMT 431.25");
		madgraph2.readString(" set YMB 1.88");
		madgraph2.readString(" set YMC 3.55");
		madgraph2.readString(" set YMTAU 0.7108");
		madgraph2.readString(" set YMM 0.042264");
		madgraph2.readString(" set nevents "+sEvents);
		madgraph2.readString(" set ebeam1 6500");
		madgraph2.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph2);
		run(pythia, SMIA, nEvents);
		delete pythia;
	}
	else if(name=="SMIA_8TeV")
	{
		// Produce leading-order gg->tt SM+A events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun_SMIA_8TeV");
		LHAupMadgraph madgraph_SMIA_8TeV(pythia, true, "madgraphrun_SMIA_8TeV", exe);
		// madgraph2.setEvents(nEvents);
		madgraph_SMIA_8TeV.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph_SMIA_8TeV.readString("generate g g > t t~ / h HIW=1 HIG=1 QED=99 QCD=99");
		// Note the need for a blank character before "set".
		// [0] tanb=0.400000 sba=1.000000 cba=0.000000 wA=143.907354 wH=80.405335 YMT=431.250000 YMB=1.880000 YMC=3.550000 YMM=0.042264 YMTAU=0.710800
        madgraph_SMIA_8TeV.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_SMIA_8TeV.readString(" set dynamical_scale_choice 3");
		madgraph_SMIA_8TeV.readString(" set MT 172.5");
		madgraph_SMIA_8TeV.readString(" set MP 500");
		madgraph_SMIA_8TeV.readString(" set WH1 143.907354");
		madgraph_SMIA_8TeV.readString(" set YMT 431.25");
		madgraph_SMIA_8TeV.readString(" set YMB 1.88");
		madgraph_SMIA_8TeV.readString(" set YMC 3.55");
		madgraph_SMIA_8TeV.readString(" set YMTAU 0.7108");
		madgraph_SMIA_8TeV.readString(" set YMM 0.042264");
		madgraph_SMIA_8TeV.readString(" set nevents "+sEvents);
		madgraph_SMIA_8TeV.readString(" set ebeam1 4000");
		madgraph_SMIA_8TeV.readString(" set ebeam2 4000");
		pythia->setLHAupPtr(&madgraph_SMIA_8TeV);
		run(pythia, SMIA_8TeV, nEvents);
		delete pythia;
	}
	else if(name=="A_8TeV")
	{
		// Produce leading-order gg->A->tt events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun_A_8TeV");
		LHAupMadgraph madgraph_A_8TeV(pythia, true, "madgraphrun_A_8TeV", exe);
		// madgraph_A_8TeV.setEvents(nEvents);
		madgraph_A_8TeV.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph_A_8TeV.readString("generate g g > h1 > t t~ / h QED=99 QCD=99");
		// Note the need for a blank character before "set".
		// [0] tanb=0.400000 sba=1.000000 cba=0.000000 wA=143.907354 wH=80.405335 YMT=431.250000 YMB=1.880000 YMC=3.550000 YMM=0.042264 YMTAU=0.710800
        madgraph_A_8TeV.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph_A_8TeV.readString(" set dynamical_scale_choice 3");
		madgraph_A_8TeV.readString(" set MT 172.5");
		madgraph_A_8TeV.readString(" set MP 500");
		madgraph_A_8TeV.readString(" set WH1 143.907354");
		madgraph_A_8TeV.readString(" set YMT 431.25");
		madgraph_A_8TeV.readString(" set YMB 1.88");
		madgraph_A_8TeV.readString(" set YMC 3.55");
		madgraph_A_8TeV.readString(" set YMTAU 0.7108");
		madgraph_A_8TeV.readString(" set YMM 0.042264");
		madgraph_A_8TeV.readString(" set nevents "+sEvents);
		madgraph_A_8TeV.readString(" set ebeam1 4000");
		madgraph_A_8TeV.readString(" set ebeam2 4000");
		pythia->setLHAupPtr(&madgraph_A_8TeV);
		pythia->particleData.addParticle(9000006,"h1",0,0,0,500.,143.907354);
		run(pythia, A_8TeV, nEvents);
		delete pythia;
	}
	else if(name=="A")
	{
		// Produce leading-order gg->tt A events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun3");
		LHAupMadgraph madgraph3(pythia, true, "madgraphrun3", exe);
		// madgraph3.setEvents(nEvents);
		madgraph3.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph3.readString("generate g g > h1 > t t~ / h QED=99 QCD=99");
		// [59] tanb=0.860000 sba=1.000000 cba=0.000000 wA=30.931584 wH=17.434106 YMT=200.581395 YMB=4.042000 YMC=1.651163 YMM=0.090868 YMTAU=1.528220
		// Note the need for a blank character before "set".
		madgraph3.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph3.readString(" set dynamical_scale_choice 3");
		madgraph3.readString(" set MT 172.5");
		madgraph3.readString(" set MP 500");
		madgraph3.readString(" set WH1 143.907354");
		madgraph3.readString(" set YMT 431.25");
		madgraph3.readString(" set YMB 1.88");
		madgraph3.readString(" set YMC 3.55");
		madgraph3.readString(" set YMTAU 0.7108");
		madgraph3.readString(" set YMM 0.042264");
		madgraph3.readString(" set nevents "+sEvents);
		madgraph3.readString(" set ebeam1 6500");
		madgraph3.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph3);
		pythia->particleData.addParticle(9000006,"h1",0,0,0,500.,143.907354);
		run(pythia,A,nEvents);
		delete pythia;
	}
	else if(name=="H")
	{
		// Produce leading-order gg->tt H events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun4");
		LHAupMadgraph madgraph4(pythia, true, "madgraphrun4", exe);
		// madgraph4.setEvents(nEvents);
		madgraph4.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph4.readString("generate g g > h > t t~ QED=99 QCD=99");
		// Note the need for a blank character before "set".
		madgraph4.readString(" set cut_decays F");
		if(dynamical_scale_choice3) madgraph4.readString(" set dynamical_scale_choice 3");
		madgraph4.readString(" set MT 172.5");
		madgraph4.readString(" set MH 500");
		madgraph4.readString(" set WH 49.63");
		madgraph4.readString(" set YMTAU 1.1351");
		madgraph4.readString(" set YMT 260.8092");
		madgraph4.readString(" set YMB 2.6825");
		madgraph4.readString(" set nevents "+sEvents);
		madgraph4.readString(" set ebeam1 6500");
		madgraph4.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph4);
		run(pythia,H,nEvents);
		delete pythia;
	}
	else if(name=="SMIH_valid")
	{
		// Produce leading-order gg->tt SM+A events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun_valid");
		LHAupMadgraph madgraph_valid(pythia, true, "madgraphrun_valid", exe);
		// madgraph_valid.setEvents(nEvents);
		madgraph_valid.readString("import model Higgs_Effective_Couplings_FormFact");
		madgraph_valid.readString("generate g g > t t~ / h1 HIW=1 HIG=1 QED=99 QCD=99");
		// Note the need for a blank character before "set".
		// [20] tanb=0.660000 sba=1.000000 cba=0.000000 wA=52.510675 wH=29.593605 YMT=-261.363636 YMB=3.102000 YMC=-2.151515 YMM=0.069736 YMTAU=1.172820
        madgraph_valid.readString(" set cut_decays F");
        if(dynamical_scale_choice3) madgraph_valid.readString(" set dynamical_scale_choice 3");
		madgraph_valid.readString(" set MT 172.5");
		madgraph_valid.readString(" set MP 500.");
		madgraph_valid.readString(" set MH 500.");
		madgraph_valid.readString(" set WH 29.593605");
		madgraph_valid.readString(" set WH1 52.510675");
		madgraph_valid.readString(" set YMTAU 1.172820");
		madgraph_valid.readString(" set YMM 0.069736");
		madgraph_valid.readString(" set YMC -2.151515");
		madgraph_valid.readString(" set YMB 3.102000");
		madgraph_valid.readString(" set YMT -261.363636"); 
		madgraph_valid.readString(" set nevents "+sEvents);
		madgraph_valid.readString(" set ebeam1 6500");
		madgraph_valid.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph_valid);
		run(pythia, SMIH_valid, nEvents);
		delete pythia;
	}
	else if(name=="SMIHjj")
	{
		// Produce leading-order gg->tt SM+A events with MadGraph 5.
		pythia = new Pythia();
		system("rm -rf madgraphrun_SMIHjj");
		LHAupMadgraph madgraph_SMIHjj(pythia, true, "madgraphrun_SMIHjj", exe);
		// madgraph_valid.setEvents(nEvents);
		madgraph_SMIHjj.readString("import model Higgs_Effective_Couplings_FormFact");
		// madgraph_SMIHjj.readString("generate g g > t t~ / h1 HIW=1 HIG=1 QED=99 QCD=99");
		madgraph_SMIHjj.readString("define j = g u c d s b u~ c~ d~ s~ b~");
		madgraph_SMIHjj.readString("define p = g u c d s b u~ c~ d~ s~ b~");
		madgraph_SMIHjj.readString("generate p p > t t~ / h1 HIW=1 HIG=1 QED=99 QCD=99 @0");
		madgraph_SMIHjj.readString("add process p p > t t~ j / h1 HIW=1 HIG=1 QED=99 QCD=99 @1");
		madgraph_SMIHjj.readString("add process p p > t t~ j j / h1 HIW=1 HIG=1 QED=99 QCD=99 @2");
		// Note the need for a blank character before "set".
		// [20] tanb=0.660000 sba=1.000000 cba=0.000000 wA=52.510675 wH=29.593605 YMT=-261.363636 YMB=3.102000 YMC=-2.151515 YMM=0.069736 YMTAU=1.172820
        madgraph_SMIHjj.readString(" set cut_decays F");
        if(dynamical_scale_choice3) madgraph_SMIHjj.readString(" set dynamical_scale_choice 3");
		madgraph_SMIHjj.readString(" set MT 172.5");
		madgraph_SMIHjj.readString(" set MP 500.");
		madgraph_SMIHjj.readString(" set MH 500.");
		madgraph_SMIHjj.readString(" set WH 29.593605");
		madgraph_SMIHjj.readString(" set WH1 52.510675");
		madgraph_SMIHjj.readString(" set YMTAU 1.172820");
		madgraph_SMIHjj.readString(" set YMM 0.069736");
		madgraph_SMIHjj.readString(" set YMC -2.151515");
		madgraph_SMIHjj.readString(" set YMB 3.102000");
		madgraph_SMIHjj.readString(" set YMT -261.363636"); 
		madgraph_SMIHjj.readString(" set nevents "+sEvents);
		madgraph_SMIHjj.readString(" set ebeam1 6500");
		madgraph_SMIHjj.readString(" set ebeam2 6500");
		pythia->setLHAupPtr(&madgraph_SMIHjj);
		run(pythia, SMIHjj, nEvents);
		delete pythia;
	}
	else
	{
		_FAT("unknown model");
	}
	
	return 0;
}