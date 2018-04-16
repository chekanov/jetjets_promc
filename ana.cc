/***************************************************************************
 *  How to use ProMC files from HepSim, and how to build anti-KT jets 
 *  S.Chekanov (ANL) chekanov@anl.gov
 *  A library for HEP events storage and processing based on Google's PB   
 *  The project web site: http://atlaswww.hep.anl.gov/hepsim/
****************************************************************************/

#include<iostream>
#include<fstream>
#include<stdlib.h>

// check directory
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include "TMath.h"
#include"time.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <map>
struct stat sb;

// promc
#include "ProMC.pb.h"
#include "ProMCBook.h"
#include "LParticle.h"
#include "CParticle.h"

const double kPI   = TMath::Pi();
const double k2PI  = 2*kPI;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

using namespace std;
using namespace promc;


// find all files inside a directory
std::vector<std::string> open(std::string name = "data.in") {
	vector<std::string> ntup;
	ifstream myfile;
	myfile.open(name.c_str(), ios::in);

	if (!myfile) {
		cerr << " -> Can't open input file:  " << name << endl;
		exit(1);
	} else {
		cout << "-> Read data file=" << name << endl;
	}

	string temp;
	while (myfile >> temp) {
		//the following line trims white space from the beginning of the string
		temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
		if (temp.find("#") == 0) continue;
		ntup.push_back(temp);
	}
	cout << "-> Number of files=" << ntup.size()  << endl;
	myfile.close();

	for (unsigned int i=0; i<ntup.size(); i++) {
		cout << ".. file to analyse="+ntup[i] << endl;
	}
	return ntup;
}


// main example
int main(int argc, char **argv)
{

	string mass("-");
	if (argc == 2) {
		mass = argv[1];
	} else if (argc != 3) {
		cerr << "Usage: " << argv[0] << "[mass]" << endl;
		exit(1);
	}


	// fastjet
	const double ptLepton=60;
	const double ptJet=20;
	const double etaJet=2.4;
	const double R = 0.4;
	const double Rlep_iso=0.4;
	const double Rjet_iso=0.2;
	cout << "min PT lepton=" << ptLepton << endl;
	cout << "min PT jet=" << ptJet << endl;
	cout << "max ETA jet=" << etaJet << endl;
	cout << "lepton dR isolation=" << Rlep_iso << endl;
	cout << "jet    dR isolation=" << Rjet_iso << endl;


	int nDiJets=0;
	// total events
	int ntot=0; // total events
	int nfiles=0; // total files
	double weight=1.0;


	std::vector<std::string> files = open("data.in");

	double cross=0;
	double xcross=0;
	int nselect=0;

	string outputfile="output.root";
	cout << "\n -> Output file is =" << outputfile << endl;
	TFile * RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
	//  RootFile->SetCompressionLevel(0);

	TH1D * h_debug = new TH1D("debug", "debug", 10, 0, 10.);
	TH1D * h_cross = new TH1D("cross", "cross,events,lumi", 5, 0, 5.);
	TH1D * h_info = new TH1D("dijet_info", "Dijet info", 5, 0, 5.);
	TH1D * h_pt = new TH1D("jet_pt", "pt",100,0,1000);
	TH1D * h_eta = new TH1D("jet_eta", "eta", 40, -10, 10);
	TH1D * h_pt_lead = new TH1D("jet_pt_lead", "pt of lead jet",1000,0,30000);
	TH1D * h_pt_lepton = new TH1D("lepton_pt", "pt of leptons",200,0,2000);
	TH1D * h_n_lepton = new TH1D("lepton_nr", "Nr of leptons",20,0,20);
	TH1D * h_iso = new TH1D("iso_energy", "isolation fraction", 100, 0, 3.0);

	double mjjBins[] = {99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156, 13407, 13659, 13912, 14166, 14421, 14677, 14934, 15192, 15451, 15711, 15972, 16234, 16497, 16761, 17026, 17292, 17559, 17827, 18096, 18366, 18637, 18909, 19182, 19456, 19731, 20007, 20284, 20562, 20841, 21121, 21402, 21684, 21967, 22251, 22536, 22822, 23109, 23397, 23686, 23976, 24267, 24559, 24852, 25146, 25441, 25737, 26034, 26332, 26631, 26931, 27232, 27534, 27837, 28141, 28446, 28752, 29059, 29367, 29676, 29986};

	const int nBins=sizeof(mjjBins)/sizeof(double);
	double xbins[nBins];
	double xbins_tev[nBins];
	for (int j=0; j<nBins; j++){xbins[j]=mjjBins[j]; xbins_tev[j]=0.001*mjjBins[j]; };
	TH1D* jetjetmass_2jet=new TH1D( "JetJetMass_2jet", "JetJet Mass > 1 jet", nBins-1, xbins);jetjetmass_2jet->Sumw2();
	TH1D * binsM = new TH1D("bins_m", "bins_m", nBins-1, xbins);
	binsM->Sumw2();
	TH1D * binsM_tev = new TH1D("bins_m_tev", "bins_m_tev", nBins-1, xbins_tev);
	binsM_tev->Sumw2();


	for (Int_t j=0; j<nBins-1; j++) {
		float x=xbins[j+1]-xbins[j];
		binsM->Fill(xbins[j]+0.5*x,x);
		float xx=xbins_tev[j+1]-xbins_tev[j];
		binsM_tev->Fill(xbins_tev[j]+0.5*xx,xx);
		// cout << j << " " << xbins[j+1] << " bin size=" << x << endl;
	}
	// set bin errors to 0 (we use it to divide the bin width!)
	for (int i=0 ; i<(binsM->GetNbinsX()); i++) {
		binsM->SetBinError(  i+1, 0.0);
		binsM_tev->SetBinError(  i+1, 0.0);
	}


	// jets
	Strategy strategy = fastjet::Best;
	JetDefinition jet_def(fastjet::antikt_algorithm, R, strategy);
	//JetDefinition jet_def(fastjet::kt_algorithm, Rparam, strategy);
	//JetDefinition jet_def(fastjet::cambridge_algorithm, Rparam, strategy);
	//JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);


	std::map<int,int> chargemap; // pad that keeps charge*3

	for(unsigned int m=0; m < files.size(); m++){
		string Rfile=files[m];
		ProMCBook*  epbook = new ProMCBook(Rfile.c_str(),"r");

		cout << "\n\n Start to read.." << endl;
		// get the version number
		int  h=epbook->getVersion();
		if (m==0) cout << "Version = " << h << endl;
		// get the description of this file
		string d=epbook->getDescription();
		if (m==0) cout << "Description = " << d << endl;
		int  nev=epbook->getEvents();
		cout << "Events = " << nev  << endl;
		// get the header file with units, cross sections etc.
		ProMCHeader header = epbook->getHeader();


		// this part reads the header information with particle data.
		// you can access names, id, charge*3, and other useful
		// information from PDG table. As an example, we make a map
		// that keeps charge*3.
		if (m==0) {
			for (int jj=0; jj<header.particledata_size(); jj++){
				ProMCHeader_ParticleData data= header.particledata(jj);
				int charge3=data.charge();
				int id=data.id();
				double mass=data.mass();
				string name=data.name();
				cout << "name=" << name << " mass=" << mass << " charge=" << charge3 << endl;
				chargemap[id]=charge3;
			}
		}


		// here are the units
		double kEV=(double)(header.momentumunit());
		//double kLe=(double)(header.lengthunit());


		// loop over all events
		for (int nn=0; nn<nev; nn++){
			if (epbook->next() !=0) continue;
			ProMCEvent eve = epbook->get();

			// get truth information
			ProMCEvent_Particles  *pa=eve.mutable_particles();
			h_debug->Fill("Events",1.0);
			ntot++;

			double xlumi=0;
			if (m>0) {
				xlumi=(double)ntot/cross; // lumi so far
			}

			if (ntot%1000==0)
				cout <<  " # Events=" << ntot  << " X-cross="<< cross << " pb" << " Lumi so far=" << xlumi/1000.0 << " fb-1" << endl;


			vector<PseudoJet> avec;
			vector<int> nrid;
			vector<LParticle> lepton_candidates;
			// fill stable and no neutrino
			for (int j=0; j<pa->pdg_id_size(); j++){
				if (pa->status(j)!=1) continue;
				int type=pa->pdg_id(j);
				if (abs(type)==12 || abs(type)==14 || abs(type)==16 ) continue;
				double px= pa->px(j)/kEV;
				double py= pa->py(j)/kEV;
				double pz= pa->pz(j)/kEV;
				double ee= pa->energy(j)/kEV;
				double pt=sqrt(px*px+py*py);
				double eta=-log(tan(atan2(pt,(double)pz)/2));


				// muon or electron
				if (abs(type) ==11 || abs(type) ==13) {
					int charge=1;
					if (type<0) charge=-1;
					LParticle p(px,py,pz,ee,charge);
					p.SetStatus(j);
					p.SetType(type);
					if (pt>ptLepton && TMath::Abs(eta)<etaJet) lepton_candidates.push_back(p);
				};
				// int charge=chargemap[pa->pdg_id(j)]; // get charge
				if ( pt < 0.1)                   continue;
				if ( fabs(eta)> 3.5 )            continue;
				avec.push_back(  PseudoJet(px,py,pz,ee)  );
				nrid.push_back(j);
			}


			if (lepton_candidates.size()<1) continue;
			h_debug->Fill("Lepton cand.",1.0);


			// isolate leptons
			vector<LParticle> leptons;
			double IsoEnergy=0.1;
			// isolate leading lepton
			for (unsigned int k1 = 0; k1<lepton_candidates.size(); k1++) {
				LParticle xlep=(LParticle)lepton_candidates.at(k1);
				TLorentzVector lep=xlep.GetP();
				double p_pt=lep.Perp();
				double p_eta=lep.PseudoRapidity();
				double p_phi=lep.Phi();
				if (p_phi<0) p_phi=k2PI+p_phi;
				if (p_pt<ptLepton) continue;
				double esumP=0;
				for (unsigned int k2 = 0; k2<avec.size(); k2++) {
					PseudoJet part = avec.at(k2);
					double pt=part.perp();
					double eta=part.pseudorapidity();
					double phi=part.phi();
					if (phi<0) phi=k2PI+phi;
					double deta    = p_eta - eta;
					double dphi    = p_phi - phi;
					double adphi=TMath::Abs(dphi);
					if (adphi>kPI) adphi=k2PI-adphi;
					double ddr = TMath::Sqrt(deta*deta+adphi*adphi);
					if (ddr<Rlep_iso) esumP=esumP+pt;
				}
				double isoFrac=esumP/p_pt;
				//cout << "Esum=" << esumP << " p_pt = " << p_pt << " iso=" << isoFrac << endl;
				h_iso->Fill(isoFrac);
				if (isoFrac<1.0+IsoEnergy) {
					leptons.push_back(xlep);
				};
			}


			unsigned int nLeptons=leptons.size();
			if (nLeptons>1) std::sort(leptons.begin(), leptons.end(), greater<LParticle>() ) ;
			if (nLeptons<1) continue;
			LParticle pL=leptons.at(0);
			TLorentzVector LL=pL.GetP();
			h_pt_lepton->Fill(LL.Perp());
			h_n_lepton->Fill(nLeptons);
			h_debug->Fill("Iso leptons",1.0);


			// remove isolated lepton from vectors used for jets
			vector<PseudoJet> hadrons;
			for (unsigned int k = 0; k<avec.size(); k++) {
				PseudoJet part = avec.at(k);
				int id=nrid.at(k);
				int isLep=false;
				for (unsigned int ll=0; ll<leptons.size(); ll++){
					LParticle LL=leptons.at(ll);
					int id_lep=LL.GetStatus();
					if (id_lep == id) isLep=true;
				}
				if (!isLep) hadrons.push_back(part);
			}


			// make jets
			ClusterSequence clust_seq(hadrons, jet_def);
			vector<PseudoJet> jets_truth = clust_seq.inclusive_jets(ptJet);
			vector<PseudoJet> sorted_jets = sorted_by_pt(jets_truth);
			vector<LParticle> jets;

			for (unsigned int k = 0; k<sorted_jets.size(); k++) {
				double eta=sorted_jets[k].pseudorapidity();
				if ( fabs(eta)> etaJet )            continue;
				double phi=sorted_jets[k].phi();
				double pt = sorted_jets[k].perp();
				double e = sorted_jets[k].e();
				h_pt->Fill(pt);
				h_eta->Fill(eta);
				TLorentzVector l;
				l.SetPtEtaPhiE(pt, eta, phi, e);
				LParticle p;
				p.SetP(l);
				p.SetType(1);
				jets.push_back(p);

			} // end loop


			unsigned int nJets=jets.size();
			if (nJets>1) std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
			if (nJets<2) continue;
			h_debug->Fill("2 jets.",1.0);


			// overlap removal
			unsigned int nMaxJet=2;
			if (jets.size()<nMaxJet) nMaxJet=jets.size();


			vector<LParticle> leptons_iso;
			for (unsigned int ll=0; ll<nLeptons; ll++){
				LParticle LL=leptons.at(ll);
				TLorentzVector LP=LL.GetP();
				double phi_lep=LP.Phi();
				double eta_lep=LP.PseudoRapidity();
				double y_lep=LP.Rapidity();

				bool found=false; // isolated from 2 leading jets
				for (unsigned int ii=0; ii<nMaxJet; ii++){
					LParticle LPP=jets.at(ii);
					TLorentzVector LP=LPP.GetP();
					double phi_jet=LP.Phi();
					double eta_jet=LP.PseudoRapidity();
					double y_jet=LP.Rapidity();
					double deta=TMath::Abs(y_jet-y_lep);
					double dphi=TMath::Abs(phi_jet - phi_lep);
					if (dphi>kPI) dphi=k2PI-dphi;
					double dR=TMath::Sqrt(deta*deta+dphi*dphi);
					if (dR<Rjet_iso) found=true;
				}
				if (found) continue;
				leptons_iso.push_back(LL);
			}


			unsigned int nLeptonsIso=leptons_iso.size();
			if (nLeptonsIso<1) continue;
			h_debug->Fill("dR=0.4 Lep-jet ",1.0);


			// consider more than 2 jets
			if (jets.size()>1) {
				h_debug->Fill("2 jets",1.);
				LParticle p1=jets.at(0);
				LParticle p2=jets.at(1);
				TLorentzVector LP1=p1.GetP();
				h_pt_lead->Fill(LP1.Perp(),weight);
				TLorentzVector LP2=p2.GetP();
				TLorentzVector PP=LP1+LP2;
				double mass_jj=PP.M();
				jetjetmass_2jet->Fill(mass_jj,weight);
				nDiJets++;
			}


		} // end event loop


		ProMCStat stat = epbook->getStatistics();
		cross=stat.cross_section_accumulated();
		epbook->close(); // close
		nfiles++;
		xcross=xcross+cross;


	} // end loop over all files



	xcross=xcross/(double)nfiles; // average cross for all files
	cout << "Total events=" << ntot << endl;
	cout << "Total files=" << nfiles << endl;
	cout << "Total events =" << ntot << endl;

	cout << "Average cross section for all files = " << xcross << " pb"<< endl;
	double width=h_pt->GetBinWidth(1);
	double lumi=(double)ntot/xcross;
	cout << "Lumi for all files  =" << lumi << " pb-1" << endl;
	double norm=width*lumi;
	h_pt->Scale(1.0/norm);

	h_cross->Fill(0.0,(double)xcross); // 1
	h_cross->Fill(1.0,(double)ntot); // 2
	h_cross->Fill(2.0,(double)lumi); // 3
	h_cross->Fill(3.0,(double)nfiles);
	cout << "Nr of selected jets=" << nselect << endl;

	h_info->Fill("Nr of dijets",nDiJets); // calibration check
	h_info->Fill("Dijet section",nDiJets/lumi); // calibration check
	cout << " Nr of dijets=" << nDiJets << endl;
	cout <<" Observed cross section=" << nDiJets/lumi << " pb " << endl;

	ofstream myfile;
	myfile.open ("cross.txt");
	myfile <<  mass << " " << nDiJets/lumi << endl;
	myfile.close();

	/*
	// calculate cross sections
	   TAxis *xaxis = h_pt->GetXaxis(); 
	   xaxis->SetTitle("p_{T}(jet) GeV");
	   TAxis *yaxis = h_pt->GetYaxis();
	   yaxis->SetTitle("d #sigma / dp_{T} [pb / GeV]");

	    norm=lumi*(h_pt->GetBinWidth(1));
	    h_pt->Scale(1.0/norm); 
	    TAxis *xaxis1 = h_pt->GetXaxis();
	    xaxis1->SetTitle("p_{T}(jet) GeV");
	    TAxis *yaxis1 = h_pt->GetYaxis();
	    yaxis1->SetTitle("d #sigma / dp_{T} [pb / GeV]");
	*/

	//m_ntuple->Fill();
	RootFile->Write();
	RootFile->Print();
	RootFile->Close();

	cout << "Writing ROOT file "+ outputfile << endl;

	return 0;
}
