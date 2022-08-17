
/*********************
This version has a file find option enables to read multiple files from a location
**********************/
//Reconstruction of D0->K- Pi+
#include <TLorentzVector.h>
#include <iostream>
#include <vector>
using namespace std;

void D1massrecreation1()  
{

	char *outputpath = "/home/siddharth/ATHENA/minQ2=1000/EICOutputFiles/OutputK-Pi+/";
	char *dirname="/home/siddharth/ATHENA/minQ2=1000/EICInputFiles/10x100/";
	char *ext=".root";
	TString o = "OUTPUT";
	TSystemDirectory dir(dirname, dirname);

	TList *files = dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next())) {
		  fname = file->GetName();
		  if (!file->IsDirectory() && fname.EndsWith(ext)) {
		     cout << "READING FILE : " << fname.Data() << endl;
		     TString fullfname = dirname + fname;
		     cout << "SOURCE DIRECTORY : " << fullfname.Data() << endl;
		     
				 TFile * myFile = new TFile(fullfname.Data());



	//Open File
	//TFile *myFile = new TFile("pythia8NCDIS_10x100_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0068.root");	//Open the root file
	TTreeReader myReader("events", myFile); 								//Pass the Tree's name and the TFile name in TTreeReader


	//Get Dynamical Variables
	TTreeReaderArray<Float_t> px(myReader, "ReconstructedParticles.p.x");
	TTreeReaderArray<Float_t> py(myReader, "ReconstructedParticles.p.y");
	TTreeReaderArray<Float_t> pz(myReader, "ReconstructedParticles.p.z");
	TTreeReaderArray<Float_t> energy(myReader, "ReconstructedParticles.energy");
	TTreeReaderArray<Int_t> pid(myReader, "ReconstructedParticles.pid");

	//Dmass Recreation: (range is from 0 to 6)
	TH1F *dmass0 = new TH1F("MD0","Combinational Mass with no conditions",100, 0, 6);		//Combinational Mass with no conditions
	TH1F *dmass1 = new TH1F("MD1","Combinational Mass with Pseudorapidity between -2 and 4; Mass (GeV/c^2) ;Events",100, 0, 6);
	TH1F *dmass2 = new TH1F("MD2","Combinational Mass with Pseudorapidity between -2 and 7 (Extended Pt Cut); Mass (GeV/c^2) ;Events",100, 0, 6);
	TH1F *dmass3 = new TH1F("MD3","Combinational Mass with Pt > 0.5; Mass (GeV/c^2) ;Events",100, 0, 6);
	
	TH1F *dpseudorapidity = new TH1F("DPseudo","Pseudorapidity of Combinations; Pseudorapidity ;Events",100, 0, 6);
	TH1F *dpt = new TH1F("DPt","pT of the Combinational Masses;Transverse Momentum (GeV/c);Events",100, 0, 25);
	
	//K:
	TH1F *khpx = new TH1F("Kx_Kaon","Momentum in X Direction;Px (GeV/c);Events",100, -15, 15);
	TH1F *khpy = new TH1F("Ky_Kaon","Momentum in Y Direction;Py (GeV/c);Events",100, -15, 15);
	TH1F *khpz = new TH1F("Kz_Kaon","Momentum in Z Direction;Pz (GeV/c);Events",100, -15, 15);
	TH1F *ken = new TH1F("Energy_Kaon","Energy of Kaons; Energy (GeV);Events",100, 0, 25);
	//TH1F *kpseudorapidity = new TH1F("Pseudorapidity_K","Pseudorapidity",100, -5, 5);
	TH1F *kpt = new TH1F("KPt","Transverse Momentum;Pt (GeV/c);Events",100, 0, 6);
	TH1F *k_mass = new TH1F("Mass_K","Invariant Mass Distribution of Kaon; Mass (GeV/c^2); Events",100, 0.49366, 0.49370);											//Mass of Kaon
	
	//Pi:
	TH1F *pihpx = new TH1F("Px_Pion","Momentum in X Direction;Px (GeV/c);Events",100, -10, 10);
	TH1F *pihpy = new TH1F("Py_Pion","Momentum in Y Direction;Py (GeV/c);Events",100, -10, 10);
	TH1F *pihpz = new TH1F("Pz_Pion","Momentum in Z Direction;Pz (GeV/c);Events",100, -10, 10);
	TH1F *pien = new TH1F("Energy_Pion","Energy; Energy (GeV);Events",100, 0, 25);
	//TH1F *pipseudorapidity = new TH1F("Pseudorapidity_Pi","Pseudorapidity",100, -5, 5);
	TH1F *pipt = new TH1F("PiPt","Transverse Momentum;Pt (GeV/c);Events",100, 0, 6);
	TH1F *pi_mass = new TH1F("Mass_Pi","Invariant Mass Distribution of Pion; Mass (GeV/c^2); Events",100, 0.139554, 0.139580);				//Mass of Pion
	
	double dmass, dpseudo, dptcut;
	
	//We are creating here TLorentzVector to store 4 vectors for d, k and pi
	TLorentzVector d;
	TLorentzVector pi;
	TLorentzVector k;

	//We will now store the 4vectors of pi and k in pivalues and kvalues respectively
	std::vector<TLorentzVector> *pivalues = new std::vector<TLorentzVector>();
	std::vector<TLorentzVector> *kvalues = new std::vector<TLorentzVector>();

	int p=0, q=0, r=0, s=0, t=0, u=0;		//counters

	while (myReader.Next())		//while loop to access events
	{
		p++;  //counters
		for (int i = 0; i < pid.GetSize(); i++)								//for loop from 0 to pid.GetSize()
		{
			if(pid[i] == 211)																	//recreate mass of Pi+ PID = 221
			{
				while (true)
				{
					r++;	//counter
					pi.SetPxPyPzE(px[i],py[i],pz[i],energy[i]);			//Creating the 4 vector for pi
					pihpx->Fill(pi.Px());														//Momentum in x direction
					pihpy->Fill(pi.Py());														//Momentum in y direction
					pihpz->Fill(pi.Pz());														//Momentum in z direction
					pien->Fill(pi.E());															//Energy 
					pi_mass->Fill(pi.M());													//Mass
					//pipseudorapidity->Fill(pi.PseudoRapidity());		//PseudoRapidity
					pipt->Fill(pi.Pt());
					pivalues->push_back(pi);												//pivalues stores the 4 vector of pi 
					break;
				}
			}
			
			
			else if(pid[i] == -321)															//recreate mass of K- PID = -321
			{
				while (true)
				{
					s++;	//counters
					k.SetPxPyPzE(px[i],py[i],pz[i],energy[i]);			//Creating the 4 vector for k
					khpx->Fill(k.Px());															//Momentum in x direction
					khpy->Fill(k.Py());															//Momentum in y direction
					khpz->Fill(k.Pz());															//Momentum in z direction
					ken->Fill(k.E());																//Energy
					k_mass->Fill(k.M());														//Mass
					//kpseudorapidity->Fill(pi.PseudoRapidity());			//PseudoRapidity
					kpt->Fill(k.Pt());
					kvalues->push_back(k);													//kvalues stores the 4 vector of k
					break;
				}
			}
		}
	}


	//Temporary variables xx and yy used to call the 4 vectors to add them
	TLorentzVector xx;
	TLorentzVector yy;

	int xcount =0, ycount=0;	//Creating variabkes to act as pointers for pivalues and kvalues
	
	//Nested for loops to create pairings of Pi and K
	for (auto x = pivalues->begin(); x != pivalues->end(); ++x)
	{	
		xx = (*pivalues)[xcount];
		for (auto y = kvalues->begin(); y != kvalues->end(); ++y)
		{
			yy = xx + (*kvalues)[ycount];		//Adding the 4 vectors of k and pi to give a pair of D
			dmass = yy.M();									//dmass is the mass of the pair of pi and k 
			
			dmass0->Fill(dmass);						//Filling the mass of the 4 vector calculated, gives a distribution of Combinations
			
			dpseudo = yy.PseudoRapidity();
			dpseudorapidity->Fill(dpseudo);	//Pseudorapidity of the combinational mass

			dptcut=yy.Pt();
			dpt->Fill(dptcut);							//Pt of the combinational mass
			
			if (dmass>=1.82 && dmass<=1.90)
			q++; //Counter to check for D Mass by pairing

			if(dpseudo<=4 && dpseudo>=-2 && dptcut>0.5)	//Final 
				dmass1->Fill(yy.M());
			if(dpseudo<=7 && dpseudo>=-2)
				dmass2->Fill(yy.M());
			if(dptcut>0.5)
				dmass3->Fill(yy.M());
			
				
			++ycount;
		}
		++xcount;
		ycount = 0;
	}


	cout << " p:" << p << " r: " << r << " s: " << s << endl;	//Counters
	cout << "Number of D0 Mass Recreations : " << q << endl;
	
	TFile * f = new TFile(fullfname.Data());
	TString outputfullfname = outputpath + fname;
	TFile * f1 = new TFile(outputfullfname.Data(),"recreate");

	
	pihpx->SetLineWidth(3); pihpx->Write();
	pihpy->SetLineWidth(3); pihpy->Write();
	pihpz->SetLineWidth(3); pihpz->Write();
	pien->SetLineWidth(3); pien->Write();
	pi_mass->SetLineWidth(3); pi_mass->Write();
	//pipseudorapidity->SetLineWidth(3); pipseudorapidity->Write();
	pipt->SetLineWidth(3); pipt->Write();
	
	khpx->SetLineWidth(3); khpx->Write();
	khpy->SetLineWidth(3); khpy->Write();
	khpz->SetLineWidth(3); khpz->Write();
	ken->SetLineWidth(3); ken->Write();
	k_mass->SetLineWidth(3); k_mass->Write();
	//kpseudorapidity->SetLineWidth(3); kpseudorapidity->Write();
	kpt->SetLineWidth(3); kpt->Write();
	
	dmass0->SetLineWidth(3); dmass0->Write();	
	dmass1->SetLineWidth(3); dmass1->Write();	
	dmass2->SetLineWidth(3); dmass2->Write();	
	dmass3->SetLineWidth(3); dmass3->Write();	
	
	dpseudorapidity->SetLineWidth(3); dpseudorapidity->Write();
	dpt->SetLineWidth(3); dpt->Write();
	
	gStyle->SetLineWidth(2);
	}
}
}
}
