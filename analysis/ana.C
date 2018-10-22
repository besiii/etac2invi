#include "/afs/ihep.ac.cn/users/l/liukai/.rootInclude.h"

void ana(string filelist = "", char* nfname = "", bool isIncl = false){
	bool ifCut = true;


	double const EMASS = 0.511*0.001;
	double const MUMASS= 0.105658;
	double const PI0MASS= 0.134977;

	//new file
	TFile *newf=new TFile(nfname,"recreate");
	TTree* ntree=new TTree("ana","");

	TChain chain("ana","");

	string filename;
	ifstream is(filelist.c_str());
	while(getline(is,filename))
	{
		chain.Add(filename.c_str());
		cout<<" filename "<<filename<<endl;

	}
	cout<<"entries "<<chain.GetEntries()<<endl;






	Int_t run,event,ntrk, nshw, ncharged, nneutral;
	Double_t p4trk[4][6];//0,3 p4; 4 charge; 5 type;
	Double_t p4shw[4][4];
	Double_t p4truth[100][6];//0,3 p4;4 pdgid; 5 motherid;
	Int_t indexmc;

	chain.SetBranchAddress("run",&run);
	chain.SetBranchAddress("event",&event);
	chain.SetBranchAddress("ntrk",&ntrk);
	chain.SetBranchAddress("nshw",&nshw);
	chain.SetBranchAddress("p4trk",&p4trk);
	chain.SetBranchAddress("p4shw",&p4shw);
	chain.SetBranchAddress("ncharged",&ncharged);
	chain.SetBranchAddress("nneutral",&nneutral);

	if(isIncl){
		chain.SetBranchAddress("indexmc", &indexmc);
		chain.SetBranchAddress("p4truth",&p4truth);
	}




	Double_t n_mpipi, n_rmpipi;
	Double_t n_mpipigam, n_rmpipigam;
	Double_t n_egam, n_cosgam;

	Int_t n_indexmc, n_motheridx[100], n_pdgid[100];
	Double_t n_p4truth[100][6];//0,3 p4;4 pdgid; 5 motherid;

	ntree->Branch("egam", &n_egam,"egam/D");
	ntree->Branch("cosgam", &n_cosgam,"cosgam/D");
	ntree->Branch("mpipi", &n_mpipi,"mpipi/D");
	ntree->Branch("rmpipi", &n_rmpipi,"rmpipi/D");
	ntree->Branch("mpipigam", &n_mpipigam,"mpipigam/D");
	ntree->Branch("rmpipigam", &n_rmpipigam,"rmpipigam/D");

	if(isIncl){
		ntree->Branch("indexmc", &n_indexmc, "indexmc/I");
		ntree->Branch("motheridx", n_motheridx, "motheridx[100]/I");
		ntree->Branch("pdgid", n_pdgid, "pdgid[100]/I");
		ntree->Branch("p4truth",n_p4truth,"p4truth[100]/D");
	}

		bool debug = 1;



	TLorentzVector cms_p4(0.011*3.686 ,0,0, 3.686);
	for(Int_t i=0; i<chain.GetEntries(); i++)
	{
		chain.GetEntry(i);


		if(ntrk != 2) continue;
		if(nshw != 1) continue;
		if(ncharged != ntrk) continue;
		if(nneutral != nshw) continue;

		int Npip(0), Npim(0);
		TLorentzVector p4pip, p4pim;
		for(int a=0; a< ntrk; a++)
		{
			double charge = p4trk[a][4];
			if(charge >0){ Npip ++;   
				for(int b=0; b<4; b++) p4pip[b]	= p4trk[a][b];
			}
			if(charge <0) {
				Npim ++;
				for(int b=0; b<4; b++) p4pim[b]	= p4trk[a][b];
			}
		}
		if(Npip != 1 || Npim != 1) continue;




		double pipiAngle = p4pip.Angle(p4pim.Vect());
		double cospipi = TMath::Cos(pipiAngle);
		if(cospipi >= 0.95 ) continue;

		TLorentzVector p4pipi= p4pip + p4pim;
		n_mpipi = p4pipi.M();

		double cospipi_sys = fabs(p4pipi.CosTheta());
		if(cospipi_sys >=0.9) continue;


		double rmpipi = (cms_p4 - p4pip - p4pim).M();
		//3.08351  3.1107
		if(rmpipi<3.0835 || rmpipi>3.1107) continue;

		n_rmpipi = rmpipi;



		TLorentzVector p4gam;
		for(int c=0; c<4; c++) p4gam[c] = p4shw[0][c];
		n_egam = p4gam.E();
		n_cosgam = p4gam.CosTheta();

		n_mpipigam = (p4pip + p4pim + p4gam).M();
		n_rmpipigam= (cms_p4 - p4pip - p4pim - p4gam).M();









		if(isIncl){
			for(int c=0; c<indexmc; c++)
			{
				n_motheridx[c] = p4truth[c][5];
				n_pdgid[c] = p4truth[c][4];
				for(int d=0; d<6; d++) n_p4truth[c][d] = p4truth[c][4];
			}
			n_indexmc = indexmc;
		}
		ntree->Fill();

	}
	newf->cd();
	ntree->Write();
	newf->Close();
	cout<<"finished ^.^"<<endl;

}
