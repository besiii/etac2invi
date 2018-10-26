{
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.5);
	TChain chain("ana");
	chain.Add("../../selection/12data/data_cut.root");
	TChain chainmc("ana");
	chainmc.Add("../../selection/12incl/cut.root");

	//1
	TH1F* hmass = new TH1F("hmass","",50,2.9,3.05);
	chain.Project("hmass","rmpipigam","");
	cout<<"entries of hmass"<<hmass->GetEntries()<<endl;

	TH1F* hmassmc = new TH1F("hmassmc","",50,2.9,3.05);
	chainmc.Project("hmassmc","rmpipigam","");
	cout<<"entries of hmass"<<hmassmc->GetEntries()<<endl;

	//2
	TH1F* h1 = new TH1F("h1","",50,-1,1);
	chain.Project("h1","cosgam","");

	TH1F* h1mc = new TH1F("h1mc","",50,-1,1);
	chainmc.Project("h1mc","cosgam","");
	//3
	TH1F* h2 = new TH1F("h2","",50,3.07,3.13);
	chain.Project("h2","rmpipi","");
	TH1F* h2mc = new TH1F("h2mc","",50,3.07,3.13);
	chainmc.Project("h2mc","rmpipi","");


	TCanvas* can = new TCanvas("c","fit Dmass to get Cuts",1250,800);
	can->Divide(2,2);
	//Double_t Nmcnorm = hmassmc->GetEntries() * (341./400);

	can->cd(1);
	Double_t Nmcnorm = hmass->GetEntries();
	hmassmc->SetLineColor(2);
	hmassmc->SetLineWidth(4);
	hmassmc->GetXaxis()->SetTitle("RM(#pi#pi#gamma) GeV");
	hmassmc->DrawNormalized("e",Nmcnorm);
	hmass->SetLineWidth(2);
	hmass->Draw("same");


	can->cd(2);
	Nmcnorm = h1->GetEntries();
	h1mc->SetLineColor(2);
	h1mc->SetLineWidth(4);
	h1mc->GetXaxis()->SetTitle("cos(#gamma)");
	h1->GetXaxis()->SetTitle("cos(#gamma)");
	h1->SetLineWidth(2);
	h1->Draw();
	h1mc->DrawNormalized("esame",Nmcnorm);

	can->cd(3);
	Nmcnorm = h2->GetEntries();
	h2mc->SetLineColor(2);
	h2mc->SetLineWidth(4);
	h2mc->GetXaxis()->SetTitle("RM(#pi#pi) GeV");
	h2mc->DrawNormalized("e",Nmcnorm);
	h2->SetLineWidth(2);
	h2->Draw("same");


	can->Print("incldatashow.eps");




}
