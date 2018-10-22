#include "EtacInv/EtacInv.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
const int PSI2S_PDG_ID = 100443;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};//e,mu,pi,k,p

EtacInv::EtacInv(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name,pSvcLocator){
		declareProperty("Vr0cut",   m_vr0cut=1.0);
		declareProperty("Vz0cut",   m_vz0cut=10.0);
		declareProperty("cha_costheta_cut", m_cha_costheta_cut=0.93);
		declareProperty("total_number_of_charged_max", m_total_number_of_charged_max =10 );
		declareProperty("MinEstCut", m_min_emctime=0.0);
		declareProperty("MaxEstCut", m_max_emctime=14.0);
		declareProperty("CosthetaBarrelMax", m_costheta_barrel_max=0.8);
		declareProperty("CosthetaEndcapMin", m_costheta_endcap_min=0.86);
		declareProperty("CosthetaEndcapMax", m_costheta_endcap_max=0.92);
		declareProperty("ZChi_AnaCondition", m_isZCcondition=false);
		declareProperty("EnergyBarrelMin", m_energy_barrel_min=0.025); 
		declareProperty("EnergyEndcapMin", m_energy_endcap_min=0.050); 
		declareProperty("GammaCosCut",  m_gammaCosCut= 0.93); 
		declareProperty("PhotonIsoAngleMin", m_photon_iso_angle_min=20);

}
EtacInv::~EtacInv(){
	//add your code for deconstructor
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode EtacInv::initialize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"EtacInv::initialize()"<<endreq;
	//add your code here
	StatusCode status;
	NTuplePtr nt1(ntupleSvc(), "FILE1/ana");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/ana", CLID_ColumnWiseTuple, "ana N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("run",   m_run);
			status = m_tuple1->addItem ("event", m_event);
			status = m_tuple1->addItem("indexmc",          m_idxmc, 0, 100);
			//status = m_tuple1->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			//status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			status = m_tuple1->addIndexedItem("p4truth", m_idxmc, 6, m_p4truth);

			status = m_tuple1->addItem ("ntrk",  m_ntrk, 0, 10);//good charged track
			status = m_tuple1->addItem ("nshw",  m_nshw, 0, 4);//good neutrual track
			status = m_tuple1->addItem ("ncharged", m_ncharged, 0, 10);
			status = m_tuple1->addItem ("nneutral", m_nneutral, 0, 4);
			status = m_tuple1->addIndexedItem("p4trk",     m_ntrk,6, m_p4trk);
			status = m_tuple1->addIndexedItem("p4shw",     m_nshw,4, m_p4shw);
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode EtacInv::beginRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"EtacInv::beginRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode EtacInv::execute(){

	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"EtacInv::execute()"<<endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	m_run = eventHeader->runNumber();
	m_event = eventHeader->eventNumber();
	if(m_run%10000==0)log << MSG::WARNING <<"run, evtnum = " << m_run << " , " << m_event <<endreq;


	// Begin to select events
	// to retrieve RecEvent
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
	if ( ! evtRecEvent ) { 
		cout << MSG::FATAL << "Could not find EvtRecEvent" << endl;
		exit(1);
	}   
	// to retrieve RecTrackCol
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol( eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
	if ( ! evtRecTrkCol ) { 
		cout << MSG::FATAL << "Could not find EvtRecTrackCol" << endl;
		exit(1);
	}



	//charged
	int totalCharged = evtRecEvent->totalCharged();
	if (totalCharged< 2 || totalCharged >4) return StatusCode::SUCCESS;
	m_ncharged = totalCharged;
	std::vector<int> iGood;
	int ngoodtrk = selectGoodChargedTracks(evtRecEvent, evtRecTrkCol, iGood);
	if(ngoodtrk <2 || ngoodtrk>4) return StatusCode::SUCCESS;

	//pi
	int pipId = -2;
	int pimId = -2;
	int pipCount =0;
	int pimCount =0;
	for(int a = 0; a < iGood.size(); a++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[a];
		RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack(); 
		double p = mdcTrk->p();
		if(p<0.45)
		{
			if(!isGoodPion(itTrk)) continue;
			if(mdcTrk->charge() >0) { pipId = a;pipCount++;} 
			if(mdcTrk->charge() <0) { pimId = a;pimCount++;}
		}
	}
	if(pipCount ==0 || pimCount ==0) return StatusCode::SUCCESS;

	for(int a = 0; a < iGood.size(); a++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[a];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack(); 
		double mass = 0;
		if(a == pipId || a == pimId )
		{
			mdcKalTrk->setPidType(RecMdcKalTrack::pion);
			mass = xmass[2];
			m_p4trk[a][5]=2;//pid type
			for ( int ll = 0; ll < 4; ll++ ) m_p4trk[a][ll]=mdcKalTrk->p4(mass)[ll];
			m_p4trk[a][4]=mdcKalTrk->charge();
		}

	}

	m_ntrk=iGood.size();


	int nphoton = evtRecEvent->totalNeutral();
	if(nphoton ==0||nphoton>4) return StatusCode::SUCCESS;
	m_nneutral = nphoton;
	std::vector<int> iGoodGam;   iGoodGam.clear();
	selectNeutralTracks(evtRecEvent, evtRecTrkCol, iGoodGam);
	if ( iGoodGam.size() ==0 ||  iGoodGam.size()>4  ) return StatusCode::SUCCESS;

	for(int i=0;i<iGoodGam.size();i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGoodGam[i];
		RecEmcShower* emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();
		double phi = emcTrk->phi();
		double theta = emcTrk->theta();
		HepLorentzVector p4 = HepLorentzVector(eraw * sin(theta) * cos(phi), eraw * sin(theta) * sin(phi),eraw * cos(theta),eraw );
		for(int a=0; a<4; a++) m_p4shw[i][a] = p4[a];
	}
	m_nshw = iGoodGam.size();


	if(m_run <0) mctruth();

	m_tuple1->write();

	return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode EtacInv::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"EtacInv::endRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode EtacInv::finalize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"EtacInv::finalize()"<<endreq;
	//add your code here
	std::cout<<"00-00-00"<<endl;
	std::cout<<"save p4 of pi pi gam"<<endl;
	return StatusCode::SUCCESS;
}



//--------------------------------------------------------------------------------------------
//--------------------add your code here,for other member-functions---------------------------
//--------------------------------------------------------------------------------------------

int EtacInv::selectGoodChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
		SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
		std::vector<int> & iGood) {

	CLHEP::Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double *dbv = vtxsvc->PrimaryVertex(); 
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}

	iGood.clear();
	int np=0;
	int nm=0;

	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();

		if (!passVertexSelection(xorigin, mdcTrk) ) continue; 
		if(fabs(cos(mdcTrk->theta())) > m_cha_costheta_cut) continue;

		iGood.push_back((*itTrk)->trackId());

	} // end charged tracks
	return iGood.size();
}

bool EtacInv::passVertexSelection(CLHEP::Hep3Vector xorigin,
		RecMdcKalTrack* mdcTrk ) {
	HepVector a = mdcTrk->helix();
	HepSymMatrix Ea = mdcTrk->err();
	HepPoint3D point0(0.,0.,0.);

	VFHelix helixip(point0,a,Ea);
	HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	helixip.pivot(IP);
	HepVector vecipa = helixip.a();

	double m_vz0 = vecipa[3];
	double m_vr0 = vecipa[0];

	if(fabs(m_vz0) >= m_vz0cut) return false;
	if(fabs(m_vr0) >= m_vr0cut) return false;

	return true;
}

bool EtacInv::isGoodPion(EvtRecTrackIterator itTrk){
	//pi particle identification
	ParticleID * pidp = ParticleID::instance();
	pidp->init();
	pidp->setMethod(pidp->methodProbability());
	pidp->setChiMinCut(4);
	pidp->setRecTrack(*itTrk);
	// use PID sub-system
	pidp->usePidSys(pidp->useDedx() | pidp->useTof1() | pidp->useTof2());
	pidp->identify(pidp->onlyPionKaonProton());
	pidp->calculate();
	if(pidp->IsPidInfoValid()) {
		double probpi = pidp->probPion();
		double probk  = pidp->probKaon();
		if(probpi > probk && probpi > 0.001) return true;
		else
			return false;
	}
	else
		return false;
}



void EtacInv::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
		SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,  vector<int> &iGam) {


	for(int i=evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
		if (i > m_total_number_of_charged_max) break;

		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i ;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();

		// TDC window
		if ( !(emcTrk->time() >= m_min_emctime && emcTrk->time() <= m_max_emctime) )
			continue; 

		// Energy threshold
		double abs_costheta(fabs(cos(emcTrk->theta())));
		bool barrel = (abs_costheta < m_costheta_barrel_max); 
		bool endcap = (abs_costheta > m_costheta_endcap_min
				&& abs_costheta < m_costheta_endcap_max);
		double eraw = emcTrk->energy();

		if (!m_isZCcondition){     // Cut by "costheta"
			if ( !( (barrel && eraw > m_energy_barrel_min)
						|| (endcap && eraw > m_energy_endcap_min)))  continue; 
		}
		else{                      // Cut by "module"
			int module = emcTrk->module();
			if( module == 1 ){  if( !(eraw > m_energy_barrel_min) ) continue; }
			else{ if( !(eraw > m_energy_endcap_min) ) continue; }
		}

		// photon isolation: the opening angle between a candidate shower
		// and the closest charged track should be larger than 10 degree 
		CLHEP::Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		// EMC costheta cut 
		double costhe = cos(emcpos.theta());
		if ( fabs(costhe) >= m_gammaCosCut) continue;

		// find the nearest charged track
		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.; 
		for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			CLHEP::Hep3Vector extpos = extTrk->emcPosition();
			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);
			thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;

			if(fabs(thed) < fabs(dthe)) dthe = thed;
			if(fabs(phid) < fabs(dphi)) dphi = phid;
			if(angd < dang) dang = angd;	    
		}

		if(dang>=200) continue;
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);
		if (dang < m_photon_iso_angle_min ) continue; 

		iGam.push_back((*itTrk)->trackId());
	} // end loop neutral tracks     


}


void EtacInv::mctruth(){

	SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

	int m_numParticle = 0;
	bool Decay = false;
	int rootIndex = -1;
	bool strange = false;
	Event::McParticleCol::iterator iter_mc_topo = mcParticleCol->begin();
	for (; iter_mc_topo != mcParticleCol->end(); iter_mc_topo++){
		if ((*iter_mc_topo)->primaryParticle()&&(*iter_mc_topo)->particleProperty()==11&&((*iter_mc_topo)->mother()).particleProperty()==11) {strange=true;}
		if ((*iter_mc_topo)->primaryParticle()) continue;
		if (!(*iter_mc_topo)->decayFromGenerator()) continue;
		if ((*iter_mc_topo)->particleProperty() == PSI2S_PDG_ID){
			Decay = true;
			rootIndex = (*iter_mc_topo)->trackIndex();
		}
		if (!Decay) continue;
		int mcidx = ((*iter_mc_topo)->mother()).trackIndex() - rootIndex;
		int pdgid = (*iter_mc_topo)->particleProperty();
		if(strange&&((*iter_mc_topo)->mother()).particleProperty()!=PSI2S_PDG_ID) mcidx--;
		HepLorentzVector p4iter = (*iter_mc_topo)->initialFourMomentum();
		for(int a=0; a<4; a++) m_p4truth[m_numParticle][a] = p4iter[a];
		m_p4truth[m_numParticle][4] = pdgid * 1.0;
		m_p4truth[m_numParticle][5] = mcidx * 1.0;

		m_numParticle++;
	}
	m_idxmc = m_numParticle;

}


