#ifndef EtacInv_Header
#define EtacInv_Header

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "McDecayModeSvc/McDecayModeSvc.h"
//#include "BestDTagSvc/BestDTagSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include <vector>
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
//using Event::McParticle;
#include "ParticleID/ParticleID.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/VertexFit.h"
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTrkPara;
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"

#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "DTagTool/DTagTool.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"


extern VertexFit * vtxfit;
extern KinematicFit * kmfit;
//
//namespace
//
class EtacInv:public Algorithm {
  public:
    EtacInv(const std::string& name, ISvcLocator* pSvcLocator);
    ~EtacInv();
    StatusCode initialize();
    StatusCode beginRun();   
    StatusCode execute();
    StatusCode endRun();
    StatusCode finalize();

int selectGoodChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent, SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,std::vector<int> & iGood) ;
bool passVertexSelection(CLHEP::Hep3Vector xorigin,	RecMdcKalTrack* mdcTrk ) ;
bool isGoodPion(EvtRecTrackIterator itTrk);
void selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent, SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,  vector<int> &iGam) ;
void mctruth();


  private:
    double m_vr0cut;
    double m_vz0cut;
	double m_cha_costheta_cut, m_total_number_of_charged_max, m_min_emctime,m_max_emctime, m_costheta_barrel_max, m_costheta_endcap_min,m_costheta_endcap_max, m_energy_barrel_min,m_energy_endcap_min,m_gammaCosCut,m_photon_iso_angle_min;
	bool m_isZCcondition;


    NTuple::Tuple* m_tuple1;
    NTuple::Item<int> m_run, m_event;
    NTuple::Item<int> m_idxmc;
    //NTuple::Array<int> m_pdgid, m_motheridx;
	NTuple::Item<int> m_ntrk, m_nshw, m_ncharged, m_nneutral;
	NTuple::Matrix<double> m_p4trk, m_p4shw, m_p4truth;


};

//add your inline methods
//?

#endif//EtacInv_Header
