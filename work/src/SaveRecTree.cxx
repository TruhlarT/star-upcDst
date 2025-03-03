#include "RecTree.h"
#include "StPicoPhysicalHelix.h"


void RecTree::SaveEventInfo(const StUPCEvent *upcEvt)
{
   this->setTofMult( upcEvt->getTOFMultiplicity());
   this->setEventNumber( upcEvt->getEventNumber());
   this->setFillNumber( upcEvt->getFillNumber());
   this->setBunchCrossId( upcEvt->getBunchCrossId());
   this->setBunchCrossId7bit( upcEvt->getBunchCrossId7bit());
   this->setRunNumber( upcEvt->getRunNumber() );
   this->setNVertecies( upcEvt->getNumberOfVertices() );
}

void RecTree::SaveTrackInfo(const StUPCTrack *trk, unsigned int iTrack)
{
   this->setDEdxInKevCm( trk->getDEdxSignal()*1000000 , iTrack); // convert to KeV/cm
   TVector3 momentum;
   trk->getMomentum(momentum);
   this->setMomentumInGev(  momentum.Mag() , iTrack );
   this->setPtInGev(  trk->getPt() , iTrack );
   this->setPxInGev(  momentum.X() , iTrack );
   this->setPyInGev(  momentum.Y() , iTrack );
   this->setPzInGev(  momentum.Z() , iTrack );
   this->setTofTimeInNs(  trk->getTofTime() , iTrack );
   this->setTofLengthInCm(  trk->getTofPathLength() , iTrack );
   this->setCharge(  trk->getCharge() , iTrack );
   this->setDcaXYInCm(  trk->getDcaXY() , iTrack );
   this->setDcaZInCm(  trk->getDcaZ() , iTrack );
   this->setNHitsFit(  trk->getNhitsFit() , iTrack );
   this->setNHitsDEdx(  trk->getNhitsDEdx() , iTrack );
   this->setEta(  trk->getEta() , iTrack );
   this->setPhi(  trk->getPhi() , iTrack );
   this->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kPion) , iTrack, PION );
   this->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kKaon) , iTrack, KAON );
   this->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kProton) , iTrack, PROTON );
}

void RecTree::SaveRPConfigInfo()
{ 
   double mRPDeltaPhi = TMath::Abs(mPhiRp[East] - mPhiRp[West])*convertToDegree;
   if(mRPDeltaPhi > 180)
      mRPDeltaPhi = 360 - mRPDeltaPhi;

   mIsElastic = mRPDeltaPhi > 90 ? true : false ; 
}//SaveRPConfigInfo



void RecTree::SaveRPinfo(const StUPCRpsTrack *trackRP, unsigned int iSide)
{
   const StUPCRpsTrackPoint* trackPoint = trackRP->getTrackPoint(RP1);

   this->setRpX( trackPoint->x() , iSide );
   this->setRpY( trackPoint->y() , iSide );
   this->setRpZ( trackPoint->z() , iSide );
   this->setThetaRp( trackRP->thetaRp(2) , iSide );
   this->setPhiRp( trackRP->phiRp() , iSide );
   this->setTimeRp( trackRP->time() , iSide );
   this->setT( trackRP->t(mUtil->p0()) , iSide );
   this->setPRp( trackRP->p() , iSide );
   this->setPtRp( trackRP->pt() , iSide );
   this->setEtaRp( trackRP->eta() , iSide );
   this->setXi( trackRP->xi(mUtil->p0()) , iSide );
   this->setPx( trackRP->pVec().X() , iSide );
   this->setPy( trackRP->pVec().Y() , iSide );
   this->setPz( trackRP->pVec().Z() , iSide );

/*
   // subtract offset
   int rnNumber, day, run, rpID;
   rpID = trackPoint->rpId();
   rnNumber = upcEvt->getRunNumber() - 18000000;
   day = rnNumber / 1000 - 55;
   run = rnNumber % 1000;
   rpYinner[iSide] = rpY[iSide] - mOffSet[rpID][day][run]; 
*/
}

void RecTree::SaveStateInfo(TLorentzVector state){

   this->setInvMass( state.M() );
   this->setTheta( state.Theta() );
   this->setPhi( state.Phi() );
   this->setP( state.P() );
   this->setPt( state.Pt() );
   this->setRap( state.Rapidity() );
}

void RecTree::SaveVertexInfo(const StUPCVertex *vtx)
{
   mVertexIsSet = true;
   this->setVertexZInCm( vtx->getPosZ() );
   this->setVertexYInCm( vtx->getPosY() );
   this->setVertexXInCm( vtx->getPosX() );
}


void RecTree::SaveZdcInfo(const StUPCEvent *upcEvt)
{

   for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
   {
      this->setZdcAdcEastPmt( upcEvt->getZDCEastADC(iPmt+1), iPmt);
      this->setZdcAdcWestPmt( upcEvt->getZDCWestADC(iPmt+1), iPmt);
   }
   this->setZdcTdcEast( upcEvt->getZDCEastTDC() );
   this->setZdcTdcWest( upcEvt->getZDCWestTDC() );
   this->setZdcTimeDiff( upcEvt->getZDCTimeDiff() );
   this->setZdcVertexZ( upcEvt->getZdcVertexZ() );
   this->setZdcEastRate( upcEvt->getZDCEastRate() );
   this->setZdcWestRate( upcEvt->getZDCWestRate() );
   this->setZdcUnAttEast( upcEvt->getZDCUnAttEast() );
   this->setZdcUnAttWest( upcEvt->getZDCUnAttWest() );

}

void RecTree::SaveBbcInfo(const StUPCEvent *upcEvt)
{
   this->setBbcSmallEast( upcEvt->getBBCSmallEast());
   this->setBbcSmallWest( upcEvt->getBBCSmallWest());
   this->setBbcLargeEast( upcEvt->getBBCLargeEast());
   this->setBbcLargeWest( upcEvt->getBBCLargeWest());
}

void RecTree::SetDeltaTrackInfo(const StUPCTrack *trk1, const StUPCTrack *trk2, unsigned int s)
{
   mDeltaPt[s] = trk1->getPt() - trk2->getPt(); 
   mDeltaEta[s] = trk1->getEta() - trk2->getEta(); 
   mDeltaPhi[s] = trk1->getPhi() - trk2->getPhi(); 
   mDeltaCurvature[s] = trk1->getCurvature() - trk2->getCurvature(); 
   mDeltaDipAngle[s] = trk1->getDipAngle() - trk2->getDipAngle(); 
   mDeltaPhase[s] = trk1->getPhase() - trk2->getPhase(); 
}

void RecTree::SetDeltaTrackInfo(unsigned int s)
{
   mDeltaPt[s] = -999; 
   mDeltaEta[s] = -999; 
   mDeltaPhi[s] = -999; 
   mDeltaCurvature[s] = -999; 
   mDeltaDipAngle[s] = -999; 
   mDeltaPhase[s] = -999; 
}


void RecTree::SaveTracksDcas(const StUPCTrack* particle1, const StUPCTrack* particle2, const StUPCEvent *upcEvt)
{
   pair<double,double> dca;
   StPicoPhysicalHelix p1Helix = StPicoPhysicalHelix (particle1->getCurvature(), 
                                                     particle1->getDipAngle(), 
                                                     particle1->getPhase(), 
                                                     particle1->getOrigin(),
                                                     particle1->getCharge());

   StPicoPhysicalHelix p2Helix =  StPicoPhysicalHelix (particle2->getCurvature(), 
                                                     particle2->getDipAngle(), 
                                                     particle2->getPhase(), 
                                                     particle2->getOrigin(),
                                                     particle2->getCharge());

   p1Helix.moveOrigin(p1Helix.pathLength(0,0));
   p2Helix.moveOrigin(p2Helix.pathLength(0,0));


/*
   pair<double, double> const ss = p1Helix.pathLengths(p2Helix);
   TVector3 const p1AtDcaToP2 = p1Helix.at(ss.first);
   TVector3 const p2AtDcaToP1 = p2Helix.at(ss.second);
*/
   TVector3 p1AtDcaToP2 = p1Helix.origin();
   TVector3 p2AtDcaToP1 = p2Helix.origin();

   TVector3 mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;

   double x = upcEvt->getBeamXPosition() + upcEvt->getBeamXSlope()*mDecayVertex[Z];
   double y = upcEvt->getBeamYPosition() + upcEvt->getBeamYSlope()*mDecayVertex[Z];
   TVector3 beam( x, y, mDecayVertex[Z]);

   p1Helix.moveOrigin(p1Helix.pathLength(beam));
   p2Helix.moveOrigin(p2Helix.pathLength(beam));

   p1AtDcaToP2 = p1Helix.origin();
   p2AtDcaToP1 = p2Helix.origin();

   this->setVertexStudyHypoZ(mDecayVertex[Z]);
   // calculate DCA of particle1 to particle2 at their DCA
   this->setVertexStudyDcaParticles((p1AtDcaToP2 - p2AtDcaToP1).Mag());
   // distance between the "vertex and beam line"
   this->setVertexStudyDcaBeamline((mDecayVertex - beam).Mag());
}


void RecTree::SaveTriggerInfo(const StUPCEvent *upcEvt, const StRPEvent *rpEvt)
{
   this->SaveBbcInfo(upcEvt);
   this->SaveZdcInfo(upcEvt);
   this->setBbcTrigBit( upcEvt->getBBCSmallEast() > 20 || upcEvt->getBBCSmallWest() > 20 || upcEvt->getBBCLargeEast() > 50 || upcEvt->getBBCLargeWest() > 50);

   this->setZdcETrigBit(false);
   this->setZdcWTrigBit(false);
   if( (upcEvt->getZDCEastTDC() > 500 && upcEvt->getZDCEastTDC() < 2700) && (upcEvt->getZDCEastADC(1) > 25 || upcEvt->getZDCEastADC(2) > 25 || upcEvt->getZDCEastADC(3) > 25))
   {
      this->setZdcTrigBit(true);      
      this->setZdcETrigBit(true);
   }   
   if( (upcEvt->getZDCWestTDC() > 500 && upcEvt->getZDCWestTDC() < 2700) && (upcEvt->getZDCWestADC(1) > 25 || upcEvt->getZDCWestADC(2) > 25 || upcEvt->getZDCWestADC(3) > 25))
   {
      this->setZdcTrigBit(true);
      this->setZdcWTrigBit(true);
   }    

   this->setTofTrigBit( !( upcEvt->getTOFMultiplicity() < 2 || upcEvt->getTOFMultiplicity() > 10) );
   saveRpTrigBit(rpEvt);

   bool RpEt = upcEvt->getLastDSM0() & (1 << 2);
   bool TofMult0 = upcEvt->getLastDSM0() & (1 << 4);
   bool TofMult2 = upcEvt->getLastDSM0() & (1 << 6);
   bool RpIt = upcEvt->getLastDSM0() & (1 << 12);
   bool BbcE = upcEvt->getLastDSM1() & (1 << 1);
   bool BbcW = upcEvt->getLastDSM1() & (1 << 2);
   bool BbcLE = upcEvt->getLastDSM1() & (1 << 3);
   bool BbcLW = upcEvt->getLastDSM1() & (1 << 4);
   bool ZdcE = upcEvt->getLastDSM1() & (1 << 7);
   bool ZdcW = upcEvt->getLastDSM1() & (1 << 10);

   this->setRpItDsmBit(RpIt);
   this->setRpEtDsmBit(RpEt);
   this->setTofDsmABit( TofMult0 );
   this->setTofDsmBBit( TofMult2 );
   this->setTofDsmBit( ( TofMult0 && !TofMult2 ) );
   this->setRpDsmBit( ( (RpEt && !RpIt) || (!RpEt && RpIt) ) );
   this->setBbcSmallEDsmBit(BbcE);
   this->setBbcSmallWDsmBit(BbcW);
   this->setBbcLargeEDsmBit(BbcLE);
   this->setBbcLargeWDsmBit(BbcLW);
   this->setBbcDsmBit( ( BbcE || BbcW || BbcLE || BbcLW ) );
   this->setZdcEDsmBit(ZdcE);
   this->setZdcWDsmBit(ZdcW);
   this->setZdcDsmBit( ( ZdcE || ZdcW ) );
}

void RecTree::saveRpTrigBit(const StRPEvent *rpEvt)
{
   for (unsigned int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      this->setRpTrigBits( IsRpTrigBit(rpEvt, iRp), iRp);
   }

   bool EU = this->getRpTrigBits(E1U) || this->getRpTrigBits(E2U);
   bool ED = this->getRpTrigBits(E1D) || this->getRpTrigBits(E2D);
   bool WU = this->getRpTrigBits(W1U) || this->getRpTrigBits(W2U);
   bool WD = this->getRpTrigBits(W1D) || this->getRpTrigBits(W2D);

   bool RP_ET = (EU && WD) || (ED && WU);
   bool RP_IT = (EU && WU) || (ED && WD);

   this->setRpItTrigBit(RP_IT);
   this->setRpEtTrigBit(RP_ET);
   this->setRpTrigBit((RP_ET && !RP_IT) || (!RP_ET && RP_IT));

}

bool RecTree::IsRpTrigBit(const StRPEvent *rpEvt, unsigned int iRp)
{
   for (unsigned int iPmt = 0; iPmt < 2; ++iPmt)
      if( rpEvt->tac(iRp, iPmt) > 200 && rpEvt->tac(iRp, iPmt) < 1750 && rpEvt->adc(iRp, iPmt) > 30)
         return true;

   return false;   
}