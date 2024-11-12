#include "EmbedMaker.h"
#include "RunDef.h"

EmbedMaker::EmbedMaker() : kMaxPitchesToMatch(3.0) 
{ 
  // --- default values for pp run17 (in case it's not read from database) --- 
  mXYZ_IP[kX] = 0.0;
  mXYZ_IP[kY] = 0.0;
  mXYZ_IP[kZ] = 0.0;
  mThetaXY_tilt[kX] = 0.0;
  mThetaXY_tilt[kY] = 0.0;
  mDistanceFromIPtoDX[kEast] = 9.8;
  mDistanceFromIPtoDX[kWest] = 9.8;
  mLDX[kEast] = 3.7;
  mLDX[kWest] = 3.7;
  mBendingAngle[kEast] = 0.018832292;
  mBendingAngle[kWest] = 0.018826657;
  mConversion_TAC_time = 18e-12;
  mDoEmbedding = true;
  mDoAfterburner = true;
  // --- --- --- --- --- --- --- --- --- 
}


EmbedMaker::~EmbedMaker()
{
  //destructor

}//~EmbedMaker

void EmbedMaker::MakeTracks(float blue_beamenergy, float yellow_beamenergy) {
  vector<unsigned int> trackPointsVec[kMAXSEQ];
  // reconstructing track-points
  formTrackPoints( trackPointsVec );
  // reconstructing tracks
  formTracks( trackPointsVec, blue_beamenergy, yellow_beamenergy );
  // embed signal in PMTs in additive way ADC from MC + ADC from ZB
  MakePMTs();
}

void EmbedMaker::MakePMTs(){
  for(UInt_t iRP = 0; iRP < kMAXSEQ; ++iRP){
    double adc[2];
    double tac[2];
    adc[0] = MCEvent->adc(iRP, 0); adc[1] = MCEvent->adc(iRP, 1);
    tac[0] = MCEvent->tac(iRP, 0); tac[1] = MCEvent->tac(iRP, 1);
    if(mDoEmbedding){
      adc[0] += ZBEvent->adc(iRP, 0); adc[1] += ZBEvent->adc(iRP, 1);
      tac[0] += ZBEvent->tac(iRP, 0); tac[1] += ZBEvent->tac(iRP, 1);
    }
    rpEvent->setAdc(iRP, adc[0], adc[1]);
    rpEvent->setTac(iRP, tac[0], tac[1]); 
  }
}


void EmbedMaker::formTracks( const vector<unsigned  int > *trackPointVec, const float beamMomentumWest, const float beamMomentumEast ) const
{

  double beamMomentum[2];
  beamMomentum[kEast] = beamMomentumEast;
  beamMomentum[kWest] = beamMomentumWest;

  for(int branch=0; branch<kBranches; ++branch)
  { // loop over all branches in the Roman Pot system

    unsigned int side = ( branch < kBranches/2 ? kEast : kWest );
    int sign = (side == kEast ? -1 : 1 );
    int nPts[kStationsPerBranch]; // reading number of track-points found in the branch
    nPts[kRP1] = trackPointVec[ kRpInBranch[branch][kRP1] ].size();
    nPts[kRP2] = trackPointVec[ kRpInBranch[branch][kRP2] ].size();

    if( nPts[kRP1] && nPts[kRP2] )
    { // if track-points reconstructed in both stations in branch
      for(int i=0; i<nPts[kRP1]; ++i)
      { // loops over all combinations of track-points
        for(int j=0; j<nPts[kRP2]; ++j)
        {
          StUPCRpsTrack* track = rpEvent->addTrack();

          track->setEvent( rpEvent );
          track->setBranch( branch ); // setting ID of branch
          track->setType( StUPCRpsTrack::rpsGlobal ); // setting the type of the track
          track->setFirstTrackPointId( trackPointVec[ kRpInBranch[branch][kRP1] ][i]); // setting constituent track-points
          track->setSecondTrackPointId( trackPointVec[ kRpInBranch[branch][kRP2] ][j]); // setting constituent track-points
          //StUPCRpsTrackPoint *TP1 = rpEvent->getTrackPoint(trackPointVec[ kRpInBranch[branch][kRP1] ][i]);
          //StUPCRpsTrackPoint *TP2 = rpEvent->getTrackPoint(trackPointVec[ kRpInBranch[branch][kRP2] ][j]);
          
          // below calculating momentum vector
          StUPCRpsTrackPoint *trkPoint = rpEvent->getTrackPoint(trackPointVec[ kRpInBranch[branch][kRP1] ][i]);
          double localThetaX = track->thetaRp( StUPCRpsTrack::rpsAngleThetaX ) - sign*mThetaXY_tilt[kX]; // REMINDER: sensitive to changes in StUPCRpsTrack::thetaRp() !
          double localThetaY = track->thetaRp( StUPCRpsTrack::rpsAngleThetaY ) - sign*mThetaXY_tilt[kY]; // REMINDER: sensitive to changes in StUPCRpsTrack::thetaRp() !
          double x_BCS =  trkPoint->x() - sin(mThetaXY_tilt[kX])*( trkPoint->z() ); 
          // assuming vertex in [0,0,0]
          double d2 = abs( trkPoint->z() ) - mLDX[side] - mDistanceFromIPtoDX[side]; // distance from DX magnet exit to first RP station
          double thetaX_IP = ( x_BCS - (d2 + 0.5*mLDX[side])*localThetaX ) / ( mDistanceFromIPtoDX[side]  + 0.5*mLDX[side] );
          double xi = 1. / ( 1 + (mBendingAngle[side]*(mDistanceFromIPtoDX[side] + 0.5*mLDX[side])) / ( localThetaX*abs( trkPoint->z() ) - x_BCS ) );
          double momentumValue = beamMomentum[side] * (1.-xi);
          TVector3 momentumVector( 0, 0, sign*momentumValue );
          momentumVector.RotateX( -sign*localThetaY );
          momentumVector.RotateY( sign*thetaX_IP );
          track->setP( momentumVector ); // setting the momentum vector
        }
      }
    }
  }
}


void EmbedMaker::formTrackPoints(vector<unsigned  int > *trackPointVec) const
{
  for(int iRP=0; iRP<kMAXSEQ; ++iRP)
  { // loop over all Roman Pots
    int nTrackPoints = 0;
    // looking for hits in X and Y direction (necessary to determine (x,y) coordinates of track-point)
    vector<EmbedMaker::StRpsHit> hits[kCoordinates];

    hits[kY] = formHits(iRP, kY);
    if( hits[kY].size()==0 ) 
      continue; // if no hits in planes A&C => cannot reconstruct a track-point

    hits[kX] = formHits(iRP, kX);
    if( hits[kX].size()==0 ) 
      continue; // if no hits in planes B&D => cannot reconstruct a track-point

    // calculating time of detection in PMTs (invoked here to avoid multiple calclation in case of many hits)
    double time[2] = {-1, -1};
    for(unsigned int pmt=0; pmt<2; ++pmt)
    {
      if( ZBEvent->tac(iRP,pmt) < kMaxPedestalTAC ) 
        continue; // don't calculate time if TAC is at pedestal

      time[pmt] = timeFromTAC( iRP, pmt, ZBEvent->tac(iRP,pmt), ZBEvent->adc(iRP,pmt));
    }
    hHitSizeMap[iRP]->Fill(hits[kX].size(), hits[kY].size());
    // loops over all combinations of hits in X and Y directions
    for(unsigned int j=0; j<hits[kX].size(); ++j)
    {
      for(unsigned int k=0; k<hits[kY].size(); ++k)
      {
        StUPCRpsTrackPoint* trackPoint = rpEvent->addTrackPoint();

        // setting position of the track-point
        double x = hits[kX][j].mPositionXY;
        double y = hits[kY][k].mPositionXY;
        double z = (hits[kX][j].mPositionZ*hits[kX][j].mPlanesUsed + hits[kY][k].mPositionZ*hits[kY][k].mPlanesUsed) / (hits[kX][j].mPlanesUsed + hits[kY][k].mPlanesUsed);
        //cout<<"Creating TP: ["<<x<<" , "<<y<<" , "<<z<<" ]"<<endl;
        TVector3 pos( x, y, z );
        trackPoint->setPosition( pos );

        // setting ID of Roman Pot
        trackPoint->setRpId( iRP );

        // setting IDs of clusters used to form a track-point
        for(int l=0; l<kPlanesPerCoordinate; ++l)
        {
          trackPoint->setClusterId( hits[kY][k].mClusterId[l], kPlanes[kY][l] );
          trackPoint->setClusterId( hits[kX][j].mClusterId[l], kPlanes[kX][l] );
        }
        // setting time of the hit (in time units)
        for(unsigned int pmt=0; pmt<trackPoint->mNumberOfPmtsInRp; ++pmt)
        {
          trackPoint->setTime( time[pmt], pmt );
        }

        // setting flag of track-point quality
        if( hits[kX][j].mGolden && hits[kY][k].mGolden ) 
          trackPoint->setQuality( StUPCRpsTrackPoint::rpsGolden );
        else 
          trackPoint->setQuality( StUPCRpsTrackPoint::rpsNormal );

        // push back trackPoint ID to the vector for the current branch
        trackPointVec[iRP].push_back( rpEvent->getNumberOfTrackPoints() - 1);
        if(trackPoint->planesUsed() > 2)
          nTrackPoints++;
      }
    }
    hGoodTrackPoints[iRP]->Fill(nTrackPoints);
  }
}


vector<EmbedMaker::StRpsHit> EmbedMaker::formHits(const unsigned int RpId, 
  const int coordinate) const
{
  vector<EmbedMaker::StRpsHit> hitVec;
  vector<double> pos[kPlanesPerCoordinate];
  vector<unsigned int> en[kPlanesPerCoordinate];
  vector<unsigned int> len[kPlanesPerCoordinate];
  vector<unsigned int> id[kPlanesPerCoordinate];

  preselectClusters(RpId, coordinate, pos, en, len, id);
  int clCase = classifyClustersCase(pos);
  if(clCase>0)
  {

    std::vector<unsigned int> validClusters[kPlanesPerCoordinate];
    bool matched = matchClusters(coordinate, clCase, pos, validClusters);

    if( matched )
    {  // if there are pair of clusters which match - use only those
      for(unsigned int k=0; k<validClusters[kFirst].size(); ++k)
      {
        EmbedMaker::StRpsHit hit;
        hit.mPositionXY = ( pos[kFirst][validClusters[kFirst][k]] + pos[kSecond][validClusters[kSecond][k]] )/2;
        hit.mPositionZ = ( ZBEvent->z( RpId, kPlanes[coordinate][0]) + ZBEvent->z( RpId, kPlanes[coordinate][1]))/2;
        for(int j=0; j<kPlanesPerCoordinate; ++j)  
          hit.mClusterId[j] = id[j][validClusters[j][k]];
          
        hit.mPlanesUsed = kPlanesPerCoordinate;

        if(clCase==5) hit.mGolden = true; // golden hit <-- 1/1
        else hit.mGolden = false;

        hitVec.push_back( hit );
      }
    }else{ // if clusters don't match, use each one separately
      for(int j=0; j<kPlanesPerCoordinate; ++j)
      { // loop over 2 planes in given _coordinate_
        for(unsigned int k=0; k<pos[j].size(); ++k)
        {
          EmbedMaker::StRpsHit hit;
          hit.mPositionXY = pos[j][k];
          hit.mPositionZ = ZBEvent->z( RpId, kPlanes[coordinate][j]);
          for(int l=0; l<kPlanesPerCoordinate; ++l)
          {
            if(l==j) hit.mClusterId[l] = id[l][k];
            else hit.mClusterId[l] = -1;
          }
          hit.mPlanesUsed = 1;
          hit.mGolden = false;
          hitVec.push_back( hit );
        }
      }
    }
  }

  for (unsigned int i = 0; i < hitVec.size(); ++i)
  {
    EmbedMaker::StRpsHit hit = hitVec[i];
    for (int l = 0; l < 2; ++l)
      if(hit.mClusterId[l] > -1)
      {
        StUPCRpsCluster *cluster = rpEvent->getCluster(hit.mClusterId[l]);

        eneSimHits->Fill(cluster->energy(),4*cluster->romanPotId() + cluster->planeId());
        posSimHits->Fill(hit.mPositionXY,4*cluster->romanPotId() + cluster->planeId());    
      }
  }

  return hitVec;
}


void EmbedMaker::preselectClusters(const unsigned int RpId, const int coordinate, 
  vector<double>* pos, vector<unsigned int>* en, vector<unsigned int>* len, vector<unsigned int>* id) const
{
  vector<int> clusterPerPlane[kMAXSEQ][kMAXCHAIN]; // 8 sequencers/roman pots, 4 chains/planes
  vector<int> clusterPerPlaneUsed[kMAXSEQ][kMAXCHAIN]; // 8 sequencers/roman pots, 4 chains/planes

  StRPEvent *tmpEvt = new StRPEvent(*ZBEvent);
  tmpEvt->clearEvent();

  int clusterID = 0;
  for (unsigned int iClstr = 0; iClstr < MCEvent->getNumberOfClusters(); ++iClstr)
  {
    StUPCRpsCluster *MCCluster = MCEvent->getCluster(iClstr);
    if(!MCCluster)
      continue;

    if(MCCluster->romanPotId() != RpId) 
      continue;
    
    if(!(MCCluster->planeId() == kPlanes[coordinate][0] || MCCluster->planeId() == kPlanes[coordinate][1])) 
      continue;

    StUPCRpsCluster *cluster = tmpEvt->addCluster();
    cluster->setPosition( MCCluster->position() );
    cluster->setLength(MCCluster->length()); 
    cluster->setEnergy(MCCluster->energy()); 
    cluster->setXY(MCCluster->xy()); 
    cluster->setPlaneId(MCCluster->planeId());
    cluster->setRomanPotId(RpId);
    cluster->setQuality(0);
    clusterPerPlane[RpId][MCCluster->planeId()].push_back(clusterID++);
  }
  if( mDoEmbedding )
  {
    for (unsigned int iClstr = 0; iClstr < ZBEvent->getNumberOfClusters(); ++iClstr)
    {
      StUPCRpsCluster *ZBCluster = ZBEvent->getCluster(iClstr);
      if(ZBCluster->romanPotId() != RpId) 
        continue;

      if(ZBCluster->planeId() != kPlanes[coordinate][0] && ZBCluster->planeId() != kPlanes[coordinate][1]) 
        continue;

      int len = ZBCluster->length();
      int plane = ZBCluster->planeId();
      double pitch = STRIP_PITCH[RpId*4+plane];
      double posMean = ZBCluster->position();
      if(mDoAfterburner)
        posMean -= correction[RpId][coordinate];
      double E = ZBCluster->energy();
      double posMin = posMean - pitch*len/2;
      double posMax = posMean + pitch*len/2;
      bool overlap;
      double pos0, eCl;
      int len0;
      double posL, posH;
      for(unsigned int iCl = 0; iCl < clusterPerPlane[RpId][plane].size(); ++iCl)
      { // loop over clusters in this plane
        overlap = false;
        StUPCRpsCluster *clstr = tmpEvt->getCluster(clusterPerPlane[RpId][plane][iCl]);
        len0 = clstr->length();
        eCl = clstr->energy();
        pos0 = clstr->position();
        if(len0 <= 0 || eCl <= 0)
          continue;
        posL = pos0 - pitch*len0/2;
        posH = pos0 + pitch*len0/2;
        if( posL>=posMin && posH<=posMax ){ // mc cluster fully contained in new cluster 
          overlap = true;  
        }else if( posH>posMax && posL< posMin ){// new cluster fully contained in MC cluster
          overlap = true;
          posMin = posL;
          posMax = posH;
        }else if( posL < posMin && (posH>=posMin && posH<=posMax) ){ // mc/new clusters left overlap
          overlap = true;
          posMin = posL;
        }else if((posL >= posMin && posL<= posMax) && posH>posMax ){ // mc/new clusters right overlap
          overlap = true;
          posMax = posH;
        }
        if(overlap)
        {
          //cout<<"ZB cluster ("<<posMean - pitch*len/2<<", "<<posMean + pitch*len/2<<") overlaps with ("<<posL<<", "<<posH<<")"<<endl;
          posMean = (posMean*E + pos0*eCl)/(E+eCl);
          E += eCl;
          clstr->setLength(0); 
          clstr->setEnergy(0);
        }
      }
      int newLen = posMax-posMin+1;
      StUPCRpsCluster *cluster = tmpEvt->addCluster();
      cluster->setPosition(posMean);
      //clstr->setPositionRMS(); //true RMS calculated from strips in St_pp2pp_Maker/St_pp2pp_Maker.cxx line 841 
      cluster->setLength(newLen); 
      cluster->setEnergy(E); 
      cluster->setXY(ZBCluster->xy()); //clstr->setXY(clstr->xy() + orientations2[4*RpId+plane]*(posMean - clstr->position()));
      cluster->setQuality(ZBCluster->quality()); 
      cluster->setPlaneId(plane);
      cluster->setRomanPotId(RpId);
      clusterPerPlane[RpId][plane].push_back(clusterID++);

    }
  }

  clusterID = rpEvent->getNumberOfClusters();
  for (unsigned int iClstr = 0; iClstr < tmpEvt->getNumberOfClusters(); ++iClstr)
  {
    StUPCRpsCluster *clstr = tmpEvt->getCluster(iClstr);
    if(clstr->length() <= 0 || clstr->energy() <= 0)
      continue;

    StUPCRpsCluster *cluster = rpEvent->addCluster();
    cluster->setPosition(clstr->position());
    cluster->setLength(clstr->length()); 
    cluster->setEnergy(clstr->energy()); 
    cluster->setXY(clstr->xy()); 
    cluster->setPlaneId(clstr->planeId());
    cluster->setRomanPotId(RpId);
    cluster->setQuality(clstr->quality());
    clusterPerPlaneUsed[RpId][clstr->planeId()].push_back(clusterID++);
  }

  delete tmpEvt; 
  tmpEvt = 0;

  for(int j=0; j<kPlanesPerCoordinate; ++j) // j = 0, 1 (kPlanesPerCoordinate)
  { // loop over planes measuring given _coordinate_
    int planeID = kPlanes[coordinate][j];
    int nClusters = clusterPerPlaneUsed[RpId][planeID].size();

    if(nClusters >= kMaxNumberOfClusterPerPlane) // continue only if nClusters is small, otherwise plane is not used in reconstruction
      continue;

    for(int k=0; k < nClusters; ++k)
    { // loop over clusters in this plane
      StUPCRpsCluster *cluster = rpEvent->getCluster(clusterPerPlaneUsed[RpId][planeID][k]);
      int lenCluster = cluster->length();
      double enCluster = cluster->energy();
      eneClusters->Fill(enCluster);
      lenClusters->Fill(lenCluster);
      posClusters->Fill(cluster->xy(), 4*RpId + planeID);
      if(lenCluster > kMaxClusterLength || lenCluster <= 0)
        continue;

      if(enCluster < kEmin[ RpId ][lenCluster-1]) // allow using this cluster only if it passes the energy cut
        continue;

      pos[j].push_back( cluster->xy() );
      en[j].push_back( enCluster );
      len[j].push_back( lenCluster );
      id[j].push_back( clusterPerPlaneUsed[RpId][planeID][k]);
    }
  }
}


Int_t EmbedMaker::classifyClustersCase(vector<double>* pos) const
{

  int lA = pos[kFirst].size();
  int lB = pos[kSecond].size();

  if( lA==0 && lB==0 )  return -1; else
  if( lA==1 && lB==1 )  return 5; else
  if( lA==1 && lB >1 )  return 6; else
  if( lA >1 && lB==1 )  return 7; else
  if( lA==0 && lB==1 )  return 2; else
  if( lA==1 && lB==0 )  return 1; else
  if( lA>=2 && lB>=2 )  return 8; else
  if( lA==0 && lB >1 )  return 4; else
  if( lA >1 && lB==0 )  return 3; else
  return -100;
}


Bool_t EmbedMaker::matchClusters(const int coordinate, const int clCase, 
  const vector<double>* pos, std::vector<unsigned int>* validClusters) const
{
  switch(clCase)
  {
    case 5:
    { 
      if( areMatched(coordinate, pos[kFirst][0], pos[kSecond][0]) )
      {
        validClusters[kFirst].push_back( 0 );
        validClusters[kSecond].push_back( 0 );
        return true;
      } 
      return false;
    }  
    case 6:
    case 7: 
    {
      double DeltaPosition = 1e9;
      double minDeltaPosition = 1e9;
      int index[kPlanesPerCoordinate] = { -1, -1 };
      for(unsigned int c1=0; c1<pos[kFirst].size(); ++c1)
      {
        for(unsigned int c2=0; c2<pos[kSecond].size(); ++c2)
        {
          if(areMatched(coordinate, pos[kFirst][c1], pos[kSecond][c2], &DeltaPosition))
          {
            if(abs(DeltaPosition) < minDeltaPosition)
            {
              minDeltaPosition = abs(DeltaPosition);
              index[kFirst] = c1;
              index[kSecond] = c2;
            }
          }
        }
      }

      if(index[kFirst]>-1)
      {
        validClusters[kFirst].push_back( index[kFirst] );
        validClusters[kSecond].push_back( index[kSecond] );
        return true;
      } else 
        return false;
    }
    case 8: 
    {
      for(unsigned int c1=0; c1<pos[kFirst].size(); ++c1)
      {
        for(unsigned int c2=0; c2<pos[kSecond].size(); ++c2)
        {
          if(areMatched(coordinate, pos[kFirst][c1], pos[kSecond][c2]))
          {
            validClusters[kFirst].push_back( c1 );
            validClusters[kSecond].push_back( c2 );
          }
        }
      }
      if(validClusters[kFirst].size()>0) 
        return true;
      else 
        return false;
    }
    default: 
      return false;
  }
  return false;
}


Bool_t EmbedMaker::areMatched(const int coordinate, const double p1, const double p2, 
  double *deltaPitches) const
{
  if(deltaPitches) 
    *deltaPitches = (p1 - p2) / kPitch[coordinate];
  return  abs( p1 - p2 ) < kMaxPitchesToMatch*kPitch[coordinate] ? true : false;
}

Double_t EmbedMaker::timeFromTAC(const int Rp, const int pmt, const int Tac, const int Adc) const{
  return mConversion_TAC_time * (Tac + mSkew_param[Rp][pmt][0]
  + mSkew_param[Rp][pmt][1]*exp(-mSkew_param[Rp][pmt][2]*(Adc-mSkew_param[Rp][pmt][3])) / (Adc-mSkew_param[Rp][pmt][3]) );
}

void EmbedMaker::Init(TFile *outfile)
{
   if( DEBUG )
      cout<<"EmbedMaker::Init() called"<<endl;
  outfile->mkdir("EmbedMaker")->cd();
  eneSimHits = new TH2D("eneSimHits"," Energy of hits ; ADC; plane",255,-0.5,254.5,32,0.,32.);
  numSimHits = new TH2D("numSimHits"," Number of hits ; N_{hits}; plane",255,-0.5,254.5,32,0.,32.);
  posSimHits = new TH2D("posSimHits"," Position of hits ;channel ; plane",100,-0.05,0.05,32,0.,32.);
  hmatchDist  = new TH2D("hmatchDist","Cluster matching distance ; dist [mm]; plane",40,-1.,1., 16, 0., 16. );
  numClusters = new TH2D("numClusters"," Number of Clusters per plane; N_{clu}; Entries",50,0.,49.5, 32, 0., 32. );
  numClustersUsed = new TH2D("numClustersUsed"," Number of used Clusters per plane; N_{clu}; Entries",50,0.,49.5, 32, 0., 32. );
  posClusters = new TH2D("posClusters"," Cluster position; local position [m] ; PlaneName",160,-0.08,0.08, 32, 0., 32.);
  lenClusters = new TH1D("lenClusters"," Cluster length; N_{strips} ; N_{clu}",128,-0.5,127.5 );
  eneClusters = new TH1D("eneClusters"," Cluster energy; ADC ; N_{clu}",500,-0.5,499.5 );
  numPoints   = new TH2D("numPoints","Number of reco 2D points; N_{points}; RP_NAME",50,-0.5,49.5, 8,0.,8.);
  posPoints   = new TH2D("posPoints","Position of reco points; local position [m] ; RP_NAME",320,-0.08 ,0.08, 16,0.,16.);
  

  TString rpNames[kMAXSEQ] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};
  for (int i = 0; i < kMAXSEQ; ++i){
    hGoodTrackPoints[i] = new TH1D("hGoodTrackPoints_"+rpNames[i],"Number of trackpoints in "+rpNames[i], 30, -0.5, 29.5);
    hHitSizeMap[i] = new TH2D("hHitSizeMap_"+rpNames[i],"Number of hits in Y vs X for "+rpNames[i], 20, -0.5, 19.5, 20, -0.5, 19.5);
  }


  outfile->cd();
}

void EmbedMaker::setOffsetCorrections(TVector3 corr[kMAXSEQ]){
  for (int i = 0; i < kMAXSEQ; ++i)
    correction[i] =  corr[i];
}

const double EmbedMaker::kEmin[kMAXSEQ][kMaxClusterLength] =
{{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20},
{20, 20, 20, 20, 20}};

const unsigned int EmbedMaker::kPlanes[kCoordinates][kPlanesPerCoordinate] =
{{1, 3},  // kX (vertical strips)
{0, 2}}; // kY (horizontal strips)

const int EmbedMaker::kRpInBranch[kBranches][kStationsPerBranch] =
{{0, 2}, {1, 3},
{4, 6}, {5, 7}};

const double EmbedMaker::kPitch[kCoordinates]=
{ 1.050E-4, 9.74E-5 };


const double EmbedMaker::STRIP_PITCH[kMAXSEQ*kMAXCHAIN] = {
  0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03,
  0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0955E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03,
  0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03,
  0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03, 0.0974E-03, 0.1050E-03
};

const double EmbedMaker::RP_POS_Z[ kMAXSEQ*kMAXCHAIN ] = { 
  -15.76908,   -15.77578,   -15.78248,   -15.78918,
  -15.77078,   -15.77748,   -15.78418,   -15.79088,
  -17.56906,   -17.57576,   -17.58246,   -17.58916,
  -17.56906,   -17.57576,   -17.58246,   -17.58916,
  15.76912,    15.77582,    15.78252,    15.78922,
  15.76977,    15.77647,    15.78317,    15.78987,
  17.57023,    17.57693,    17.58363,    17.59033,
  17.56916,    17.57586,    17.58256,    17.58926
};

const short EmbedMaker::orientations2[kMAXCHAIN*kMAXSEQ] = {1,1,1,1,  -1,-1,-1,-1,  1,1,1,1, -1,-1,-1,-1,  1,-1,1,-1,  -1,1,-1,1,  1,-1,1,-1, -1,1,-1,1 };
