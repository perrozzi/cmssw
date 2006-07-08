/*
 * \file EBPedestalTask.cc
 *
 * $Date: 2006/07/08 07:22:02 $
 * $Revision: 1.42 $
 * \author G. Della Ricca
 *
*/

#include <DQM/EcalBarrelMonitorTasks/interface/EBPedestalTask.h>

EBPedestalTask::EBPedestalTask(const ParameterSet& ps){

  init_ = false;

  for (int i = 0; i < 36 ; i++) {
    mePedMapG01_[i] = 0;
    mePedMapG06_[i] = 0;
    mePedMapG12_[i] = 0;
    mePed3SumMapG01_[i] = 0;
    mePed3SumMapG06_[i] = 0;
    mePed3SumMapG12_[i] = 0;
    mePed5SumMapG01_[i] = 0;
    mePed5SumMapG06_[i] = 0;
    mePed5SumMapG12_[i] = 0;
    mePnPedMapG01_[i] = 0;
    mePnPedMapG16_[i] = 0;
  }

}

EBPedestalTask::~EBPedestalTask(){

}

void EBPedestalTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  DaqMonitorBEInterface* dbe = 0;

  // get hold of back-end interface
  dbe = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe ) {
    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask");
    dbe->rmdir("EcalBarrel/EBPedestalTask");
  }

}

void EBPedestalTask::setup(void){

  init_ = true;

  Char_t histo[200];

  DaqMonitorBEInterface* dbe = 0;

  // get hold of back-end interface
  dbe = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe ) {
    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask");

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain01");
    for (int i = 0; i < 36 ; i++) {
      sprintf(histo, "EBPT pedestal SM%02d G01", i+1);
      mePedMapG01_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 3sum SM%02d G01", i+1);
      mePed3SumMapG01_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 5sum SM%02d G01", i+1);
      mePed5SumMapG01_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
    }

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain06");
    for (int i = 0; i < 36 ; i++) {
      sprintf(histo, "EBPT pedestal SM%02d G06", i+1);
      mePedMapG06_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 3sum SM%02d G06", i+1);
      mePed3SumMapG06_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 5sum SM%02d G06", i+1);
      mePed5SumMapG06_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
    }

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain12");
    for (int i = 0; i < 36 ; i++) {
      sprintf(histo, "EBPT pedestal SM%02d G12", i+1);
      mePedMapG12_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 3sum SM%02d G12", i+1);
      mePed3SumMapG12_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
      sprintf(histo, "EBPT pedestal 5sum SM%02d G12", i+1);
      mePed5SumMapG12_[i] = dbe->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 4096, 0., 4096., "s");
    }

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask");

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask/Gain01");
    for (int i = 0; i < 36 ; i++) {
      sprintf(histo, "EBPDT PNs pedestal SM%02d G01", i+1);
      mePnPedMapG01_[i] =  dbe->bookProfile2D(histo, histo, 1, 0., 1., 10, 0., 10., 4096, 0., 4096., "s");
    }

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask/Gain16");
    for (int i = 0; i < 36 ; i++) {
      sprintf(histo, "EBPDT PNs pedestal SM%02d G16", i+1);
      mePnPedMapG16_[i] =  dbe->bookProfile2D(histo, histo, 1, 0., 1., 10, 0., 10., 4096, 0., 4096., "s");
    }

  }

}

void EBPedestalTask::cleanup(void){

  DaqMonitorBEInterface* dbe = 0;

  // get hold of back-end interface
  dbe = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe ) {

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain01");
    for ( int i = 0; i < 36; i++ ) {
      if ( mePedMapG01_[i] ) dbe->removeElement( mePedMapG01_[i]->getName() );
      mePedMapG01_[i] = 0;
      if ( mePed3SumMapG01_[i] ) dbe->removeElement( mePed3SumMapG01_[i]->getName() );
      mePed3SumMapG01_[i] = 0;
      if ( mePed5SumMapG01_[i] ) dbe->removeElement( mePed5SumMapG01_[i]->getName() );
      mePed5SumMapG01_[i] = 0;
    }

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain06");
    for ( int i = 0; i < 36; i++ ) {
      if ( mePedMapG06_[i] ) dbe->removeElement( mePedMapG06_[i]->getName() );
      mePedMapG06_[i] = 0;
      if ( mePed3SumMapG06_[i] ) dbe->removeElement( mePed3SumMapG06_[i]->getName() );
      mePed3SumMapG06_[i] = 0;
      if ( mePed5SumMapG06_[i] ) dbe->removeElement( mePed5SumMapG06_[i]->getName() );
      mePed5SumMapG06_[i] = 0;
    }

    dbe->setCurrentFolder("EcalBarrel/EBPedestalTask/Gain12");
    for ( int i = 0; i < 36; i++ ) {
      if ( mePedMapG12_[i] ) dbe->removeElement( mePedMapG12_[i]->getName() );
      mePedMapG12_[i] = 0;
      if ( mePed3SumMapG12_[i] ) dbe->removeElement( mePed3SumMapG12_[i]->getName() );
      mePed3SumMapG12_[i] = 0;
      if ( mePed5SumMapG12_[i] ) dbe->removeElement( mePed5SumMapG12_[i]->getName() );
      mePed5SumMapG12_[i] = 0;
    }

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask");

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask/Gain01");
    for ( int i = 0; i < 36; i++ ) {
      if ( mePnPedMapG01_[i]) dbe->removeElement( mePnPedMapG01_[i]->getName() );
      mePnPedMapG01_[i] = 0;
    }

    dbe->setCurrentFolder("EcalBarrel/EBPnDiodeTask/Gain16");
    for ( int i = 0; i < 36; i++ ) {
      if ( mePnPedMapG16_[i]) dbe->removeElement( mePnPedMapG16_[i]->getName() );
      mePnPedMapG16_[i] = 0;
    }

  }

  init_ = false;

}

void EBPedestalTask::endJob(void){

  LogInfo("EBPedestalTask") << "analyzed " << ievt_ << " events";

  this->cleanup();

}

void EBPedestalTask::analyze(const Event& e, const EventSetup& c){

  bool enable = false;
  map<int, EcalDCCHeaderBlock> dccMap;

  Handle<EcalRawDataCollection> dcchs;
  e.getByLabel("ecalEBunpacker", dcchs);

  for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin(); dcchItr != dcchs->end(); ++dcchItr ) {

    EcalDCCHeaderBlock dcch = (*dcchItr);

    dccMap[dcch.id()] = dcch;

    if ( dccMap[dcch.id()].getRunType() == EcalDCCHeaderBlock::PEDESTAL_STD ) enable = true;

  }

  if ( ! enable ) return;

  if ( ! init_ ) this->setup();

  ievt_++;

  Handle<EBDigiCollection> digis;
  e.getByLabel("ecalEBunpacker", digis);

  int nebd = digis->size();
  LogDebug("EBPedestalTask") << "event " << ievt_ << " digi collection size " << nebd;

  float xmap01[36][85][20];
  float xmap06[36][85][20];
  float xmap12[36][85][20];

  for ( int ism = 1; ism <= 36; ism++ ) {
    for ( int ie = 1; ie <= 85; ie++ ) {
      for ( int ip = 1; ip <= 20; ip++ ) {

	 xmap01[ism-1][ie-1][ip-1] = 0.;
         xmap06[ism-1][ie-1][ip-1] = 0.;
         xmap12[ism-1][ie-1][ip-1] = 0.;

      }
    }
  }

  for ( EBDigiCollection::const_iterator digiItr = digis->begin(); digiItr != digis->end(); ++digiItr ) {

    EBDataFrame dataframe = (*digiItr);
    EBDetId id = dataframe.id();

    int ic = id.ic();
    int ie = (ic-1)/20 + 1;
    int ip = (ic-1)%20 + 1;

    int ism = id.ism();

    float xie = ie - 0.5;
    float xip = ip - 0.5;

    if ( dccMap[ism-1].getRunType() != EcalDCCHeaderBlock::PEDESTAL_STD ) continue;

    LogDebug("EBPedestalTask") << " det id = " << id;
    LogDebug("EBPedestalTask") << " sm, eta, phi " << ism << " " << ie << " " << ip;

    for (int i = 0; i < 10; i++) {

      EcalMGPASample sample = dataframe.sample(i);
      int adc = sample.adc();

      MonitorElement* mePedMap = 0;

      if ( sample.gainId() == 1 ) mePedMap = mePedMapG12_[ism-1];
      if ( sample.gainId() == 2 ) mePedMap = mePedMapG06_[ism-1];
      if ( sample.gainId() == 3 ) mePedMap = mePedMapG01_[ism-1];

      float xval = float(adc);

      if ( mePedMap ) mePedMap->Fill(xie, xip, xval);

      if ( sample.gainId() == 1 ) xmap12[ism-1][ie-1][ip-1] = xval;
      if ( sample.gainId() == 2 ) xmap06[ism-1][ie-1][ip-1] = xval;
      if ( sample.gainId() == 3 ) xmap01[ism-1][ie-1][ip-1] = xval;

    }

  }

  // to be re-done using the 3x3 & 5x5 Selectors (if faster)

  for ( int ism = 1; ism <= 36; ism++ ) {
    for ( int ie = 1; ie <= 85; ie++ ) {
      for ( int ip = 1; ip <= 20; ip++ ) {

        float xie = ie - 0.5;
        float xip = ip - 0.5;

        float x3val01;
        float x3val06;
        float x3val12;

        if ( ie >= 2 && ie <= 84 && ip >= 2 && ip <= 19 ) {

          x3val01 = 0.;
          x3val06 = 0.;
          x3val12 = 0.;
          for ( int i = -1; i <= +1; i++ ) {
            for ( int j = -1; j <= +1; j++ ) {

              x3val01 = x3val01 + xmap01[ism-1][ie-1+i][ip-1+j];
              x3val06 = x3val06 + xmap06[ism-1][ie-1+i][ip-1+j];
              x3val12 = x3val12 + xmap12[ism-1][ie-1+i][ip-1+j];

            }
          }
          x3val01 = x3val01 / 9.;
          x3val06 = x3val06 / 9.;
          x3val12 = x3val12 / 9.;
          if ( mePed3SumMapG01_[ism-1] && x3val01 != 0. ) mePed3SumMapG01_[ism-1]->Fill(xie, xip, x3val01);
          if ( mePed3SumMapG06_[ism-1] && x3val06 != 0. ) mePed3SumMapG06_[ism-1]->Fill(xie, xip, x3val06);
          if ( mePed3SumMapG12_[ism-1] && x3val12 != 0. ) mePed3SumMapG12_[ism-1]->Fill(xie, xip, x3val12);

        }

        float x5val01;
        float x5val06;
        float x5val12;

        if ( ie >= 3 && ie <= 83 && ip >= 3 && ip <= 18 ) {

          x5val01 = 0.;
          x5val06 = 0.;
          x5val12 = 0.;
          for ( int i = -2; i <= +2; i++ ) {
            for ( int j = -2; j <= +2; j++ ) {

              x5val01 = x5val01 + xmap01[ism-1][ie-1+i][ip-1+j];
              x5val06 = x5val06 + xmap06[ism-1][ie-1+i][ip-1+j];
              x5val12 = x5val12 + xmap12[ism-1][ie-1+i][ip-1+j];

            }
          }
          x5val01 = x5val01 / 25.;
          x5val06 = x5val06 / 25.;
          x5val12 = x5val12 / 25.;
          if ( mePed5SumMapG01_[ism-1] && x5val01 != 0. ) mePed5SumMapG01_[ism-1]->Fill(xie, xip, x5val01);
          if ( mePed5SumMapG06_[ism-1] && x5val06 != 0. ) mePed5SumMapG06_[ism-1]->Fill(xie, xip, x5val06);
          if ( mePed5SumMapG12_[ism-1] && x5val12 != 0. ) mePed5SumMapG12_[ism-1]->Fill(xie, xip, x5val12);

        }

      }
    }
  }

  Handle<EcalPnDiodeDigiCollection> pns;
  e.getByLabel("ecalEBunpacker", pns);

  int nep = pns->size();
  LogDebug("EBPedestalTask") << "event " << ievt_ << " pns collection size " << nep;

  for ( EcalPnDiodeDigiCollection::const_iterator pnItr = pns->begin(); pnItr != pns->end(); ++pnItr ) {

    EcalPnDiodeDigi pn = (*pnItr);
    EcalPnDiodeDetId id = pn.id();

//    int ism = id.ism();
    int ism = id.iDCCId();

    int num = id.iPnId();

    if ( dccMap[ism-1].getRunType() != EcalDCCHeaderBlock::PEDESTAL_STD ) continue;

    LogDebug("EBPedestalTask") << " det id = " << id;
    LogDebug("EBPedestalTask") << " sm, num " << ism << " " << num;

    for (int i = 0; i < 50; i++) {

      EcalFEMSample sample = pn.sample(i);
      int adc = sample.adc();

      MonitorElement* mePNPed = 0;

      if ( sample.gainId() == 0 ) mePNPed = mePnPedMapG01_[ism-1];
      if ( sample.gainId() == 1 ) mePNPed = mePnPedMapG16_[ism-1];

      float xval = float(adc);

      if ( mePNPed ) mePNPed->Fill(0.5, num - 0.5, xval);

    }

  }

}

