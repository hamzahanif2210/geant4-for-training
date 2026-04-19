//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file exampleB4a.cc
/// \brief Main program of the B4a example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4a [-m macro] [-u UIsession] [-t nThreads] [-s seed]"
           << G4endl;
    G4cerr << "            [-cx cellSizeX_cm] [-cy cellSizeY_cm] [-cz cellSizeZ_cm]"
           << G4endl;
    G4cerr << "            [-mat materialType] [-emin minEnergyGeV] [-emax maxEnergyGeV]"
           << G4endl;
    G4cerr << "            [-o outputFile] [-vDefault]"
           << G4endl;
    G4cerr << "   Cell sizes are in cm. Default: 4 4 10 (= 4x4x10 cm^3 per cell)"
           << G4endl;
    G4cerr << "   Material: PbF2 (default) or PbWO4"
           << G4endl;
    G4cerr << "   Energy range in GeV. Default: emin=1 emax=5 (uniform random per event)"
           << G4endl;
    G4cerr << "   -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  G4String macro;
  G4String session;
  G4int seed = 12345;
  G4bool verboseBestUnits = true;
  G4double cellSizeX = 4.0;  // cm
  G4double cellSizeY = 4.0;  // cm
  G4double cellSizeZ = 10.0; // cm
  G4String materialType = "PbF2";
  G4double energyMin = 1.0;  // GeV
  G4double energyMax = 5.0;  // GeV
  G4String outputFile = "";
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m"   ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-s"   ) seed = std::stoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-u"   ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-cx"  ) cellSizeX = std::stod(argv[i+1]);
    else if ( G4String(argv[i]) == "-cy"  ) cellSizeY = std::stod(argv[i+1]);
    else if ( G4String(argv[i]) == "-cz"  ) cellSizeZ = std::stod(argv[i+1]);
    else if ( G4String(argv[i]) == "-mat" ) materialType = argv[i+1];
    else if ( G4String(argv[i]) == "-emin") energyMin = std::stod(argv[i+1]);
    else if ( G4String(argv[i]) == "-emax") energyMax = std::stod(argv[i+1]);
    else if ( G4String(argv[i]) == "-o"   ) outputFile = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t"  ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else if ( G4String(argv[i]) == "-vDefault" ) {
      verboseBestUnits = false;
      --i;  // this option is not followed with a parameter
    }
    else {
      PrintUsage();
      return 1;
    }
  }

  // Auto-build output filename if not explicitly provided
  if ( outputFile.empty() ) {
    auto fmtNum = [](G4double v) -> G4String {
      G4int iv = (G4int)v;
      return (G4double)iv == v ? std::to_string(iv) : std::to_string(v);
    };
    outputFile = "photons_" + fmtNum(cellSizeX) + "x" + fmtNum(cellSizeY)
               + "x" + fmtNum(cellSizeZ) + "cm_"
               + fmtNum(energyMin) + "to" + fmtNum(energyMax) + "GeV_"
               + materialType + ".root";
  }

  // Random seed
  CLHEP::HepRandom::setTheSeed(seed);
  G4Random::setTheSeed(seed);

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  auto RNG = new CLHEP::MTwistEngine(seed);
  G4Random::setTheEngine(RNG);

  // Use G4SteppingVerboseWithUnits
  if ( verboseBestUnits ) {
    G4int precision = 4;
    G4SteppingVerbose::UseBestUnit(precision);
  }

  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
#endif

  // Set mandatory initialization classes
  //
  auto detConstruction = new B4::DetectorConstruction(cellSizeX, cellSizeY, cellSizeZ, materialType);
  runManager->SetUserInitialization(detConstruction);

  auto physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  auto actionInitialization = new B4a::ActionInitialization(detConstruction, outputFile, energyMin, energyMax);
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  //
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
