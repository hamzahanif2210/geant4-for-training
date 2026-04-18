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

#include <cstdlib>
#include <fstream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4a [-m macro ] [-u UIsession] [-t nThreads] [-s seed] [-n nEvents]"
           << G4endl;
    G4cerr << "           [-e energyGeV] [-c cellSizeCm] [--cell-x cm] [--cell-y cm] [--cell-z cm]"
           << G4endl;
    G4cerr << "           [-o outputFile] [-vDefault]"
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }

  G4bool MacroContainsBeamOn(const G4String& macroPath) {
    std::ifstream macroFile(macroPath);
    if (!macroFile) {
      return false;
    }

    std::string line;
    while (std::getline(macroFile, line)) {
      const auto firstNonSpace = line.find_first_not_of(" 	");
      if (firstNonSpace == std::string::npos) {
        continue;
      }
      if (line[firstNonSpace] == '#') {
        continue;
      }
      if (line.compare(firstNonSpace, 11, "/run/beamOn") == 0) {
        return true;
      }
    }

    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  G4String macro;
  G4String session;
  G4int seed = 12345;
  G4int nEvents = -1;
  G4double primaryEnergyGeV = -1.;
  G4double isotropicCellSizeCm = -1.;
  G4double cellSizeXCm = -1.;
  G4double cellSizeYCm = -1.;
  G4double cellSizeZCm = -1.;
  G4String outputFile;
  G4bool verboseBestUnits = true;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i = 1; i < argc; ++i ) {
    G4String option = argv[i];
    auto requireValue = [&](const G4String& optName) -> char* {
      if (i + 1 >= argc) {
        G4cerr << "Missing value for " << optName << G4endl;
        PrintUsage();
        std::exit(1);
      }
      ++i;
      return argv[i];
    };

    if      ( option == "-m" ) macro = requireValue(option);
    else if ( option == "-s" ) seed = std::stoi(requireValue(option));
    else if ( option == "-u" ) session = requireValue(option);
    else if ( option == "-n" ) nEvents = std::stoi(requireValue(option));
    else if ( option == "-e" ) primaryEnergyGeV = std::stod(requireValue(option));
    else if ( option == "-c" ) isotropicCellSizeCm = std::stod(requireValue(option));
    else if ( option == "--cell-x" ) cellSizeXCm = std::stod(requireValue(option));
    else if ( option == "--cell-y" ) cellSizeYCm = std::stod(requireValue(option));
    else if ( option == "--cell-z" ) cellSizeZCm = std::stod(requireValue(option));
    else if ( option == "-o" ) outputFile = requireValue(option);
#ifdef G4MULTITHREADED
    else if ( option == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(requireValue(option));
    }
#endif
    else if ( option == "-vDefault" ) {
      verboseBestUnits = false;
    }
    else {
      PrintUsage();
      return 1;
    }
  }

  // Random seed
  CLHEP::HepRandom::setTheSeed(seed); 
  G4Random::setTheSeed(seed);

  if (primaryEnergyGeV > 0.) {
    setenv("B4_PRIMARY_ENERGY_GEV", std::to_string(primaryEnergyGeV).c_str(), 1);
  }
  if (isotropicCellSizeCm > 0.) {
    const auto cellSizeValue = std::to_string(isotropicCellSizeCm);
    setenv("B4_CELL_SIZE_CM", cellSizeValue.c_str(), 1);
  }
  if (cellSizeXCm > 0.) {
    setenv("B4_CELL_SIZE_X_CM", std::to_string(cellSizeXCm).c_str(), 1);
  }
  if (cellSizeYCm > 0.) {
    setenv("B4_CELL_SIZE_Y_CM", std::to_string(cellSizeYCm).c_str(), 1);
  }
  if (cellSizeZCm > 0.) {
    setenv("B4_CELL_SIZE_Z_CM", std::to_string(cellSizeZCm).c_str(), 1);
  }
  if (!outputFile.empty()) {
    setenv("B4_OUTPUT_FILE", outputFile.c_str(), 1);
  }

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
  auto detConstruction = new B4::DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  auto physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  auto actionInitialization = new B4a::ActionInitialization(detConstruction);
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
    const auto macroHasBeamOn = MacroContainsBeamOn(macro);
    UImanager->ApplyCommand(command+macro);
    if (nEvents > 0) {
      if (macroHasBeamOn) {
        G4cout << "[exampleB4a] '-n " << nEvents
               << "' ignored because macro already contains /run/beamOn."
               << G4endl;
      }
      else {
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(nEvents));
      }
    }
  }
  else  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    if (nEvents > 0) {
      UImanager->ApplyCommand("/run/initialize");
      UImanager->ApplyCommand("/run/beamOn " + std::to_string(nEvents));
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
