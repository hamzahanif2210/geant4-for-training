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
/// \file B4/B4a/src/SteppingAction.cc
/// \brief Implementation of the B4a::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
using namespace B4;

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // Step coordinates
  auto x = step->GetPostStepPoint()->GetPosition().x();
  auto y = step->GetPostStepPoint()->GetPosition().y();
  auto z = step->GetPostStepPoint()->GetPosition().z();
  auto timeNs = step->GetPreStepPoint()->GetGlobalTime() / ns;

  auto touchable = step->GetPreStepPoint()->GetTouchableHandle();
  auto copyNumber = touchable->GetCopyNumber();

  int nCellsXY = 5;
  int layerIndex = copyNumber / (nCellsXY * nCellsXY);


  // get event number
  auto eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // get gun energy
  auto gunEnergy = G4RunManager::GetRunManager()->GetCurrentEvent()->GetPrimaryVertex(0)->GetPrimary()->GetKineticEnergy();

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill ntuple with event ID too
  analysisManager->FillNtupleIColumn(0, eventID);
  analysisManager->FillNtupleDColumn(1, gunEnergy);
  analysisManager->FillNtupleDColumn(2, x);
  analysisManager->FillNtupleDColumn(3, y);
  analysisManager->FillNtupleDColumn(4, z);
  analysisManager->FillNtupleDColumn(5, edep);
  analysisManager->FillNtupleDColumn(6, timeNs);
  analysisManager->FillNtupleIColumn(7, layerIndex);
  analysisManager->FillNtupleSColumn(8, std::string(fDetConstruction->GetMaterialType()));
  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
