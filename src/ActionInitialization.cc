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
/// \file B4/B4a/src/ActionInitialization.cc
/// \brief Implementation of the B4a::ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

using namespace B4;

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
 : fDetConstruction(detConstruction)
{
  fGunMessenger = new G4GenericMessenger(this, "/B4/gun/", "Incident particle commands");
  fGunMessenger->DeclareMethodWithUnit("emin", "GeV",
    &ActionInitialization::CmdSetEmin, "Set minimum incident energy");
  fGunMessenger->DeclareMethodWithUnit("emax", "GeV",
    &ActionInitialization::CmdSetEmax, "Set maximum incident energy");

  fRunMessenger = new G4GenericMessenger(this, "/B4/run/", "Run control commands");
  fRunMessenger->DeclareMethod("outputFile",
    &ActionInitialization::CmdSetOutputFile, "Set output ROOT file name");
}

ActionInitialization::~ActionInitialization()
{
  delete fGunMessenger;
  delete fRunMessenger;
}

void ActionInitialization::CmdSetEmin(G4double e)
{
  PrimaryGeneratorAction::SetEmin(e);
}

void ActionInitialization::CmdSetEmax(G4double e)
{
  PrimaryGeneratorAction::SetEmax(e);
}

void ActionInitialization::CmdSetOutputFile(G4String name)
{
  RunAction::SetOutputFileName(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction);
  auto eventAction = new EventAction;
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(fDetConstruction,eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
