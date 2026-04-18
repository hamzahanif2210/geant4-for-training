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
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied,  *
// * regarding  this  software system or assume any liability for its  *
// * use.  Please see the license in the file  LICENSE  and URL above  *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its   *
// * use  in  resulting  scientific  publications,  and indicate your   *
// * acceptance of all terms of the Geant4 Software license.           *
// ********************************************************************
//
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <vector>
#include <random>
#include <cstdlib>
#include <stdexcept>

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  auto getEnvDouble = [](const char* name, G4double fallback) {
    const char* value = std::getenv(name);
    if (!value) {
      return fallback;
    }

    try {
      return std::stod(value);
    }
    catch (const std::exception&) {
      G4cout << "[B4] Invalid value for " << name << ": " << value
             << ". Using fallback " << fallback << G4endl;
      return fallback;
    }
  };

  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  // Set energy from fixed discrete values {1, 10, 50, 100, 200 GeV}
  // std::vector<G4double> energies = {1.*GeV, 10.*GeV, 20.*GeV, 30.*GeV, 40.*GeV, 50.*GeV, 60.*GeV, 70.*GeV, 80.*GeV, 90.*GeV, 100.*GeV};
  // std::random_device rd;  
  // std::mt19937 gen(rd()); // Mersenne Twister RNG seeded by rd()
  // std::uniform_int_distribution<size_t> dist(0, energies.size() - 1);
  // auto particleEnergy = energies[dist(gen)];
  // fParticleGun->SetParticleEnergy(particleEnergy);
  // fParticleGun->SetParticleEnergy(1000.*MeV); // 1000 MeV = 1 GeV
  const G4double primaryEnergyGeV = getEnvDouble("B4_PRIMARY_ENERGY_GEV", 15.0);
  fParticleGun->SetParticleEnergy(primaryEnergyGeV * GeV);
  G4cout << "[B4] Primary energy: " << primaryEnergyGeV << " GeV" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the beginning of each event

  // Get the world volume to determine the gun's starting position
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has a box shape
  G4Box* worldBox = nullptr;
  if (worldLV) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if (worldBox) {
    worldZHalfLength = worldBox->GetZHalfLength();
  } else {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be placed in the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
                "MyCode0002", JustWarning, msg);
  }

  // Set gun position
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -(worldZHalfLength + 0.*mm)));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

} // namespace B4
