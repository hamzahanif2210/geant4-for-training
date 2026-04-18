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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_PbWO4");//PbWO4
  nistManager->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");//BaF2
  G4Element* elPb = new G4Element("Lead", "Pb", 82, 207.2*g/mole);
  G4Element* elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole);

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  G4int numberOfatoms, nComp;
  G4Material *PbF2 = new G4Material("LeadFluoride", density = 8.445 * g / cm3, nComp = 2);
  PbF2->AddElement(elPb, numberOfatoms = 1);
  PbF2->AddElement(elF, numberOfatoms = 2);
 
  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
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

  auto getEnvInt = [](const char* name, G4int fallback) {
    const char* value = std::getenv(name);
    if (!value) {
      return fallback;
    }

    try {
      auto parsed = std::stoi(value);
      return parsed > 0 ? parsed : fallback;
    }
    catch (const std::exception&) {
      G4cout << "[B4] Invalid value for " << name << ": " << value
             << ". Using fallback " << fallback << G4endl;
      return fallback;
    }
  };

  // Geometry parameters
  G4double defaultCellSizeCm = getEnvDouble("B4_CELL_SIZE_CM", 4.0);
  G4double cellSizeX = getEnvDouble("B4_CELL_SIZE_X_CM", defaultCellSizeCm) * cm;
  G4double cellSizeY = getEnvDouble("B4_CELL_SIZE_Y_CM", defaultCellSizeCm) * cm;
  G4double cellSizeZ = getEnvDouble("B4_CELL_SIZE_Z_CM", defaultCellSizeCm) * cm;
  G4int nCellsX = getEnvInt("B4_NUM_CELLS_X", 5);
  G4int nCellsY = getEnvInt("B4_NUM_CELLS_Y", 5);
  G4int nCellsZ = getEnvInt("B4_NUM_CELLS_Z", 5);
  G4double calorSizeX = nCellsX * cellSizeX;
  G4double calorSizeY = nCellsY * cellSizeY;
  G4double calorSizeZ  = nCellsZ  * cellSizeZ;
  G4cout << "[B4] Geometry configuration: "
         << nCellsX << "x" << nCellsY << "x" << nCellsZ << " cells, "
         << "cell size "
         << cellSizeX / cm << "x" << cellSizeY / cm << "x" << cellSizeZ / cm
         << " cm^3" << G4endl;

  auto worldSizeX = calorSizeX;
  auto worldSizeY = calorSizeY;
  auto worldSizeZ = calorSizeZ;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto detectorMaterial = G4Material::GetMaterial("LeadFluoride");

  if ( ! defaultMaterial || ! detectorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS = new G4Box("World", worldSizeX/2, worldSizeY/2, worldSizeZ/2);

  auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");

  auto worldPV = new G4PVPlacement(
    0,                                // no rotation 
    G4ThreeVector(0.,0.,0.),          // at (0,0,0)
    worldLV,                          // its logical volume
    "World",                          // its name
    0,                                // its mother  volume
    false,                            // no boolean operation
    0,                                // copy number
    false                             // checking overlaps
  );


  auto solidDetector = new G4Box("solidDetector", cellSizeX/2, cellSizeY/2, cellSizeZ/2);
  auto logicDetector= new G4LogicalVolume(solidDetector,detectorMaterial,"logicDetector");
  G4int copyNumber = 0;
  for (G4int k = 0; k < nCellsZ; k++) {
    for (G4int j = 0; j < nCellsY; j++) {
      for (G4int i = 0; i < nCellsX; i++) {
        new G4PVPlacement(
          0,
          G4ThreeVector(
            -(calorSizeX - cellSizeX) / 2. + (cellSizeX * i),
            -(calorSizeY - cellSizeY) / 2. + (cellSizeY * j),
            -(calorSizeZ - cellSizeZ) / 2. + (cellSizeZ * k)
          ),
          logicDetector,
          "physDetector",
          worldLV,
          false,
          copyNumber,
          false
        );
        ++copyNumber;
      }
    }
  }
  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

} // end namespace
