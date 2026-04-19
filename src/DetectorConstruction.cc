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

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4double cellSizeX, G4double cellSizeY, G4double cellSizeZ,
                                           const G4String& materialType)
  : fCellSizeX(cellSizeX), fCellSizeY(cellSizeY), fCellSizeZ(cellSizeZ), fMaterialType(materialType)
{}

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
  // Geometry parameters — cell sizes provided in cm, converted to Geant4 internal units
  G4double cellSizeX = fCellSizeX * cm;
  G4double cellSizeY = fCellSizeY * cm;
  G4double cellSizeZ = fCellSizeZ * cm;
  G4int nCellsXY = 5;
  G4int nCellsZ  = 5;
  G4double calorSizeX = nCellsXY * cellSizeX;
  G4double calorSizeY = nCellsXY * cellSizeY;
  G4double calorSizeZ = nCellsZ  * cellSizeZ;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  G4String g4MatName = (fMaterialType == "PbWO4") ? "G4_PbWO4" : "LeadFluoride";
  auto detectorMaterial = G4Material::GetMaterial(g4MatName);

  if ( ! defaultMaterial || ! detectorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS = new G4Box("World", calorSizeX/2, calorSizeY/2, calorSizeZ/2);

  auto worldLV = new G4LogicalVolume(worldS, defaultMaterial, "World");

  auto worldPV = new G4PVPlacement(
    0,
    G4ThreeVector(0., 0., 0.),
    worldLV,
    "World",
    0,
    false,
    0,
    false
  );

  auto solidDetector = new G4Box("solidDetector", cellSizeX/2, cellSizeY/2, cellSizeZ/2);
  auto logicDetector = new G4LogicalVolume(solidDetector, detectorMaterial, "logicDetector");
  G4int copyNum = 0;
  for (G4int k = 0; k < nCellsZ; k++) {
    for (G4int j = 0; j < nCellsXY; j++) {
      for (G4int i = 0; i < nCellsXY; i++) {
        new G4PVPlacement(
          0,
          G4ThreeVector(
            -(calorSizeX - cellSizeX)/2. + cellSizeX*i,
            -(calorSizeY - cellSizeY)/2. + cellSizeY*j,
            -(calorSizeZ - cellSizeZ)/2. + cellSizeZ*k
          ),
          logicDetector,
          "physDetector",
          worldLV,
          false,
          copyNum++,
          false
        );
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
