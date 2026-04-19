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
/// \file B4/B4a/include/DetectorConstruction.hh
/// \brief Definition of the B4::DetectorConstruction class

#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

namespace B4
{

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    // cellSizeX, cellSizeY, cellSizeZ are in cm; materialType is "PbF2" or "PbWO4"
    DetectorConstruction(G4double cellSizeX = 4.0, G4double cellSizeY = 4.0, G4double cellSizeZ = 10.0,
                         const G4String& materialType = "PbF2");
    ~DetectorConstruction() override = default;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // get methods
    //
    const G4VPhysicalVolume* GetAbsorberPV() const;
    const G4VPhysicalVolume* GetGapPV() const;
    const G4String& GetMaterialType() const { return fMaterialType; }
    G4double GetCellSizeX() const { return fCellSizeX; } // cm
    G4double GetCellSizeY() const { return fCellSizeY; } // cm
    G4double GetCellSizeZ() const { return fCellSizeZ; } // cm
    G4int GetNCellsXY() const { return fNCellsXY; }
    G4int GetNCellsZ()  const { return fNCellsZ; }

  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;

    G4VPhysicalVolume* fAbsorberPV = nullptr;
    G4VPhysicalVolume* fGapPV = nullptr;

    G4bool fCheckOverlaps = true;

    G4double fCellSizeX; // cm
    G4double fCellSizeY; // cm
    G4double fCellSizeZ; // cm
    G4String fMaterialType; // "PbF2" or "PbWO4"
    G4int fNCellsXY = 5;
    G4int fNCellsZ  = 5;
};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const {
  return fAbsorberPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetGapPV() const  {
  return fGapPV;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

