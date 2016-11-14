/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________ 
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |                                              
      |_|                                              
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    class for immersed bodies representation.

SourceFiles
    immersedBodies.C

Contributors
    Federico Municchi (2016)
\*---------------------------------------------------------------------------*/

#ifndef immersedBody_H
#define immersedBody_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "triSurfaceMesh.H"
#include <vector>

namespace Foam
{
 enum {STATICBODY , ROTATINGBODY};
 
 class immersedBody
 {
  
  private:

   //first time update?
   bool                                       isFirstUpdate_; 
   
   //Dictionary
   dictionary                                  immersedDict_;
   
   //Mesh
   const  fvMesh&                                      mesh_;
   
   //Operation to perform
   int                                        bodyOperation_;
         
   //base geometry
   triSurfaceMesh*                             bodySurfMesh_;
   
   //surface cells
   std::vector< label >                           surfCells_;
   
   //internal cells
   std::vector< label >                            intCells_;
  
   //interpolation points, point[0] is a surface point
   std::vector< std::vector< point > >  interpolationPoints_;
   
   //interpolation cells
   std::vector< std::vector< label > >   interpolationCells_;
  
   //Create immersed body info
   void  createImmersedBody( volScalarField& body );

   //Update immersed body info
   void  updateImmersedBody( volScalarField& body );
   
   //Create interpolation points
   void  calculateInterpolationPoints(volScalarField& body,  triSurfaceSearch * ibTriSurfSearch );
   
   //Rotate immersed body
   void rotateImmersedBody();
  
   
  public:
   
   immersedBody(word fileName, const Foam::fvMesh& mesh, dictionary HFDIBDict );
   ~immersedBody();
   
  //Getters
    
   const std::vector< label > *     getSurfaceCellList() const {return &surfCells_;};
   
   const std::vector< std::vector< point > > * getInterpolationPoints() const {return &interpolationPoints_;};  
   
   const std::vector< std::vector< label > > * getInterpolationCells() const {return &interpolationCells_;};  
   
   //Update body field
   void  updateBodyField( volScalarField & body );
  
   //Update imposed vector field
   void updateVectoField(volVectorField & VS, word Vname);

   
 };
}
#endif
