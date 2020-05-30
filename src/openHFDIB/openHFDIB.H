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
    implementation of the HFDIBM method (Municchi and Radl, 2016) in OpenFOAM

SourceFiles
    openHFDIB.C

Contributors
    Federico Municchi (2018)
\*---------------------------------------------------------------------------*/

#ifndef openHFDIB_H
#define openHFDIB_H

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
#include "List.H"
#include "autoPtr.H"
#include "immersedBody.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#define ORDER 2

namespace Foam
{

class openHFDIB
{

    private:

    const  fvMesh&                   mesh_;
    IOdictionary                HFDIBDict_;
    dictionary            HFDIBinterpDict_;
    IOdictionary      transportProperties_;

    //Container for immersed bodies
    List<autoPtr<immersedBody> > immersedBodies_; 

    public:

    openHFDIB(const Foam::fvMesh& mesh);
     ~openHFDIB();

    void initialize();
    
    void update
    (
        volScalarField & body,
        volVectorField & f
    );

    //Interpolate field for IB
    template<class Type>
    void interpolateIB
    (
        GeometricField<Type,fvPatchField,volMesh>&   V,
        GeometricField<Type,fvPatchField,volMesh>& Vs,
        volScalarField& body
    );

 
 };

    #include "openHFDIBTemplates.C"


}



#endif
