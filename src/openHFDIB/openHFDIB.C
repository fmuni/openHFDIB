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

Contributors
    Federico Municchi (2018)
\*---------------------------------------------------------------------------*/
#include "openHFDIB.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"


using namespace Foam;

//---------------------------------------------------------------------------//
openHFDIB::openHFDIB(const Foam::fvMesh& mesh )
:
mesh_(mesh),
HFDIBDict_
 (
    IOobject
    (
        "HFDIBDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
 ),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
 ),
immersedBodies_(0)
{
}
//---------------------------------------------------------------------------//
openHFDIB::~openHFDIB()
{
    forAll(immersedBodies_,bodyI)
    {
        immersedBodies_[bodyI].clear();
    }  

    immersedBodies_.clear();
}
//---------------------------------------------------------------------------//
void openHFDIB::initialize()
{

    wordList stlNames( HFDIBDict_.lookup("bodyList") );

    HFDIBinterpDict_ = HFDIBDict_.subDict("interpolationSchemes");
   
    //Generate immersed objects    
    forAll(stlNames,name)
    {
        immersedBodies_.append
        ( 
            autoPtr<immersedBody>
            (
                new immersedBody
                (  
                    stlNames[name],
                    mesh_,
                    HFDIBDict_,
                    transportProperties_
                )
            )
        );
    }

}
//---------------------------------------------------------------------------//
void openHFDIB::update
(
    volScalarField & body,
    volVectorField & f
)
{

    forAll(immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId]->updateBodyField(body,f);
    }

}
//---------------------------------------------------------------------------//