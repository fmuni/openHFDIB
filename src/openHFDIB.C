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
    Federico Municchi (2016)
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

#define ORDER 2

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
            "constant",
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
           "constant",
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       )
 )
{
}
//---------------------------------------------------------------------------//
openHFDIB::~openHFDIB()
{
 for(unsigned int i=0;i<immersedBodies_.size();i++)
  delete immersedBodies_[i];

 immersedBodies_.clear();
}
//---------------------------------------------------------------------------//
void openHFDIB::initialize()
{

 wordList stlNames( HFDIBDict_.lookup("stlNames") );

 HFDIBinterpDict_ = HFDIBDict_.subDict("interpolationSchemes");
 //Generate immersed objects
 forAll(stlNames,name)
 {
  immersedBody * body_ptr =
               new immersedBody(  stlNames[name],
                                  mesh_,
                                  HFDIBDict_,
                                  transportProperties_
                                );

  immersedBodies_.push_back(body_ptr);

 }


}
//---------------------------------------------------------------------------//
void openHFDIB::update(volScalarField & body,
                       volVectorField & f
            )
{

 for(unsigned int bodyId=0;bodyId<immersedBodies_.size();bodyId++)
 {
    immersedBodies_[bodyId]->updateBodyField(body,f);
 }
}
//---------------------------------------------------------------------------//
void openHFDIB::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{

 //Create interpolator
 autoPtr<interpolation<vector> > interpV =
                   interpolation<vector>::New(HFDIBinterpDict_, V);



 vector zeros = vector::zero;
 //Reset imposed field
 forAll(Vs,cellI)
 {
  Vs[cellI] =zeros;
 }

 //Loop over all the immersed bodies
 for(unsigned int bodyId=0;bodyId<immersedBodies_.size();bodyId++)
 {
  //Update imposed field according to body
  immersedBodies_[bodyId]->updateVectorField(Vs, V.name());

  const std::vector<label> *  surCells  =
                      immersedBodies_[bodyId]->getSurfaceCellList();

  const std::vector< std::vector< point > > * intPoints =
                      immersedBodies_[bodyId]->getInterpolationPoints();

  const std::vector< std::vector< label > > * intCells  =
                      immersedBodies_[bodyId]->getInterpolationCells();

  //loop over all surface cells
  for(unsigned int scell=0;scell<surCells->size();scell++)
  {

   label cellI = (*surCells)[scell];
   //Check max order of accuracy
   bool allowedOrder[ORDER];

   for(int intPoint=0;intPoint<ORDER;intPoint++)
    if( (*intCells)[scell][intPoint] == -1 ) allowedOrder[intPoint] = false;
    else allowedOrder[intPoint] = true;

   bool  firstOrder = false;
   bool secondOrder = true;
   bool zeroOrder   = false;

   //Check is second order is possible
   if( allowedOrder[1] == false)
   {
    secondOrder = false;
    firstOrder  = true;
   }

   //Check if first order is possible
   if( allowedOrder[0] == false)
   {
    secondOrder = false;
    firstOrder  = false;
    zeroOrder   = true;
   }

   //Go for interpolation!
   if(secondOrder)
   {

     vector VP1 =  interpV->interpolate(  (*intPoints)[scell][1],
                                          (*intCells)[scell][0]
                                        ) - Vs[cellI];

     vector VP2 =  interpV->interpolate(  (*intPoints)[scell][2],
                                          (*intCells)[scell][1]
                                        ) - Vs[cellI];


    //distance between interpolation points
    double res_ = mag((*intPoints)[scell][2] -(*intPoints)[scell][1]);

    //cell center to surface distance
   double ds   = res_*(0.5-body[cellI ]) ;

    vector quadCoeff = 1/(res_*res_) * ( VP2/2 - VP1 );
    vector linCoeff  = 1/(2*res_) * ( 4*VP1 - VP2 );

    Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI]  ;
   }
   else if(firstOrder)
   {
     vector VP1 =  interpV->interpolate(  (*intPoints)[scell][1],
                                          (*intCells)[scell][0]
                                        ) - Vs[cellI];




    //distance between interpolation points
    double res_ = mag((*intPoints)[scell][1] -(*intPoints)[scell][0]);

    //cell center to surface distance
    double ds   = res_*(0.5-body[cellI]) ;

    vector linCoeff = VP1/res_;

    Vs[cellI] = linCoeff*ds + Vs[cellI];
   }
   else if(zeroOrder)
   {
    //Zero order
    Vs[cellI] = body[cellI]*Vs[cellI]  + (1.0-body[cellI])*V[cellI];
   }

  }

 }
}
