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
#include "immersedBody.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody(word fileName,
                           const Foam::fvMesh& mesh,
                           dictionary& HFDIBDict,
                           dictionary& transportProperties
                       )
:
isFirstUpdate_(true),
immersedDict_(HFDIBDict.subDict(fileName)),
mesh_(mesh),
transportProperties_(transportProperties),
M_( 0.),
CoM_( vector::zero),
Axis_(vector::one),
omega_( 0.),
Vel_( vector::zero),
I_( symmTensor::zero)
{

 //Read stl file from folder "constant"
 bodySurfMesh_ = new triSurfaceMesh
 (
    IOobject
    (
        fileName +".stl",
        "constant",
        "triSurface",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
 );

 Info<< "Read Immersed Boundary triSurface" << endl;
 bodySurfMesh_->writeStats(Info);

//Return if devlared static
 if(immersedDict_.found("staticBody") )
 {
     bodyOperation_=STATICBODY;
     Info << fileName << " is static body." << endl;
 }
 else if(immersedDict_.found("transRotatingBody"))
 {
     bodyOperation_=TRANSROTATINGBODY;

     //Get basic quantities from dict
     Axis_ = immersedDict_.subDict("rotatingBody").lookup("axis");
     CoM_  = immersedDict_.subDict("rotatingBody").lookup("center");

     omega_  = readScalar(
                    immersedDict_.subDict("rotatingBody").lookup("omega")
                );

     Vel_   = immersedDict_.subDict("rotatingBody").lookup("velocity");


     Info << fileName << " has scripted trans rotational motion." << endl;
 }
 else if(immersedDict_.found("fluidCoupling"))
 {
     bodyOperation_=FLUIDCOUPLING;
     Info << fileName << " is coupled with fluid phase." << endl;

 }
 else
 {
     Info << "No body operation was found for " << fileName << endl
          << "Assuming static body.";
     bodyOperation_=STATICBODY;
 }
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
 delete bodySurfMesh_;
}
//---------------------------------------------------------------------------//
//Update immersed body
void immersedBody::updateBodyField( volScalarField& body,
                                    volVectorField & f
                                 )
{
    if(isFirstUpdate_)
    {
        createImmersedBody( body );
        isFirstUpdate_ = false;
    }
    else
    {
        updateImmersedBody( body, f );
    }
}
//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body )
{

    triSurface ibTemp( *bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
    const pointField & pp = mesh_.points();

    intCells_.clear();
    surfCells_.clear();


    //Fill body field with first estimation
    forAll(mesh_.C(),cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );

        forAll(vertexesInside, verIn)
        {
            if(vertexesInside[verIn]==true)
            {
                body[cellI] += 0.125 ; //fraction of cell covered
                //  Info << "Found vertex inside\n";
            }
        }

        //Add to corresponding vector
        if( body[cellI]>0.9)        intCells_.push_back(cellI);
        else if (body[cellI]>0.1) surfCells_.push_back(cellI);

 }

    //refine body as stated in the dictionary
    refineBody(body,&ibTriSurfSearch, & pp);
    calculateInterpolationPoints(body,&ibTriSurfSearch);
    calculateGeometricalProperties(body);

}
//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody(volScalarField & body,
                                      volVectorField & f
                                  )
{
 //Check Operation to perform

 if(bodyOperation_==STATICBODY) return;
 else
 {
     if(bodyOperation_==FLUIDCOUPLING)
     {
         updateCoupling(body,f);
     }

     resetBody(body);
     moveImmersedBody();
 }

 //TODO:Very inefficient find other algorithm
 createImmersedBody(body);

}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling(volScalarField & body,
                                  volVectorField & f
                              )
{
  const uniformDimensionedVectorField g =
   mesh_.lookupObject<uniformDimensionedVectorField>("g");

  const dimensionedScalar rhof =
   dimensionedScalar(
          transportProperties_.lookup("rho")
      );

  dimensionedScalar rho_ =
       dimensionedScalar(immersedDict_.lookup("rho"));

  vector F(vector::zero);
  vector T(vector::zero);

  //Calcualate viscous force and torque
  for(unsigned int cell=0;cell<surfCells_.size();cell++)
  {
     label cellI = surfCells_[cell];

     F +=  f[cellI] * mesh_.V()[cellI] * rhof.value();
     T +=  ( (mesh_.C()[cellI]-CoM_) ^ f[cellI] )
                 *mesh_.V()[cellI]* rhof.value();

  }

  //Update body linear velocity
  Vel_ += mesh_.time().deltaT().value()
            * (
                  F / M_
                + (1.0-rhof.value()/rho_.value())*g.value()
              );

  //Update body angular velocity
  vector Omega_(vector::zero);

  Omega_ =   Axis_*omega_
           + mesh_.time().deltaT().value() * ( inv(I_) & T );

  //Split Omega_ into Axis_ and omega_
  omega_ = mag(Omega_);

  if(omega_ < 1e-16)
  {
     Axis_ = vector::one;
  }
  else
  {
     Axis_ =  Omega_/omega_;
  }


}
//---------------------------------------------------------------------------//
//Create interpolation points
void
immersedBody::calculateInterpolationPoints(volScalarField& body,
                                           triSurfaceSearch * ibTriSurfSearch
                                         )
{
  double sqrtThree_ = sqrt(3.0);
  meshSearch search_(mesh_);
  //clear previous
  interpolationPoints_.clear();
  interpolationCells_.clear();

  //Create temporary surface normals
  volVectorField surfNorm(-fvc::grad(body));

  for(unsigned int cell=0;cell<surfCells_.size();cell++)
  {
   //get surface cell label
   label scell = surfCells_[cell];

   //create vector for points and cells and add to main vectors
   std::vector<point> intPoints;
   interpolationPoints_.push_back(intPoints);

   std::vector<label> intCells;
   interpolationCells_.push_back(intCells);

   //Get interpolation distance
   double intDist = sqrtThree_*std::pow( mesh_.V()[scell] , 1.0/3.0 );
   vector intVec = surfNorm[scell] * (intDist)/(mag(surfNorm[scell]));
  // Info << "\n intDist: " << intDist << " cellV: " << mesh_.V()[scell];

   //Approximate distance using body
   point surfPoint = mesh_.C()[scell] + intVec*(0.5-body[scell]);

   //Add to list
   interpolationPoints_[cell].push_back(surfPoint);
//    Info << "\nsurfDist = " << mag(surfPoint-mesh_.C()[scell]);

   //Add other interpolation points
   for(int order=0;order<ORDER;order++)
   {
    surfPoint = surfPoint + intVec;
    interpolationPoints_[cell].push_back(surfPoint);
   }

   //Get cells
   for(int order=0;order<ORDER;order++)
   {
    //Use findNearestCell()...it is faster
    label cellI =
        search_.findNearestCell(interpolationPoints_[cell][order+1]);
    //Check if outside the domain (if not, then set to -1)
    if(!(search_.isInside(interpolationPoints_[cell][order+1])) ) cellI = -1;
    interpolationCells_[cell].push_back(cellI);
   }

  }

}
//---------------------------------------------------------------------------//
void immersedBody::calculateGeometricalProperties( volScalarField& body )
{
  //Get density
  dimensionedScalar rho =
       dimensionedScalar(immersedDict_.lookup("rho"));

  //Evaluate center of mass
  M_ = 0.;
  vector tmpCom(vector::zero);
  CoM_ = vector::zero;
  I_ = symmTensor::zero;

  for(unsigned int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
  {
     label cellI;

     if(cell<intCells_.size())
     {
         cellI = intCells_[cell];
     }
     else
     {
         cellI = surfCells_[cell-intCells_.size()];
     }

     M_ += body[cellI] * rho.value() * mesh_.V()[cellI];
     tmpCom  += body[cellI] * rho.value() * mesh_.V()[cellI] * mesh_.C()[cellI];

     I_.xx() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                   + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
               );

     I_.yy() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
                   + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
               );
     I_.zz() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                   + mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
               );

     I_.xy() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].y()
                 );

     I_.xz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].z()
                 );

     I_.yz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].z()
                 );

  }

  //Collect from processors
  reduce(M_, sumOp<scalar>());
  reduce(tmpCom,  sumOp<vector>());
  reduce(I_,  sumOp<symmTensor>());

  CoM_ = tmpCom / M_;
}
//---------------------------------------------------------------------------//
//Rotate immersed body
void immersedBody::moveImmersedBody()
{


 //Rotation angle
 scalar angle = omega_*mesh_.time().deltaT().value();
 vector transIncr = Vel_*mesh_.time().deltaT().value();

 pointField bodyPoints (bodySurfMesh_->points());

 //Move points
 forAll(bodyPoints,p)
 {
  vector normV = (bodyPoints[p]-CoM_)^Axis_;
  scalar axisMod = ((bodyPoints[p]-CoM_)&Axis_);
  scalar magDist = mag( (bodyPoints[p]-CoM_) - Axis_*axisMod/mag(Axis_) );

  //Move in tangential direction
  bodyPoints[p] = bodyPoints[p] +angle*normV;

  //Rescale
  bodyPoints[p] = (bodyPoints[p] - Axis_*((bodyPoints[p]-CoM_)&Axis_))
                  * magDist /mag(
                         bodyPoints[p] - Axis_*((bodyPoints[p]-CoM_)&Axis_)
                     ) +
                   Axis_*((bodyPoints[p]-CoM_)&Axis_);

  //Translate point
  bodyPoints[p] += transIncr;

 }

 //move mesh
 bodySurfMesh_->movePoints(bodyPoints);


}
//---------------------------------------------------------------------------//
//Update imposed vector field
void immersedBody::updateVectorField(volVectorField & VS, word Vname)
{
  //Check dictionary for parameters (only non slip allowed)

  word BC = immersedDict_.subDict(Vname).lookup("BC");

  if(BC=="noSlip")
  {
   //If STATICBODY set to zero
   if( bodyOperation_==STATICBODY)
   {

    for(unsigned int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
        label cellI;

        if(cell<intCells_.size())
        {
            cellI = intCells_[cell];
        }
        else
        {
            cellI = surfCells_[cell-intCells_.size()];
        }

        VS[cellI] = Vel_;

    }

   }
   else
   {

    //Apply
    for(unsigned int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
         label cellI;

         if(cell<intCells_.size())
         {
             cellI = intCells_[cell];
         }
         else
         {
             cellI = surfCells_[cell-intCells_.size()];
         }

         vector planarVec       =  mesh_.C()[cellI] - CoM_
                                 - Axis_*(
                                      (mesh_.C()[cellI]-CoM_)&Axis_
                                     );

         vector VSvalue = (planarVec^Axis_)*omega_ + Vel_;
         VS[cellI] = VSvalue;
    }


   }

  }

}
//---------------------------------------------------------------------------//
//Reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body)
{

   //Simply loop over all the cells and set to zero
   for(unsigned int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
      label cellI;

      if(cell<intCells_.size())
      {
          cellI = intCells_[cell];
      }
      else
      {
          cellI = surfCells_[cell-intCells_.size()];
      }

      body[cellI] = 0.0;
    }

}
//---------------------------------------------------------------------------//
//Refine body field for this immersed object using MC-like algorithm
//Cells are assumed to be hexahedral at the particle surface
//(but can have different edge length)
void immersedBody::refineBody(volScalarField& body,
                              triSurfaceSearch * ibTriSurfSearch,
                              const pointField * pp
                           )
{
    if(!immersedDict_.found("refineMC"))
     return;


    scalar nPointsEdge = readScalar(immersedDict_.lookup("refineMC"));

    //loop over all the surface cells
    for(unsigned int cell=0;cell<surfCells_.size();cell++)
    {

     label cellI = surfCells_[cell];

     scalar deltaV = 1.0/(nPointsEdge*nPointsEdge*nPointsEdge);

     //Get cell center
     point centerC = mesh_.C()[cellI];

     //Get one node
     //Check if partially or completely inside
     const labelList& vertexLabels = mesh_.cellPoints()[cellI];
     const pointField vertexPoints(*pp,vertexLabels);
     point baseNode = vertexPoints[0];

     //create vector representing 3d diagonal of the cell
     vector edgesC = 2*(centerC - baseNode);

     //create list of points
     List<point> pointsMC;

     //create deltas
     scalar delta_i = edgesC[0]/nPointsEdge;
     scalar delta_j = edgesC[1]/nPointsEdge;
     scalar delta_k = edgesC[2]/nPointsEdge;

     //add points to list
     for(int i=0;i<nPointsEdge;i++)
     {
      //point i-coordinate
      scalar icoord = baseNode[0] + delta_i*(i+0.5);

      for(int j=0;j<nPointsEdge;j++)
      {
       //point j-coordinate
       scalar jcoord = baseNode[1] + delta_j*(j+0.5);

       for(int k=0;k<nPointsEdge;k++)
       {
        //point k-coordinate
        scalar kcoord = baseNode[2] + delta_k*(k+0.5);

        //create point
        point p(icoord,jcoord,kcoord);

        //add to list
        pointsMC.append(p);

       }
      }
     }

     //Check who is inside
     pointField pField(pointsMC);
     boolList pInside = ibTriSurfSearch->calcInside( pField );

     //Calculate new body
     scalar newbody = 0.0;
     forAll(pInside,p)
     {
      if(pInside[p])
       newbody+=deltaV;


     }


     body[cellI] = newbody;

    }


}
