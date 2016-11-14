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

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody(word fileName, const Foam::fvMesh& mesh,dictionary HFDIBDict )
:
isFirstUpdate_(true),
immersedDict_(HFDIBDict.subDict(fileName)),
mesh_(mesh)
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
 if(immersedDict_.found("staticBody") )            bodyOperation_=STATICBODY;
 else if(immersedDict_.found("rotatingBody"))   bodyOperation_=ROTATINGBODY;
 else
 {
  Info << "No body operation was found for " << fileName <<". Assuming static body.";
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
void immersedBody::updateBodyField(volScalarField& body )
{
 if(isFirstUpdate_)
 {
  createImmersedBody( body );
  isFirstUpdate_ = false;
 }
 else
  updateImmersedBody( body ); 
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


 //Fill body field 
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
 
 calculateInterpolationPoints(body,&ibTriSurfSearch);

}
//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody(volScalarField& body )
{
 //Check Operation to perform
 
 if(bodyOperation_==STATICBODY) return;
 else if(bodyOperation_==ROTATINGBODY) rotateImmersedBody();

 //TODO:Very inefficient find other algorithm
 createImmersedBody(body);
 
}
//---------------------------------------------------------------------------//
//Create interpolation points
void immersedBody::calculateInterpolationPoints(volScalarField& body,
                                                triSurfaceSearch * ibTriSurfSearch
                                               )
{
  double sqrtThree_ = sqrt(3.0); 
  meshSearch search_(mesh_);
  //clear previous
  interpolationPoints_.clear();
  interpolationCells_.clear();
  
  //Create temporary surface normals
  volVectorField surfNorm = -fvc::grad(body);
     
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
   point surfPoint = mesh_.C()[scell] + intVec*(0.5-body[scell])/2.0;
   
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
    label cellI =  search_.findNearestCell(interpolationPoints_[cell][order+1]);
    //Check if outside the domain (if not, then set to -1)
    if(!(search_.isInside(interpolationPoints_[cell][order+1])) ) cellI = -1;
    interpolationCells_[cell].push_back(cellI);
   }
  
  }
 
}
//---------------------------------------------------------------------------//
//Rotate immersed body
void immersedBody::rotateImmersedBody()
{
 
 //Get basic quantities from dict
 vector axis   = immersedDict_.subDict("rotatingBody").lookup("axis");
 point  center = immersedDict_.subDict("rotatingBody").lookup("center");
 scalar omega  = readScalar(immersedDict_.subDict("rotatingBody").lookup("omega"));
 
 //get the angle 
 scalar angle  = omega*mesh_.time().deltaT().value();
 
 pointField bodyPoints = bodySurfMesh_->points();
 
 //Move points
 forAll(bodyPoints,p)
 {
  vector normV = (bodyPoints[p]-center)^axis;
  scalar axisMod = ((bodyPoints[p]-center)&axis);
  scalar magDist = mag( (bodyPoints[p]-center) - axis*axisMod/mag(axis) );
  
  //Move in tangential direction
  bodyPoints[p] = bodyPoints[p] +angle*normV;
  
  //Rescale
  bodyPoints[p] = (bodyPoints[p] - axis*((bodyPoints[p]-center)&axis))
                  * magDist/mag(bodyPoints[p] - axis*((bodyPoints[p]-center)&axis)) + 
                   axis*((bodyPoints[p]-center)&axis);

 
 }
 
 //move mesh
 bodySurfMesh_->movePoints(bodyPoints);
 

}
//---------------------------------------------------------------------------//
//Update imposed vector field
void immersedBody::updateVectoField(volVectorField & VS, word Vname)
{
  //Check dictionary for parameters (only non slip allowed)
  
  word BC = immersedDict_.subDict(Vname).lookup("BC"); 
  
  if(BC=="noSlip")
  {
   //If STATICBODY set to zero
   if( bodyOperation_==STATICBODY)
   {
     vector VSvalue=vector::zero;
    for(unsigned int cell=0;cell<intCells_.size();cell++)
    {
    
     label cellI = intCells_[cell];
     VS[cellI] = VSvalue;
     
    }

    for(unsigned int cell=0;cell<surfCells_.size();cell++)
    {
     label cellI = surfCells_[cell];
     VS[cellI] = VSvalue;
     
    }
   
   }
   
  //If ROTATINGBODY apply periferical velocity
   if( bodyOperation_==ROTATINGBODY)
   {

    //Get basic quantities from dict
    vector axis   = immersedDict_.subDict("rotatingBody").lookup("axis");
    point  center = immersedDict_.subDict("rotatingBody").lookup("center");
    scalar omega  = readScalar(immersedDict_.subDict("rotatingBody").lookup("omega"));
     
    //Apply for internal cells   
    for(unsigned int cell=0;cell<intCells_.size();cell++)
    {
     label cellI            = intCells_[cell];
     vector planarVec       = mesh_.C()[cellI]-center - axis*((mesh_.C()[cellI]-center)&axis);
     vector VSvalue = (planarVec^axis)*omega;
     VS[cellI] = VSvalue;
     
    }

    //Apply for surface cells (here should apply the surface value)   
    for(unsigned int cell=0;cell<surfCells_.size();cell++)
    {
     label cellI            = surfCells_[cell];
     point surfPoint        = interpolationPoints_[cell][0];
     vector planarVec       = surfPoint - center - axis*((surfPoint-center)&axis);
     vector VSvalue = planarVec^axis*omega;
     VS[cellI] = VSvalue;
     
    }

   
   }
  
  }

}
