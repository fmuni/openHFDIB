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
#include "UPstream.H"
#include "transformGeometricField.H"


#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word bodyName,
    const Foam::fvMesh& mesh,
    dictionary& HFDIBDict,
    dictionary& transportProperties
)
:
bodyName_(bodyName),
isFirstUpdate_(true),
immersedDict_(HFDIBDict.subDict(bodyName)),
mesh_(mesh),
transportProperties_(transportProperties),
M_(0.),
CoM_(vector::zero),
Axis_(vector::zero),
omega_(0.),
Vel_(vector::zero),
I_(symmTensor::zero),
bodySurfMesh_
(
    new triSurfaceMesh
    (
        IOobject
        (
            immersedDict_.lookup("fileName"),
            mesh_.time().constant(),
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
)
{

    //Process stl file
    Info<< "Read Immersed Boundary triSurface "
        << immersedDict_.lookup("fileName")
        << " for body " << bodyName << endl;

    bodySurfMesh_->writeStats(Info);
    Info << endl;

    if(immersedDict_.found("transform"))
    {
        dictionary transformDict = immersedDict_.subDict("transform");

        transformBody(transformDict);
    }

    //Return if declared static
    if(immersedDict_.found("staticBody") )
    {
        bodyOperation_=STATICBODY;
        Info << bodyName << " is static body." << endl;
    }
    else if(immersedDict_.found("transRotatingBody"))
    {
        bodyOperation_=TRANSROTATINGBODY;

        //Get basic quantities from dict
        Axis_ = immersedDict_.subDict("transRotatingBody").lookup("axis");
        CoM_  = immersedDict_.subDict("transRotatingBody").lookup("center");

        omega_  = readScalar
        (
            immersedDict_.subDict("transRotatingBody").lookup("omega")
        );

        Vel_ = immersedDict_.subDict("transRotatingBody").lookup("velocity");

        Info << bodyName << " has scripted trans rotational motion." << endl;
    }
    else if(immersedDict_.found("fluidCoupling"))
    {
        bodyOperation_=FLUIDCOUPLING;
        Info << bodyName << " is coupled with fluid phase." << endl;

    }
    else
    {
        Info << "No body operation was found for " << bodyName << endl
             << "Assuming static body.";
        bodyOperation_=STATICBODY;
    }
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
    bodySurfMesh_.clear();
}
//---------------------------------------------------------------------------//
void immersedBody::transformBody(dictionary& transformDict)
{

    Info << "Transforming immersed body " << bodyName_
         << " using dictionary" << endl;

    pointField bodyPoints = bodySurfMesh_->points();

    vector transVec = transformDict.lookupOrDefault<vector>
    (
        "translate",
        vector::zero
    );

//- TODO: Implement!
//    vector rotVec = transformDict.lookupOrDefault<vector>
//    (
//        "rotate",
//        vector::zero
//    );
//
//    vector scaleVec = transformDict.lookupOrDefault<vector>
//    (
//        "scale",
//        vector::one
//    );

    forAll(bodyPoints,pointI)
    {
        bodyPoints[pointI] += transVec;
    }

     bodySurfMesh_->movePoints(bodyPoints);

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

    body.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body )
{

    triSurface ibTemp(bodySurfMesh_());
    triSurfaceSearch ibTriSurfSearch(ibTemp);
    const pointField & pp = mesh_.points();

    intCells_.clear();
    surfCells_.clear();
    surfNorm_.clear();

    //Info<<"\n Creating immersed body\n";

    //Fill body field with first estimation
    forAll(mesh_.C(),cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside(vertexPoints);

        bool bodyCell(false);
        body[cellI] = 0.;
        
        label nIn(0);
        label nOut(0);
        vector centerIn(vector::zero);
        vector centerOut(vector::zero);

        forAll(vertexesInside, verIn)
        {
            if(vertexesInside[verIn]==true)
            {
                //fraction of cell covered
                body[cellI] += 1.0/(vertexPoints.size());
               // Info << "Found vertex inside\n";
                bodyCell = true;
                
                centerIn += vertexPoints[verIn];
                nIn++;
            }
            else
            {
                centerOut += vertexPoints[verIn];
                nOut++;
            }
            
        }

        //Add to corresponding vector
        if(bodyCell)
        {
            if( body[cellI]>0.9)
            {
                intCells_.append(cellI);
            }        
            else if (body[cellI]>0.1)
            {
                surfCells_.append(cellI);
                
                //- Now evaluate the normal
                vector surfN((centerOut/max(nOut,1)) - (centerIn/max(nIn,1)));
                //Info<<"\n mag(surfN)\n";
 /*               surfNorm_.append
                (
                    surfN/(mag(surfN)+small)
                ); */
            }   
        }
    }

    //refine body as stated in the dictionary
    refineBody(body,ibTriSurfSearch,pp);
    calculateInterpolationPoints(body,ibTriSurfSearch);
    calculateGeometricalProperties(body);

}
//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody
(
    volScalarField & body,
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

    //TODO:inefficient algorithm...
    createImmersedBody(body);

}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField & body,
    volVectorField & f
)
{

    const uniformDimensionedVectorField g =
    mesh_.lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar rhof(transportProperties_.lookup("rho"));

    dimensionedScalar rho_(immersedDict_.lookup("rho"));

    vector F(vector::zero);
    vector T(vector::zero);

    //Calcualate viscous force and torque
    forAll(surfCells_,cell)
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

    if(omega_ < 1e-32)
    {
        Axis_ = vector::zero;
    }
    else
    {
       Axis_ =  Omega_/omega_;
    }


}
//---------------------------------------------------------------------------//
//Create interpolation points
void
immersedBody::calculateInterpolationPoints
(
    volScalarField& body,
    triSurfaceSearch& ibTriSurfSearch
                                         )
{

    double sqrtThree_ = sqrt(3.0);
    meshSearch search_(mesh_);

    //clear previous
    interpolationPoints_.clear();
    interpolationCells_.clear();

    //Create temporary surface normals
    volVectorField surfNorm(-fvc::grad(body));
    
    //Create local list of remote nodes
    List<point> localRemotePoints;

    forAll(surfCells_,cell)
    {
        //get surface cell label
        label scell = surfCells_[cell];

        //create vector for points and cells and add to main vectors
        pointField intPoints;
        interpolationPoints_.append(intPoints);

        labelList intCells;
        interpolationCells_.append(intCells);

        //Get interpolation distance
        scalar intDist = sqrtThree_*std::pow( mesh_.V()[scell] , 1.0/3.0 );
         vector intVec = surfNorm[scell] * (intDist)/(mag(surfNorm[scell]));

        //Approximate distance using body
        point surfPoint = mesh_.C()[scell] + intVec*(0.5-body[scell]);

        //Add to list
        interpolationPoints_[cell].append(surfPoint);

       //Add other interpolation points
        for(int order=0;order<ORDER;order++)
        {
            surfPoint = surfPoint + intVec;
            interpolationPoints_[cell].append(surfPoint);
        }

        //Get cells
        for(int order=0;order<ORDER;order++)
        {

            label cellI;

            //Check if inside the domain (if not, then set to -1)
            if
            (
                !search_.isInside
                (
                    interpolationPoints_[cell][order+1]
                )
            )
            {
                cellI = -1;
                localRemotePoints.append(interpolationPoints_[cell][order+1]);
            }
            else
            {
                cellI =
                search_.findCell(interpolationPoints_[cell][order+1]);
            }

            interpolationCells_[cell].append(cellI);
        }

    }

    //- Parallel communication

    if (!UPstream::parRun())
    {
        return;
    }

    //- Standard MPI algorithm in OF language, nothing special.
    //  Uninterested reader can skip the code.
    labelList displ(UPstream::nProcs(),0);
    labelList numfrags(UPstream::nProcs(),0);
    label   remoteSize = 0;

    numfrags[UPstream::myProcNo()] = localRemotePoints.size();

    reduce(numfrags,sumOp<labelList>());

    forAll(displ,proc)
    {
        displ[proc] = remoteSize;
        remoteSize += numfrags[proc];
    }

    remoteDispl_ = displ[UPstream::myProcNo()];

    remotePoints_.clear();
    remotePoints_.resize(remoteSize,vector::zero);

    forAll(localRemotePoints,lrpI)
    {
        remotePoints_[remoteDispl_+lrpI] = localRemotePoints[lrpI];
    }

    reduce(remotePoints_,sumOp<List<point> >());

    //- Now every processor has a list of unfound (-1) intepolation points.
    //  At this point the algorithm establishes who holds what.

    //- Only poistion in remotePoints_ and corresponding cell are stored.
    remoteCells_.clear();

    forAll(remotePoints_,rpI)
    {
        //- Continue loop if outside the domain.
        //  Perhaps someone else has it.
        if(!search_.isInside(remotePoints_[rpI]))
        {
            continue;
        }

        //-  It should be here, run findCell
        label cellI =
        search_.findCell(remotePoints_[rpI]);

        //- Check if valid (it should be), then add to list
        if(cellI>-1)
        {
            labelList tmpList(2,zero());
            tmpList[0] = rpI;
            tmpList[1] = cellI;

            remoteCells_.append(tmpList);
        }

    }

}
//---------------------------------------------------------------------------//
void immersedBody::calculateGeometricalProperties( volScalarField& body )
{

    //Get density
    dimensionedScalar rho(immersedDict_.lookup("rho"));

    //Evaluate center of mass
    M_ = 0.;
    vector tmpCom(vector::zero);
    CoM_ = vector::zero;
    I_ = symmTensor::zero;

    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
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

        tmpCom  +=   body[cellI] * rho.value()
                   * mesh_.V()[cellI] * mesh_.C()[cellI];

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
//Move immersed body according to body operation
void immersedBody::moveImmersedBody()
{


    //Rotation angle
    scalar angle = omega_*mesh_.time().deltaT().value();
    vector transIncr = Vel_*mesh_.time().deltaT().value();

    pointField bodyPoints (bodySurfMesh_->points());
    
    //- Just take one point to compute the reference direction
    vector n1(bodyPoints[0]) ;
    
    //- Create vector normal to axis and radius
    vector normV = (n1-CoM_)^Axis_;

    //Move in tangential direction
    point n2 = bodyPoints[0] +angle*normV;
    
    //- Rescale
    n1 /= mag(n1);
    n2 /= mag(n2);
    
    //- Rotate points from axis n1 to axis n2
    tensor T(rotationTensor(n1, n2));
    
    //- Rotate points
    bodyPoints = transform(T, bodyPoints);

    //- Translate points
    forAll(bodyPoints,p)
    {
        bodyPoints[p] += transIncr;
    }

    //move mesh
    bodySurfMesh_->movePoints(bodyPoints);

}
//---------------------------------------------------------------------------//
//Reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body)
{

    //Simply loop over all the cells and set to zero
    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
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
void immersedBody::refineBody
(
    volScalarField& body,
    triSurfaceSearch& ibTriSurfSearch,
    const pointField& pp
)
{
    if(!immersedDict_.found("refineMC"))
    {
        return;
    }


    scalar nPointsEdge = readScalar(immersedDict_.lookup("refineMC"));

    //loop over all the surface cells
    forAll(surfCells_,cell)
    {

        label cellI = surfCells_[cell];

        scalar deltaV = 1.0/(nPointsEdge*nPointsEdge*nPointsEdge);

        //Get cell center
        point centerC = mesh_.C()[cellI];

        //Get one node
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        point baseNode = vertexPoints[0];

        //create vector representing 3d diagonal of the cell
        vector edgesC = 2.0*(centerC - baseNode);

        //create list of points
        pointField pointsMC;

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
        boolList pInside = ibTriSurfSearch.calcInside( pField );

        //Calculate new body and new normal
        scalar newbody = 0.0;
        label nIn(0);
        label nOut(0);
        vector centerIn(vector::zero);
        vector centerOut(vector::zero);
        
        forAll(pInside,p)
        {
            if(pInside[p])
            {
                newbody+=deltaV;
                centerIn += pointsMC[p];
                nIn++;
            }
            else
            {
                centerOut += pointsMC[p];
                nOut++;         
            }
            
        }


        body[cellI] = newbody;
        
        //- Now re-evaluate normal
        vector surfN( (centerOut/max(nOut,1)) - (centerIn/max(nIn,1)));
        surfNorm_[cell] =
        (
            surfN/(mag(surfN)+small)
        );

    }

}
