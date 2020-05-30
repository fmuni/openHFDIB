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
    Federico Municchi (2018)
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
#include "polyMesh.H"
#include "fvCFD.H"
#include "List.H"
#include "autoPtr.H"
#include "labelList.H"
#include "UPstream.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"


namespace Foam
{
 enum {
                STATICBODY, //Body does not move
         TRANSROTATINGBODY, //Trans-rotational imposed motion
             FLUIDCOUPLING  //Fluid driven motion
       };

 class immersedBody
 {

    private:

        //- Name of immersed body
        word          bodyName_;                  

        //- First time update?
        bool          isFirstUpdate_;

        //- Dictionary
        dictionary&   immersedDict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Dictionary
        dictionary&   transportProperties_;


        //- Body mass
        scalar M_;

        //- Location of center of mass
        vector CoM_;

        //- Rotation axis
        vector Axis_;

        //- Angular velocity
        scalar omega_;

        //- Body velocity
        vector Vel_;

        //- Moment of intertia tensor
        symmTensor  I_;

        //- Operation to perform
        int   bodyOperation_;

        //- Base geometry
        autoPtr<triSurfaceMesh>             bodySurfMesh_;

        //- Surface cells   
        labelList           surfCells_;

        //- Internal cells
        labelList            intCells_;
        
        //- Surface normals
        vectorField         surfNorm_;

        //- Interpolation points, point[0] is a surface point
        List<pointField>     interpolationPoints_;

        //- Interpolation cells
        List<labelList>              interpolationCells_;

        //- List of remote points (global)
        List<point>        remotePoints_;
    
        //- List of remote cells
        List<labelList>    remoteCells_; 

        //- Processor displacement in remotePoints_
        scalar      remoteDispl_;


        //- Create immersed body info
        void  createImmersedBody( volScalarField& body );

        //- Update immersed body info
        void  updateImmersedBody
        (
            volScalarField& body,
            volVectorField& f
        );

        //- Create interpolation points
        void  calculateInterpolationPoints
        (
            volScalarField& body,
            triSurfaceSearch& ibTriSurfSearch
        );

        //- Evaluate center of mass and inertia tensor
        void calculateGeometricalProperties(volScalarField& body);

        //- Rotate immersed body
        void moveImmersedBody();

        //- Reset body field
        void resetBody(volScalarField& body);

        //- Refine body field using MC-like integration
        void refineBody
        (
            volScalarField& body,
            triSurfaceSearch& ibTriSurfSearch,
            const pointField& pp
        );

        //- Update Fluid-Solid coupling
        void updateCoupling
        ( 
            volScalarField & body,
            volVectorField & f
        );

        //- Tranform body using dictionary
        void transformBody(dictionary& transformDict);


    public:

        immersedBody
        (
            word fileName,
            const Foam::fvMesh& mesh,
            dictionary& HFDIBDict,
            dictionary& transportProperties
        );

        ~immersedBody();

        const labelList& getSurfaceCellList() const 
        {
            return surfCells_;
        };

        const List<pointField>& getInterpolationPoints() const
        {
            return interpolationPoints_;
        };

        const List<labelList>& getInterpolationCells() const
        {
            return  interpolationCells_;
        };

        //- Update body field
        void updateBodyField
        (
            volScalarField & body,
            volVectorField & f
        );

        //- Update imposed field
        template<class Type>
        void updateImposedField(Type&, word);

        //- Remote interpolation
        template<class Type>
        void remoteInterpolateField
        (
            autoPtr<interpolation<Type> >&  interpField,
            List<Type>&  interpValues      
        );

 };

#   include "immersedBodyTemplates.C"

}
#endif
