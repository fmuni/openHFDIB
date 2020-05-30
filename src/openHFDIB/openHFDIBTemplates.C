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

template<class Type>
void openHFDIB::interpolateIB
(
    GeometricField<Type,fvPatchField,volMesh>&   V,
    GeometricField<Type,fvPatchField,volMesh>& Vs,
    volScalarField& body
)
{
    //Create interpolator
    autoPtr<interpolation<Type> > interpV =
        interpolation<Type>::New(HFDIBinterpDict_, V);


    Type zeros = pTraits<Type>::zero;
 
    //Reset imposed field
    forAll(Vs,cellI)
    {
        Vs[cellI] =zeros;
    }

    //Loop over all the immersed bodies
    forAll(immersedBodies_,bodyId)
    {
        
        //Update imposed field according to body
        immersedBodies_[bodyId]->updateImposedField(Vs, V.name());

        const labelList&  surCells  =
            immersedBodies_[bodyId]->getSurfaceCellList();

        const List<pointField>& intPoints =
            immersedBodies_[bodyId]->getInterpolationPoints();

        const List<labelList>& intCells  =
            immersedBodies_[bodyId]->getInterpolationCells();

         //loop over all surface cells
        forAll(surCells,scell)
        {

            label cellI = surCells[scell];
            
            //Check max order of accuracy
            bool allowedOrder[ORDER];

            for(int intPoint=0;intPoint<ORDER;intPoint++)
            {
                if( intCells[scell][intPoint] == -1 )
                {
                     allowedOrder[intPoint] = false;
                }                
                else 
                {
                    allowedOrder[intPoint] = true;
                }
            }

            bool  firstOrder = false;
            bool secondOrder = true;
            bool   zeroOrder = false;

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

                Type VP1 =
                    interpV->interpolate
                    (  
                        intPoints[scell][1],
                        intCells[scell][0]
                    ) 
                    - Vs[cellI];

                Type VP2 =
                   interpV->interpolate
                    (  
                        intPoints[scell][2],
                        intCells[scell][1]
                    )
                    - Vs[cellI];


                //distance between interpolation points
                scalar res_ =
                    mag(intPoints[scell][2] - intPoints[scell][1]);

                //cell centre to surface distance
                scalar ds   = res_*(0.5-body[cellI ]) ;

                //Coefficient of the quadratic term
                Type quadCoeff = 1/(res_*res_) * ( VP2/2 - VP1 );

                //Coefficient of the linear term
                Type linCoeff  = 1/(2*res_) * ( 4*VP1 - VP2 );

                Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI];
            }
            else if(firstOrder)
            {
               
                Type VP1 =
                    interpV->interpolate
                    (  
                        intPoints[scell][1],
                        intCells[scell][0]
                    ) 
                    - Vs[cellI];

                //distance between interpolation points
                scalar res_ =
                    mag(intPoints[scell][2] -intPoints[scell][1]);

                //cell centre to surface distance
                scalar ds   = res_*(0.5-body[cellI ]) ;

                //Coefficient of the linear term
                Type linCoeff = VP1/res_;

                Vs[cellI] = linCoeff*ds + Vs[cellI];
            }
            else if(zeroOrder)
            {
                
                //Zero order is pure Fictitious Domain
                Vs[cellI] = body[cellI]*Vs[cellI] + (1.0-body[cellI])*V[cellI];
            }

        }

    }
}
