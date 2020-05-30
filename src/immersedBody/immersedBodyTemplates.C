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

//Update imposed field
template<class Type>
void immersedBody::updateImposedField
(
    Type& VS,
    word  Vname
)
{
    //- Check dictionary for parameters 
    word BC = immersedDict_.subDict(Vname).lookup("BC");

    //- Implementation of the noSlip boundary condition
    if(BC=="noSlip")
    {
   
        //- Vector fields only!!
        if( Type::typeName != Foam::vector::typeName)
        {
            FatalErrorIn("void immersedBody::updateImposedField(...)")
                << "Attempted to cast noSlip boundary condition "
                << " on a " << Type::typeName << " field." << endl
                << " noSlip can only be applied to vector fields!" << endl;

        }
        
        //- Rotation and translation of the immersed body according to 
        //  centre-of-mass velocity and body angular velocity/axis-of-rotation
        for
        (   
            int cell=0;
            cell<intCells_.size()+surfCells_.size();
            cell++
        )
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

            vector planarVec  =  mesh_.C()[cellI] - CoM_
                                - Axis_
                                *(
                                    (mesh_.C()[cellI]-CoM_)&Axis_
                                 );

            vector VSvalue = (planarVec^Axis_)*omega_ + Vel_;
            VS[cellI] = VSvalue;

        }

    }

}

/*---------------------------------------------------------------------------*/

template<class Type>
void immersedBody::remoteInterpolateField
(
    autoPtr<interpolation<Type> >&  interpField,
    List<Type>&  interpValues      
)
{
    //- Loop over all the remote cells, perform interpolation
    //  and insert into interpValues
    forAll(remoteCells_,rCellI)
    {
        label remotePointId = remoteCells_[rCellI][0];
        label remoteCellId  = remoteCells_[rCellI][1];

        interpValues(remotePointId) =
            interpField->interpolate
            (
                remotePoints_[remotePointId],
                remoteCellId
            );
    }

    //- Communicate remote values
    reduce(interpValues,sumOp<List<Type> >());

}