/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::patchInteractionData

Description
    Helper class for the LocalInteraction patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef myPatchInteractionData_H
#define myPatchInteractionData_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class patchInteractionData Declaration
\*---------------------------------------------------------------------------*/

// Forward declaration of classes
class myPatchInteractionData;

// Forward declaration of friend functions
Istream& operator>>
(
    Istream& is,
    myPatchInteractionData& pid
);


class myPatchInteractionData
{
    // Private data

        //- Interaction type name
        word interactionTypeName_;

        //- Patch name
        word patchName_;

        //- characteristic Elasticity coefficient
        scalar eH_;

        //- statuc Restitution coefficient
        scalar e0_;

        //- characteristic Restitution coefficient
        scalar muH_;

	    //- static Restitution coefficient
	    scalar mu0_;

	    //- Elasticity angle
	    scalar alphaE_;

	    //- Restitution angle
	    scalar alphaMu_;

        //- Wall roughness coefficient
        scalar gamma_;

        //- Coefficients for the function of e
        scalar e1_;
        scalar e2_;
        scalar e3_;

        //- Coefficients for the function of mu
        scalar mu1_;


public:

    // Constructor

        //- Construct null
        myPatchInteractionData();


    // Member functions

        // Access

            //- Return const access to the interaction type name
            const word& interactionTypeName() const;

            //- Return const access to the patch name
            const word& patchName() const;

            //- Return const access to the characteristic elasticity coefficient
            scalar eH() const;

            //- Return const access to the static elasticity coefficient
            scalar e0() const;

            //- Return const access to the characteristic restitution coefficient
            scalar muH() const;

	        //- Return const access to the static restitution coefficient
	        scalar mu0() const;

	        //- Return const access to the elasticity angle
	        scalar alphaE() const;
	    
	        //- Return const access to the restitution angle
	        scalar alphaMu() const;

            //- Return const access to the roughness coefficient
            scalar gamma() const;

            //- Return const access to the coefficients for e
            scalar e1() const;
            scalar e2() const;
            scalar e3() const;

            //- Return const access to the coefficients for mu
            scalar mu1() const;


        // I-O

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                myPatchInteractionData& pid
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
