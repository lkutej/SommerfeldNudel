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

\*---------------------------------------------------------------------------*/

#include "myPatchInteractionData.H"
#include "dictionaryEntry.H"
#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::myPatchInteractionData::myPatchInteractionData()
:
    interactionTypeName_("unknownInteractionTypeName"),
    patchName_("unknownPatch"),
    eH_(0.0),
    e0_(0.0),
    muH_(0.0),
    mu0_(0.0),
    alphaE_(0.0),
    alphaMu_(0.0),
    gamma_(0.0),
    e1_(0.0),
    e2_(0.0),
    e3_(0.0),
    mu1_(0.0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::myPatchInteractionData::interactionTypeName() const
{
    return interactionTypeName_;
}


const Foam::word& Foam::myPatchInteractionData::patchName() const
{
    return patchName_;
}


Foam::scalar Foam::myPatchInteractionData::eH() const
{
    return eH_;
}

Foam::scalar Foam::myPatchInteractionData::e0() const
{
    return e0_;
}

Foam::scalar Foam::myPatchInteractionData::muH() const
{
    return muH_;
}

Foam::scalar Foam::myPatchInteractionData::mu0() const
{
    return mu0_;
}

Foam::scalar Foam::myPatchInteractionData::alphaE() const
{
    return alphaE_;
}

Foam::scalar Foam::myPatchInteractionData::alphaMu() const
{
    return alphaMu_;
}

Foam::scalar Foam::myPatchInteractionData::gamma() const
{
    return gamma_;
}

Foam::scalar Foam::myPatchInteractionData::e1() const
{
    return e1_;
}

Foam::scalar Foam::myPatchInteractionData::e2() const
{
    return e2_;
}

Foam::scalar Foam::myPatchInteractionData::e3() const
{
    return e3_;
}

Foam::scalar Foam::myPatchInteractionData::mu1() const
{
    return mu1_;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    myPatchInteractionData& pid
)
{
    is.check("Istream& operator>>(Istream&, myPatchInteractionData&)");

    const dictionaryEntry entry(dictionary::null, is);

    pid.patchName_ = entry.keyword();
    entry.lookup("type") >> pid.interactionTypeName_;
    pid.eH_ = entry.lookupOrDefault<scalar>("eH", 1.0);
    pid.e0_ = entry.lookupOrDefault<scalar>("e0", 1.0);
    pid.muH_ = entry.lookupOrDefault<scalar>("muH", 0.0);
    pid.mu0_ = entry.lookupOrDefault<scalar>("mu0", 0.0);
    pid.alphaE_ = entry.lookupOrDefault<scalar>("alphaE", 0.0);
    pid.alphaMu_ = entry.lookupOrDefault<scalar>("alphaMu", 0.0);
    pid.gamma_ = entry.lookupOrDefault<scalar>("gamma",0.0);
    pid.e1_ = entry.lookupOrDefault<scalar>("e1",0.0);
    pid.e2_ = entry.lookupOrDefault<scalar>("e2",0.0);
    pid.e3_ = entry.lookupOrDefault<scalar>("e3",0.0);
    pid.mu1_ = entry.lookupOrDefault<scalar>("mu1",0.0);

//    pid.eH_ = readScalar(entry.lookup("eH"));
//    pid.muH_ = readScalar(entry.lookup("muH"));
//    pid.mu0_ = readScalar(entry.lookup("mu0"));
//    pid.alphaE_ = readScalar(entry.lookup("alphaE"));
//    pid.alphaMu_ = readScalar(entry.lookup("alphaMu"));

    return is;
}


// ************************************************************************* //
