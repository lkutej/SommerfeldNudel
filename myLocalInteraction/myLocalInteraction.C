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

#include "myLocalInteraction.H"
#include "fvcCurl.H"
#include "Random.H"
#include <fstream>
#include <math.h>
#include <iostream>
using namespace std;

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myLocalInteraction<CloudType>::myLocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchData_(cloud.mesh(), this->coeffDict()),
    nEscape_(patchData_.size(), 0),
    massEscape_(patchData_.size(), 0.0),
    nStick_(patchData_.size(), 0),
    massStick_(patchData_.size(), 0.0),
    writeFields_(this->coeffDict().lookupOrDefault("writeFields", false)),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{
    if (writeFields_)
    {
        word massEscapeName(this->owner().name() + ":massEscape");
        word massStickName(this->owner().name() + ":massStick");
        Info<< "    Interaction fields will be written to " << massEscapeName
            << " and " << massStickName << endl;

        (void)massEscape();
        (void)massStick();
    }
    else
    {
        Info<< "    Interaction fields will not be written" << endl;
    }

    // check that interactions are valid/specified
    forAll(patchData_, patchi)
    {
        const word& interactionTypeName =
            patchData_[patchi].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchData_[patchi].patchName();
            FatalErrorInFunction
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }
    }
}


template<class CloudType>
Foam::myLocalInteraction<CloudType>::myLocalInteraction
(
    const myLocalInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    patchData_(pim.patchData_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    writeFields_(pim.writeFields_), 
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myLocalInteraction<CloudType>::~myLocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::volScalarField& Foam::myLocalInteraction<CloudType>::massEscape()
{
    if (!massEscapePtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massEscapePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massEscape",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    return massEscapePtr_();
}


template<class CloudType>
Foam::volScalarField& Foam::myLocalInteraction<CloudType>::massStick()
{
    if (!massStickPtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massStickPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massStick",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    return massStickPtr_();
}


template<class CloudType>
bool Foam::myLocalInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    label patchi = patchData_.applyToPatch(pp.index());

    if (patchi >= 0)
    {
        vector& U = p.U();
        vector& parcelCurl = p.parcelCurl();
	    scalar& dParcel = p.d();      
        bool& active = p.active();

        typename PatchInteractionModel<CloudType>::interactionType it =
        this->wordToInteractionType
        (
            patchData_[patchi].interactionTypeName()
        );

        switch (it)
        {
            case PatchInteractionModel<CloudType>::itNone:
            {
                return false;
            }
            case PatchInteractionModel<CloudType>::itEscape:
            {
                scalar dm = p.mass()*p.nParticle();

                keepParticle = false;
                active = false;
                U = Zero;
                nEscape_[patchi]++;
                massEscape_[patchi] += dm;
                if (writeFields_)
                {
                    label pI = pp.index();
                    label fI = pp.whichFace(p.face());
                    massEscape().boundaryFieldRef()[pI][fI] += dm;
                }
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                scalar dm = p.mass()*p.nParticle();

                keepParticle = true;
                active = false;
                U = Zero;
                nStick_[patchi]++;
                massStick_[patchi] += dm;
                if (writeFields_)
                {
                    label pI = pp.index();
                    label fI = pp.whichFace(p.face());
                    massStick().boundaryFieldRef()[pI][fI] += dm;
                }
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                active = true;

                vector nw;
                vector Up;

                this->owner().patchData(p, pp, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

		        // Betrag und Richtung der Geschwindigkeit bestimmen
		        const scalar magU = mag(U);
                /*if (magU == 0)
                {
                    cout << "magU hat den Wert 0 angenommen." << nl;
                }*/
		        const vector Udir = U/max(magU,1e-9);
		        const scalar alpha = mathematical::pi/2.0 - acos(nw & Udir);

                Random& rnd = this->owner().rndGen();

                const scalar alphaCorr = patchData_[patchi].gamma() * rnd.scalarNormal();

                const scalar alpha1 = max((alpha + alphaCorr),0);

                scalar Un = U & nw;
                //scalar UnOld = Un;
                vector Ut = U - Un*nw;
		        scalar Utabs = mag(Ut);
		        scalar e;
		        scalar mu;
		        
                scalar epsilon0 = sign(mag(Ut)-0.5*dParcel*mag(parcelCurl));
                /*
                ofstream file4("epsilon0.dat", ios::out|ios::app);
    		    file4 << epsilon0 << nl;
    		    file4.close();
                
                ofstream file("alpha.dat", ios::out|ios::app);
    		    file << alpha << nl;
    		    file.close();
                */
		        if (alpha1 == 0)
                {
                    e = patchData_[patchi].e0();
                    //cout << "Alpha hat den Wert 0 angenommen." << nl;
                }
                else if (alpha1 < patchData_[patchi].alphaE())
		        {
			        e = patchData_[patchi].e0() - ((patchData_[patchi].e0()-patchData_[patchi].eH())/patchData_[patchi].alphaE())*alpha1;
                    /*ofstream file("e.dat", ios::out|ios::app);
	                file << e << nl;
        	        file.close();*/
		        }
		        else 
		        {	
			        e = patchData_[patchi].eH();
                    /*ofstream file("e.dat", ios::out|ios::app);
                    file << e << nl;
                    file.close();*/
		        }

                // e = patchData_[patchi].e0() - (patchData_[patchi].e1() * alpha1) + (patchData_[patchi].e2() * sqr(alpha1)) - (patchData_[patchi].e3() * pow(alpha1 ,3));
                // mu = max(patchData_[patchi].muH(), patchData_[patchi].mu0() - (patchData_[patchi].mu1() * alpha1));

                if (alpha1 == 0)
                {
                    mu = patchData_[patchi].mu0();
                }
		        else if (alpha1 < patchData_[patchi].alphaMu())
		        {
			        mu = patchData_[patchi].mu0() - ((patchData_[patchi].mu0() - patchData_[patchi].muH())/patchData_[patchi].alphaMu()) * alpha1;
 			        /*ofstream file("mu.dat", ios::out|ios::app);
                    file << mu << nl;
                    file.close();*/                   
		        }
		        else
		        {
			        mu = patchData_[patchi].muH();
                    /*ofstream file("mu.dat", ios::out|ios::app);
                    file << mu << nl;
                    file.close();*/
		        }

		        if (Un > 0)
		        {
			        //U -= (1.0 + patchData_[patchi].eH())*Un*nw;
			        U -= (1.0 + e)*Un*nw;
		        }

                //U -= patchData_[patchi].muH()*(1.0 + patchData_[patchi].eH())*Un*Ut/Utabs;
		        /*if (Utabs == 0)
                {
                    cout << "Utabs hat den Wert 0 angenommen." << nl;
                }*/
                U -= mu*(1.0 + e)*epsilon0*Un*Ut/max(Utabs,1e-9);

		        parcelCurl += 5*mu*(1.0 + e)*Un*1/dParcel*epsilon0*(Ut/max(Utabs,1e-9) ^ nw);

                /*
                const scalar magUNew = mag(U);
                const vector UdirNew = U/magUNew;
                const scalar alphaNew = mathematical::pi/2.0 - acos(nw & UdirNew);
                */
		        //scalar fricCoeff = (abs((cos(alpha) * magU) - (cos(alphaNew) * magUNew)))/((1+e)*UnOld);
                U += Up;

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type "
                    << patchData_[patchi].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchData_[patchi].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::myLocalInteraction<CloudType>::info(Ostream& os)
{
    // retrieve any stored data
    labelList npe0(patchData_.size(), 0);
    this->getModelProperty("nEscape", npe0);

    scalarList mpe0(patchData_.size(), 0.0);
    this->getModelProperty("massEscape", mpe0);

    labelList nps0(patchData_.size(), 0);
    this->getModelProperty("nStick", nps0);

    scalarList mps0(patchData_.size(), 0.0);
    this->getModelProperty("massStick", mps0);

    // accumulate current data
    labelList npe(nEscape_);
    Pstream::listCombineGather(npe, plusEqOp<label>());
    npe = npe + npe0;

    scalarList mpe(massEscape_);
    Pstream::listCombineGather(mpe, plusEqOp<scalar>());
    mpe = mpe + mpe0;

    labelList nps(nStick_);
    Pstream::listCombineGather(nps, plusEqOp<label>());
    nps = nps + nps0;

    scalarList mps(massStick_);
    Pstream::listCombineGather(mps, plusEqOp<scalar>());
    mps = mps + mps0;


    forAll(patchData_, i)
    {
        os  << "    Parcel fate (number, mass)      : patch "
            <<  patchData_[i].patchName() << nl
            << "      - escape                      = " << npe[i]
            << ", " << mpe[i] << nl
            << "      - stick                       = " << nps[i]
            << ", " << mps[i] << nl;
    }

    if (this->writeTime())
    {
        this->setModelProperty("nEscape", npe);
        nEscape_ = 0;

        this->setModelProperty("massEscape", mpe);
        massEscape_ = 0.0;

        this->setModelProperty("nStick", nps);
        nStick_ = 0;

        this->setModelProperty("massStick", mps);
        massStick_ = 0.0;
    }
}


// ************************************************************************* //
