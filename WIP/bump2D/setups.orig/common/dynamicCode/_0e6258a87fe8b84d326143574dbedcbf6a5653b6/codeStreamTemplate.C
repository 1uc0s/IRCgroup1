/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 70 "/Users/lucasertugrul/Documents/Github/IRCgroup1/WIP/bump2D/setups.orig/common/system/blockMeshDict/#codeStream"
#include "pointField.H"
        #include "mathematicalConstants.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_0e6258a87fe8b84d326143574dbedcbf6a5653b6(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 76 "/Users/lucasertugrul/Documents/Github/IRCgroup1/WIP/bump2D/setups.orig/common/system/blockMeshDict/#codeStream"
constexpr scalar xMin = 0.3;
        constexpr scalar xMax = 1.2;
        constexpr label nPoints = 100;
        constexpr scalar dx = (xMax - xMin)/scalar(nPoints - 1);
        constexpr scalar pi = constant::mathematical::pi;

        os  << "(" << nl << "spline 2 3" << nl;
        pointField profile(nPoints);

        for (label i = 0; i < nPoints; ++i)
        {
            scalar x = xMin + i*dx;
            profile[i].x() = x;
            profile[i].y() = 0.5*pow4(sin(pi*x/0.9 - (pi/3)));
//changed from 0.05 to 0.5 
            profile[i].z() = 1;
        }
        os << profile << nl;

        os << "spline 14 15" << nl;
        profile.replace(2, -1);
        os << profile << nl;

        os  << ");" << nl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

