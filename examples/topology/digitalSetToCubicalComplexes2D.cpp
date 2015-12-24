/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file digitalSetToCubicalComplexes2D.cpp
 * @ingroup Examples
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/11/03
 *
 * An example file named digitalSetToCubicalComplexes2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <map>
#include "ConfigExamples.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
// Cellular grid
#include "DGtal/topology/CubicalComplex.h"
#include "DGtal/topology/CubicalComplexFunctions.h"
// Shape construction
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"
#include "DGtal/shapes/parametric/Ball3D.h"
// Drawing
//#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Color.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functions;
using namespace Z3i;
using  namespace ccops;


/////////////////////////////////////////////////////
typedef Ball3D< Space > MyEuclideanShape;
MyEuclideanShape shape( RealPoint( 0.0, 0.0, 0. ), 2);
typedef GaussDigitizer< Space, MyEuclideanShape > MyGaussDigitizer;
MyGaussDigitizer digShape;
typedef map<Cell, CubicalCellData> Map;
typedef CubicalComplex< KSpace, Map > CC;
typedef Viewer3D<Space, KSpace> MyViewer;
CC::DefaultCellMapIteratorPriority P;
typedef CC::CellMapConstIterator CellMapConstIterator;
///////////////////////////////////////////////////////////////////////////////

int get_direction(const Cell& F, const Cell& G, const KSpace& K)
{
    Vector V=K.uKCoords(F)-K.uKCoords(G);

    //loop over V coordinates
    for (int i=0; i<K.dimension; i++)
        if (V[i]!=0) return i;
}

int get_orientation(const Cell& F, const Cell& G, const KSpace& K)
{
    Vector V=K.uKCoords(F)-K.uKCoords(G);

    //loop over V coordinates
    for (int i=0; i<K.dimension; i++)
        if (V[i]!=0) return V[i];
}

void assign_values(const KSpace& K , CC::Iterator begin, CC& complex , int d,int orientation, int dim , std::vector<Cell> &SUB , int priority)
{
    Cells faces = K.uUpperIncident ( *begin );
    for ( int i = 0; i < faces.size(); i++ )
    {//test the cell if it's in the complex
        if ( complex.findCell ( d+1, faces[i] ) != complex.end(d+1) )
        {
            if (get_orientation(*begin, faces[i], K)==orientation &&
                    get_direction(*begin, faces[i], K)==dim &&
                    K.uDim(faces[i])==d+1)
            {//assigning value with cells
                SUB.push_back( faces[i] );
                complex.insertCell( SUB.back(), priority );
                SUB.push_back( *begin );
                complex.insertCell( SUB.back(), priority );
            }
            break;
        }
    }
}


void colorShape(CC& complex , MyViewer& board)
{
    for ( Dimension d = 0; d <= 2; ++d )
        for ( CellMapConstIterator it = complex.begin( d ), itE = complex.end( d );
              it != itE; ++it )
        {
            if ( d == 0 )
                board << CustomColors3D(Color(0, 255,0),Color(0, 255,0));//green : point
            else if ( d == 1 )
                board << CustomColors3D(Color(255, 0,0),Color(255, 0,0));//red : line
            else
                board << CustomColors3D(Color(0, 0,255 ),Color(0, 0,255));//bleu : surface
            board << it->first;
        }
}

int main( int argc, char** argv )
{
    QApplication application(argc,argv);



    trace.beginBlock ( "Example digitalSetToCubicalComplexes2D" );
    trace.beginBlock ( "Generate a 2D shape." );


    digShape.attach( shape );
    digShape.init ( shape.getLowerBound(), shape.getUpperBound(), 1.0 );
    Domain domainShape = digShape.getDomain();
    DigitalSet aSet( domainShape );
    Shapes<Domain>::digitalShaper( aSet, digShape );

    trace.endBlock();

    trace.beginBlock ( "Generate a 2D cubical representation." );

    KSpace K;//the complexe
    K.init (  Point(-20,-20,-20), Point(20,20,20), true );
    MyViewer board(K);
    board.show();
    CC complex ( K );
    complex.construct< DigitalSet >( aSet );


    std::vector<Cell> SUB;//subcomplex of k
    for (int Coll_iteration=0;Coll_iteration<1;Coll_iteration++)
    {
        CC boundary = complex.boundary();
        int priority = 0;
        for (int dim=0;dim<K.dimension;dim++)//dim===direction
        {
            for (int orientation=-1;orientation<=1;orientation+=2)
            {
                for (int d = K.dimension - 1; d >= 0; d--)//d: dimension
                {
                    for ( CC::Iterator begin = boundary.begin(); begin != boundary.end(); ++begin, ++priority )
                    {
                        if ( K.uDim (*begin ) == d )
                        {
                            assign_values(K ,begin, complex ,  d,orientation, dim , SUB , priority);
                        }
                    }
                    std::cout << "Removed " << collapse( complex, SUB.begin(), SUB.end(), P, true, true, true ) << std::endl;
                    SUB.clear();
                    priority=0;
                }
            }
        }
    }
    colorShape(complex, board);
    trace.endBlock();
    trace.endBlock();
    board<< MyViewer::updateDisplay;
    application.exec();
    return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*TODO
 set to be preserve
 std::vector<cell>w;
 w.push_back(...,cell[1],...);
 if(find(w.begin(),w.enf(),l)=v.end())  find will work for all the element of w
 {
   we to do what we want
 }
  refatorization of code
  two different algorithm
  bug report
 */


