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
#include "DGtal/io/viewers/Viewer3D.h"
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
#include "DGtal/io/Color.h"

#include "DGtal/topology/ParDirCollapse.h"
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

//    aSet.insert(*it);

    trace.endBlock();
    trace.beginBlock ( "Generate a 2D cubical representation." );
    KSpace K;
    K.init (  Point(-20,-20,-20), Point(20,20,20), true );
    MyViewer board(K);
    board.show();
    CC complex ( K );
    complex.construct< DigitalSet >( aSet );
    ParDirCollapse<CC, Space> dirCollape( K );
    dirCollape.init ( &complex );
    //dirCollape.exec ( 1 );
    dirCollape.collapseSurface();
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
 if(find(w.begin(),w.end(),l)=SUB.end())  find will work for all the element of w
 {
   we can do what we want
 }
  refatorization of code
  two different algorithm
  bug report
 */
