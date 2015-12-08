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
// Shape construction
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"
#include "DGtal/shapes/parametric/Flower2D.h"
// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z2i;

class PredefinedCubicalComplex {
public:
    typedef Flower2D< Space > MyEuclideanShape;
    typedef GaussDigitizer< Space, MyEuclideanShape > MyGaussDigitizer;
    typedef map<Cell, CubicalCellData>   Map;
    typedef CubicalComplex< KSpace, Map >     CC;

    KSpace K;
    CC * complex;

public:
    PredefinedCubicalComplex(Point p, double r , double smallr, unsigned int k, double phi){
      MyEuclideanShape shape = MyEuclideanShape(p, r, smallr, k, phi);
      MyGaussDigitizer digShape;
      digShape.attach( shape );
      digShape.init ( shape.getLowerBound(), shape.getUpperBound(), 1.0 );
      Domain domainShape = digShape.getDomain();
      DigitalSet aSet = DigitalSet( domainShape );
      Shapes<Domain>::digitalShaper( aSet, digShape );

      K.init ( domainShape.lowerBound(), domainShape.upperBound(), true );
      complex = new CC ( K );
      complex->construct( aSet);
    }

    CC * getCubicalFlower () const { return complex; }
//    ~PredefinedCubicalComplex() { delete  complex; }

};

///////////////////////////////////////////////////////////////////////////////

int main( int , char** )
{
  trace.beginBlock ( "Example digitalSetToCubicalComplexes2D" );
  trace.beginBlock ( "Generate a 2D shape." );
  PredefinedCubicalComplex PCC (  RealPoint( 0.0, 0.0 ), 16, 5, 5, M_PI_2/2.  );// = new PredefinedCubicalComplex();
  Board2D board;

  Color dorange ( 255,  136,  0,  220 );
  board.saveEPS("pixel-flower.eps");
  trace.endBlock();
  
  trace.beginBlock ( "Generate a 2D cubical representation." );
  typedef map<Cell, CubicalCellData>   Map;
  typedef CubicalComplex< KSpace, Map >     CC;

  
 /* typedef CC::CellMapConstIterator CellMapConstIterator;
  for ( Dimension d = 0; d <= 2; ++d )
    for ( CellMapConstIterator it = complex.begin( d ), itE = complex.end( d );
	 it != itE; ++it )
	 {
	   if ( d == 0 )
	     board << CustomStyle( it->first.className(),
				   new CustomColors( Color( 0, 0, 0 ),
						     Color( 0, 0, 0 ) ) );
	  else if ( d == 1 )
	       board << CustomStyle( it->first.className(),
				     new CustomColors( Color( 200, 0, 0 ),
						       Color( 100, 255, 100 ) ) );
	  else
		 board << CustomStyle( it->first.className(),
				       new CustomColors( Color( 0, 0, 200 ),
							 Color( 100, 255, 100 ) ) );
		 board << it->first;
	 }*/
  board << *PCC.getCubicalFlower();
  board.saveEPS ( "cubicalComplexes.eps" );
  trace.endBlock();
  trace.endBlock();
  PCC.~PredefinedCubicalComplex();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
