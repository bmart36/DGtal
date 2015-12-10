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
#include "DGtal/shapes/parametric/Flower2D.h"
// Drawing
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functions;
using namespace Z2i;
using  namespace ccops;

///////////////////////////////////////////////////////////////////////////////

int direction(const Cell& F, const Cell& G, const KSpace& K)
{
    Vector V=K.uKCoords(F)-K.uKCoords(G);

    //loop over V coordinates
    for (int i=0; i<K.dimension; i++)
        if (V[i]!=0) return i;
}

int orientation(const Cell& F, const Cell& G, const KSpace& K)
{
    Vector V=K.uKCoords(F)-K.uKCoords(G);

    //loop over V coordinates
    for (int i=0; i<K.dimension; i++)
        if (V[i]!=0) return V[i];
}

int main( int , char** )
{
  trace.beginBlock ( "Example digitalSetToCubicalComplexes2D" );
  trace.beginBlock ( "Generate a 2D shape." );
  typedef Flower2D< Space > MyEuclideanShape;
  MyEuclideanShape shape( RealPoint( 0.0, 0.0 ), 16, 5, 5, M_PI_2/2. );

  typedef GaussDigitizer< Space, MyEuclideanShape > MyGaussDigitizer;
  MyGaussDigitizer digShape;
  digShape.attach( shape );
  digShape.init ( shape.getLowerBound(), shape.getUpperBound(), 1.0 );
  Domain domainShape = digShape.getDomain();
  DigitalSet aSet( domainShape );
  Shapes<Domain>::digitalShaper( aSet, digShape );

  Board2D board;
  board << SetMode( domainShape.className(), "Paving" ) << domainShape;
  Color dorange ( 255,  136,  0,  220 );
  board << CustomStyle( aSet.className(), new CustomFillColor ( dorange ) );
  board << aSet;
  trace.endBlock();

  trace.beginBlock ( "Generate a 2D cubical representation." );
  typedef map<Cell, CubicalCellData> Map;
  typedef CubicalComplex< KSpace, Map > CC;

  KSpace K;
  K.init (  domainShape.lowerBound(), domainShape.upperBound(), true );
  CC complex ( K );
  complex.construct< DigitalSet >( aSet );

  board << SetMode( domainShape.className(), "Paving" ) << domainShape;

  typedef CC::CellMapConstIterator CellMapConstIterator;
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
     }

  board.clear();
  board.saveEPS ( "cubicalComplexes.eps" );

  CC::DefaultCellMapIteratorPriority P;
  std::vector<Cell> SUB;
  CC boundaryBefore = complex.boundary();
  for (int q=0;q<3;q++)
  {
      CC boundary = complex.boundary();
      int prio = 0;
      for (int t=0;t<K.dimension;t++)
      {
          for (int s=-1;s<=1;s+=2)
          {
              for (int d = 1; d >= 0; d--)
              {
                  for ( CC::Iterator begin = boundary.begin(); begin != boundary.end(); ++begin, ++prio )
                  {
                      if ( K.uDim (*begin ) == d )
                      {
                          Cells faces = K.uUpperIncident ( *begin );
                          for ( int i = 0; i < faces.size(); i++ )
                          {
                              if ( complex.findCell ( d+1, faces[i] ) != complex.end(d+1) )
                              {
                                  if (/**/orientation(*begin, faces[i], K)==s &&
                                          direction(*begin, faces[i], K)==t &&
                                          K.uDim(faces[i])==d+1)
                                  {
                                      SUB.push_back( faces[i] );
                                      complex.insertCell( SUB.back(), prio );
                                      SUB.push_back( *begin );
                                      complex.insertCell( SUB.back(), prio );
                                  }
                                  break;
                              }
                          }
                      }
                  }

                  collapse( complex, SUB.begin(), SUB.end(), P, true, true, true );

                  SUB.clear();
                  prio=0;
              }
          }
      }
      board << CustomStyle( complex.className(), new CustomColors( Color::Magenta, Color::Magenta ) )
            << boundary;
  }

//    board << complex;

//  board << CustomStyle( complex.className(), new CustomColors( Color::Magenta, Color::Magenta ) )
//        << boundaryBefore;

  board.saveEPS ( "cubicalCollapsedComplexes.eps" );
  trace.endBlock();
  trace.endBlock();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
