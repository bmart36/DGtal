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


template <typename TKSpace, typename TCellContainer,
          typename CellConstIterator,
          typename CellMapIteratorPriority >
DGtal::uint64_t collapse( CubicalComplex< TKSpace, TCellContainer >& K,
                            CellConstIterator S_itB, CellConstIterator S_itE,
                            const CellMapIteratorPriority& priority,
                            bool hintIsSClosed, bool hintIsKClosed,
                            bool verbose )
{

  using namespace std;
  typedef CubicalComplex< TKSpace, TCellContainer > CC;
  typedef typename CC::Cell                         Cell;
  typedef typename CC::CellType                     CellType;
  typedef typename CC::CellMapIterator              CellMapIterator;
  typedef vector< CellMapIterator >                 CMIVector;
  typedef typename CMIVector::const_iterator        CMIVectorConstIterator;
  // NB : a maximal k-cell is collapsible if it has a free incident k-1-cell.
  Dimension n = K.dim();
  CMIVector S;            // stores the cells to process
  CMIVector Q_low;        // stores the iterators on direct faces of the maximal cell.
  CMIVector Q_collapsible;// stores collapsible cells in order to clean them at the end.
  CellMapIterator        it_cell;   // generic iterator on a cell.
  CellMapIterator        it_cell_c; // points on c in the free pair (c,d)
  CellMapIterator        it_cell_d; // points on d in the free pair (c,d)
  CMIVectorConstIterator itlow;
  CMIVectorConstIterator itqup;

  if ( verbose ) trace.info() << "[CC::collapse]-+ tag collapsible elements... " << flush;
  // Restricts the set of elements that are collapsible.
  if ( hintIsSClosed )
    for ( CellConstIterator S_it = S_itB; S_it != S_itE; ++S_it )
      {
        Cell c           = *S_it;
        Dimension k      = K.dim( c );
        it_cell          = K.findCell( k, c );
    uint32_t& ccdata = it_cell->second.data;
        ASSERT( it_cell != K.end( k ) );
        S.push_back( it_cell );
    if ( ! ( ccdata & (CC::FIXED | CC::COLLAPSIBLE ) ) )
      {
        ccdata |= CC::COLLAPSIBLE;
        Q_collapsible.push_back( it_cell );
      }
      }
  else // not ( hintIsSClosed )
    for ( CellConstIterator S_it = S_itB; S_it != S_itE; ++S_it )
      {
        Cell c           = *S_it;
        Dimension k      = K.dim( c );
        it_cell          = K.findCell( k, c );
        uint32_t& ccdata = it_cell->second.data;
        ASSERT( it_cell != K.end( k ) );
        S.push_back( it_cell );
        if ( ! ( ccdata & (CC::FIXED | CC::COLLAPSIBLE ) ) )
          {
            ccdata |= CC::COLLAPSIBLE;
            Q_collapsible.push_back( it_cell );
          }
        vector<Cell> cells;
        back_insert_iterator< vector<Cell> > back_it( cells );
        K.faces( back_it, c, hintIsKClosed );
        for ( typename vector<Cell>::const_iterator
                it = cells.begin(), itE = cells.end(); it != itE; ++it )
          {
            it_cell           = K.findCell( *it );
            uint32_t& ccdata2 = it_cell->second.data;
            if ( ! ( ccdata2 & (CC::FIXED | CC::COLLAPSIBLE ) ) )
              {
                ccdata2 |= CC::COLLAPSIBLE;
                Q_collapsible.push_back( it_cell );
              }
          }
      }
  if ( verbose ) trace.info() << " " << Q_collapsible.size() << " found." << endl;

  // Fill queue
  priority_queue<CellMapIterator, CMIVector, CellMapIteratorPriority> PQ( priority );

  if ( verbose ) trace.info() << "[CC::collapse]-+ entering collapsing loop. " << endl;
  uint64_t nb_pass     = 0;
  uint64_t nb_examined = 0;

    uint64_t nb_removed  = 0;

  while ( ! S.empty() )
    {
      for ( CMIVectorConstIterator it = S.begin(), itE = S.end();
            it != itE; ++it )
    {
      PQ.push( *it );
      (*it)->second.data |= CC::USER1;
    }
      S.clear();
      if ( verbose ) trace.info() << "[CC::collapse]---+ Pass " << ++nb_pass
                                  << ", Card(PQ)=" << PQ.size() << " elements, "
                                  << "nb_exam=" << nb_examined << endl;

      // Try to collapse elements according to priority queue.
      while ( ! PQ.empty() )
    {
      // Get top element.
      CellMapIterator itcur = PQ.top();
          uint32_t& cur_data    = itcur->second.data;
      PQ.pop();
      ++nb_examined;
      // Check if the cell is removable
      if ( ( cur_data & CC::REMOVED ) || ( ! ( cur_data & CC::COLLAPSIBLE ) ) )
            continue;
          // Check if the cell was several time in the queue and is already processed.
          if ( ! ( cur_data & CC::USER1 ) )
            continue;
      ASSERT( cur_data & CC::USER1 );
          cur_data             &= ~CC::USER1;

      // Cell may be removable.
      // Check if it is a maximal cell
      CellMapIterator itup;
      const Cell & cur_c  = itcur->first;
          CellType cur_c_type = K.computeCellType( cur_c, itup, n );
      bool found_pair     = false;
          // trace.info() << "  - Cell " << cur_c << " Dim=" << dim( cur_c ) << " Type=" << cur_c_type << std::endl;
      if ( cur_c_type == CC::Maximal )
        { // maximal cell... must find a free face
          // check faces to find a free face.
              back_insert_iterator< CMIVector > back_it( Q_low );
              K.directFacesIterators( back_it, cur_c );
          CellMapIterator best_free_face_it = K.end( 0 );
          for ( CMIVectorConstIterator it = Q_low.begin(), itE = Q_low.end();
                    it != itE; ++it )
        {
                  CellMapIterator low_ic = *it;
                  uint32_t& data         = low_ic->second.data;
                  // trace.info() << "    + Cell " << low_ic->first << " data=" << data << std::endl;
          if ( ( data & CC::REMOVED ) || ! ( data & CC::COLLAPSIBLE ) ) continue;
          const Cell& cur_d   = low_ic->first;
                  CellType cur_d_type = K.computeCellType( cur_d, itup, n );
                  // trace.info() << "      + Type=" << cur_d_type << std::endl;
          if ( cur_d_type == CC::Free )
            { // found a free n-1-face ic
              if ( ( best_free_face_it == K.end( 0 ) )
               || ( ! priority( low_ic, best_free_face_it ) ) )
            best_free_face_it = low_ic;
            }
        }
          if ( best_free_face_it != K.end( 0 ) )
        {
          // delete c and ic.
          found_pair = true;
          it_cell_c  = itcur;
          it_cell_d  = best_free_face_it;
          // Q_low already contains cells that should be
          // checked again
        }
        }
      else if ( cur_c_type == CC::Free )
        { // free face... check that its 1-up-incident face is maximal.
          CellMapIterator it_up_up;
          const Cell& cur_d   = itup->first;
              CellType cur_d_type = K.computeCellType( cur_d, it_up_up, n );
          if ( cur_d_type == CC::Maximal )
        { // found a maximal face.
          found_pair = true;
          it_cell_c  = itup;
          it_cell_d  = itcur;
          // Q_low will contain cells that should be checked
          // again
                  back_insert_iterator< CMIVector > back_it( Q_low );
                  K.directFacesIterators( back_it, it_cell_c->first );
        }
        }
      if ( found_pair )
        { // If found, remove pair from complex (logical removal).
          it_cell_c->second.data |= CC::REMOVED;
          it_cell_d->second.data |= CC::REMOVED;
              nb_removed             += 2;
          // Incident cells have to be checked again.
          for ( CMIVectorConstIterator it = Q_low.begin(), itE = Q_low.end();
                    it != itE; ++it )
                {
                  it_cell             = *it;
          uint32_t& data_qlow = it_cell->second.data;
          if ( ( ! ( data_qlow & CC::REMOVED ) )
               && ( data_qlow & CC::COLLAPSIBLE )
               && ( ! ( data_qlow & CC::USER1 ) ) )
            {
              S.push_back( it_cell );
            }
        }
        }
          Q_low.clear();
    } // while ( ! PQ.empty() )
    } // while ( ! S.empty() )

  if ( verbose ) trace.info() << "[CC::collapse]-+ cleaning complex." << std::endl;

  // Now clean the complex so that removed cells are effectively
  // removed and no more cell is tagged as collapsible.
  for ( CMIVectorConstIterator it = Q_collapsible.begin(), itE = Q_collapsible.end();
        it != itE; ++it )
    {
      CellMapIterator cmIt  = *it;
      uint32_t& cur_data    = cmIt->second.data;
      if ( cur_data & CC::REMOVED ) K.eraseCell( cmIt );
      else                          cur_data &= ~CC::COLLAPSIBLE;
      // if ( (*it)->second.data & CC::REMOVED )
      //   K.eraseCell( *it );
      // else
      //   (*it)->second.data &= ~CC::COLLAPSIBLE;
    }

  return nb_removed;
}




///////////////////////////////////////////////////////////////////////////////

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
  typedef map<Cell, CubicalCellData>   Map;
  typedef CubicalComplex< KSpace, Map >     CC;  
  
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
  
//  // TO PART MOU
    CC::DefaultCellMapIteratorPriority P;

    Integer m=40;

    std::vector<Cell> SUB;

//    for ( Integer x = 0; x <= m; ++x )
//      for ( Integer y = 0; y <= m; ++y )
//        for ( Integer z = 0; z <= m; ++z )
//          {
//           //Point k1 = Point( x, y, z );
//            Point k1 = Point( x, y );

//            SUB.push_back( K.uCell( k1 ) );

//            double d1 = Point::diagonal( 1 ).dot( k1 ) / (double) KSpace::dimension; // sqrt( (double) KSpace::dimension );

//            //RealPoint v1( k1[ 0 ], k1[ 1 ], k1[ 2 ] );

//            RealPoint v1( k1[ 0 ], k1[ 1 ]);

//            v1 -= d1 * RealPoint::diagonal( 1.0 );

//            double n1 = v1.norm();
//            bool fixed = ( ( x == 0 ) && ( y == 0 ) )
//              || ( ( x == 0 ) && ( y == m ) )
//              || ( ( x == m ) && ( y == 0 ) )
//              || ( ( x == m ) && ( y == m ) )
//              || ( ( x == m/3 ) && ( y == 2*m/3 ) )
//              || ( ( x == 0 ) && ( y == 0 ) )
//              || ( ( x == 0 ) && ( y == m ) )
//              || ( ( x == m ) && ( y == 0 ) )
//              || ( ( x == m ) && ( y == m ) )
//              || ( ( x == 0 ) && ( y == m ) )
//              || ( ( x == m ) && ( y == m ) )
//              || ( ( z == 0 ) && ( y == m ) )
//              || ( ( z == m ) && ( y == m ) );
//            complex.insertCell( SUB.back(), (uint32_t) floor(64.0 * n1 ) );

//          }

    complex.fillData ( 0, CC::FIXED );
    complex.fillData ( 1, CC::FIXED );
    complex.fillData ( 2, CC::FIXED );

    std::cout << "test " << std::endl;
    board.clear();

    CC::Iterator begin = complex.boundary().begin();

    for (; begin != complex.boundary().end(); ++begin )
    {
        if ( K.uDim (*begin ) == 1 )
        {
            board << *begin;
            complex.insertCell( *begin, 10 );
            SUB.push_back( *begin );
            Cells faces = K.uUpperIncident ( *begin );

            for (int i = 0; i < faces.size(); i++)
            {
                if ( complex.find ( faces[i] ) != complex.end() )
                {
                    SUB.push_back( faces[i] );
                    complex.insertCell( faces[i], 10 );
                    board << faces[i];
                    break;
                }
            }

            Cells facesl = K.uLowerIncident ( *begin );

            for (int i = 0; i < facesl.size(); i++)
            {
                if ( complex.find ( facesl[i] ) != complex.end() )
                {
                    SUB.push_back( facesl[i] );
                    complex.insertCell( facesl[i], 11 );
                    board << facesl[i];
//                    break;
                }
            }

        }
    }

    std::cout << "collapse " << std::endl;

    uint64_t removed = collapse( complex, SUB.begin(), SUB.end(), P, true, true, true );

    std::cout << "Number " << removed << std::endl;

//    board << complex;

  board.saveEPS ( "cubicalComplexes.eps" );
  trace.endBlock();
  trace.endBlock();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
