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
 * @file testGeometricalDSS.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/09/26
 *
 * @brief Functions for testing class GeometricalDSS.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/geometry/2d/GridCurve.h"

#include "DGtal/geometry/CBidirectionalSegmentComputer.h"
#include "DGtal/io/boards/CDrawableWithBoard2D.h"

#include "DGtal/geometry/2d/GeometricalDSS.h"

#include "DGtal/geometry/2d/GreedySegmentation.h"
#include "DGtal/geometry/2d/SaturatedSegmentation.h"


#include "ConfigTest.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GeometricalDSS.
///////////////////////////////////////////////////////////////////////////////
/**
 * Basic methods
 */
template <typename TCurve>
bool testGeometricalDSS(const TCurve& curve)
{

  typedef typename TCurve::IncidentPointsRange Range; //range
  typedef typename Range::ConstIterator ConstIterator; //iterator
  typedef typename Range::ConstReverseIterator ConstReverseIterator; //reverse iterator

  unsigned int nbok = 0;
  unsigned int nb = 0;
  
  trace.beginBlock ( "Constructors, copy, assignement" );
  {
    Range r = curve.getIncidentPointsRange(); //range

    GeometricalDSS<ConstIterator> s1, s2, s3;
    s2.init(r.begin()); 
    s3.init(++r.begin()); 
    GeometricalDSS<ConstIterator> s4(s2); 
    GeometricalDSS<ConstIterator> s5(s3);
    s3 = s1; 
    
    trace.info() << s1.isValid() << s1 << endl; 
    trace.info() << s2.isValid() << s2 << endl; 
    trace.info() << s3.isValid() << s3 << endl; 
    trace.info() << s4.isValid() << s4 << endl; 
    trace.info() << s5.isValid() << s5 << endl; 

    bool myFlag = (!s1.isValid())&&(!s3.isValid())
    &&(s2.isValid())&&(s4.isValid())&&(s5.isValid())
    &&(s2 == s4)&&(s3 != s5)&&(s1 == s3)&&(s2 != s5);

    nbok += myFlag ? 1 : 0; 
    nb++;
  }
  trace.endBlock();
  
  trace.beginBlock ( "Extension operations" );
  {
    Range r = curve.getIncidentPointsRange(); //range

    GeometricalDSS<ConstIterator> s, t;

    trace.info() << "forward extension " << endl; 
    ConstIterator itBegin (r.begin()); 
    ConstIterator itEnd (r.end()); 
    s.init( itBegin+1 );
    while ( (s.end() != itEnd) && (s.isExtendable()) && (s.extend()) ) {}
    trace.info() << s << endl; 
    double a, b, c; 
    s.getParameters(a,b,c); 
    trace.info() << a << " " << b << " " << c << endl; 

    t.init( (itBegin + (itEnd - itBegin)/2) ); 
    while ( (t.end() != itEnd) && (t.extend()) 
         && (t.begin() != itBegin) && (t.extendOppositeEnd()) ) {}
    trace.info() << t << endl; 

    trace.info() << "backward extension " << endl; 
    typename GeometricalDSS<ConstIterator>::Reverse rs = s.getReverse(); 
    ConstReverseIterator ritBegin (r.rbegin()); 
    ConstReverseIterator ritEnd (r.rend()); 
    rs.init( ritBegin+1 );
    while ( (rs.end() != ritEnd) && (rs.isExtendable()) && (rs.extend()) ) {}
    trace.info() << rs << endl; 
    double ap, bp, cp; 
    rs.getParameters(ap,bp,cp); 
    trace.info() << ap << " " << bp << " " << cp << endl; 

    typename GeometricalDSS<ConstIterator>::Reverse rt = t.getReverse(); 
    rt.init( (ritBegin + (ritEnd - ritBegin)/2) ); 
    while ( (rt.end() != ritEnd) && (rt.extend()) 
         && (rt.begin() != ritBegin) && (rt.extendOppositeEnd()) ) {}
    trace.info() << rt << endl; 

    trace.info() << "comparison... " << endl; 
    bool myFlag = ( (s == t)&&(rs == rt) )
    && ( s.getUf() == rs.getUf() )
    && ( s.getUl() == rs.getUl() )
    && ( s.getLf() == rs.getLf() )
    && ( s.getLl() == rs.getLl() )
    && (a == ap)
    && (b == bp)
    && (c == cp)
    ; 

    nbok += myFlag ? 1 : 0; 
    nb++;
  }
  trace.endBlock();
  
  trace.info() << "(" << nbok << "/" << nb << ") " << endl;
  return nbok == nb;
}

/*
* simple drawing
*/
template <typename TCurve>
bool drawingTestGeometricalDSS(const TCurve& curve)
{

  typedef typename TCurve::IncidentPointsRange Range; //range
  typedef typename Range::ConstIterator ConstIterator; //iterator

  Range r = curve.getIncidentPointsRange(); //range

  GeometricalDSS<ConstIterator> s;
  ConstIterator itEnd (r.end()); 
  s.init( r.begin() );
  while ( (s.end() != itEnd) && (s.extend()) ) {}

  double a, b, c; 
  s.getParameters(a,b,c); 

  Board2D board; 
  board << r << s; 
  board.saveEPS("GeometricalDSSdrawingTest.eps"); 
  return true; 
}

void testGeometricalDSSConceptChecking()
{
   typedef std::pair<PointVector<2,int>, PointVector<2,int> > Pair; 
   typedef std::vector<Pair>::const_iterator ConstIterator; 
   typedef GeometricalDSS<ConstIterator> GeomDSS; 
   BOOST_CONCEPT_ASSERT(( CDrawableWithBoard2D<GeomDSS> ));
   BOOST_CONCEPT_ASSERT(( CBidirectionalSegmentComputer<GeomDSS> ));
}

template <typename TCurve>
bool testSegmentation(const TCurve& curve)
{

  typedef typename TCurve::IncidentPointsRange Range; //range
  Range r = curve.getIncidentPointsRange(); //range
  
  typedef typename Range::ConstIterator ConstIterator; //iterator
  typedef GeometricalDSS<ConstIterator> SegmentComputer; //segment computer
  
  unsigned int nbok = 0;
  unsigned int nb = 0;  

  
  trace.beginBlock ( "Greedy segmentation" );
  {
    typedef GreedySegmentation<SegmentComputer> Segmentation;
    Segmentation theSegmentation( r.begin(), r.end(), SegmentComputer() );
    
    Board2D board; 
    board << r; 
      
    typename Segmentation::SegmentComputerIterator it = theSegmentation.begin();
    typename Segmentation::SegmentComputerIterator itEnd = theSegmentation.end();
    unsigned int c = 0; 
    for ( ; it != itEnd; ++it, ++c) {
      board << (*it); 
    }
    
    board.saveEPS("GeometricalDSSGreedySegmentationTest.eps", Board2D::BoundingBox, 5000 ); 

    nbok += (c==10) ? 1 : 0; 
    nb++;
  }
  trace.endBlock();

  trace.beginBlock ( "Saturated segmentation" );
  {
    typedef SaturatedSegmentation<SegmentComputer> Segmentation;
    Segmentation theSegmentation( r.begin(), r.end(), SegmentComputer() );
    
    Board2D board; 
    board << r; 
    
    typename Segmentation::SegmentComputerIterator it = theSegmentation.begin();
    typename Segmentation::SegmentComputerIterator itEnd = theSegmentation.end();
    unsigned int c = 0; 
    for ( ; it != itEnd; ++it, ++c) {
      cout << c << (*it); 
      board << (*it); 
    }
    
//    board.saveEPS("GeometricalDSSSaturatedSegmentationTest.eps", Board2D::BoundingBox, 5000 ); 

    nbok += (c==1) ? 1 : 0; 
    nb++;
  }
  trace.endBlock();
  
  
  trace.info() << "(" << nbok << "/" << nb << ") " << endl;
  return (nbok == nb);
}
///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class GeometricalDSS" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res; 
  
  {//concept checking
    testGeometricalDSSConceptChecking();
  }
  
  {//basic operations
    std::string filename = testPath + "samples/DSS.dat";
    ifstream instream; // input stream
    instream.open (filename.c_str(), ifstream::in);
    
    typedef KhalimskySpaceND<2,int> KSpace; 
    GridCurve<KSpace> c; //grid curve
    c.initFromVectorStream(instream);

    res = testGeometricalDSS(c)
  && drawingTestGeometricalDSS(c); 
  }
  
  {//segmentations
    std::string filename = testPath + "samples/sinus2D4.dat";
    ifstream instream; // input stream
    instream.open (filename.c_str(), ifstream::in);
    
    typedef KhalimskySpaceND<2,int> KSpace; 
    GridCurve<KSpace> c; //grid curve
    c.initFromVectorStream(instream);
    
    res = res && testSegmentation(c); 
  }
  
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////