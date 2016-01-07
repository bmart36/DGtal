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

#pragma once

/**
 * @file ParDirCollapse.h
 * @author Mohamad ONAYSSI (\c mohamad.onayssi@edu.esiee.fr )
 * ESIEE Paris
 *
 * @date 2015/12/22
 *
 * Header file for module ParDirCollapse.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ParDirCollapse_RECURSES)
#error Recursive header files inclusion detected in ParDirCollapse.h
#else // defined(ParDirCollapse_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ParDirCollapse_RECURSES

#if !defined ParDirCollapse_h
/** Prevents repeated inclusion of headers. */
#define ParDirCollapse_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
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
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Color.h"
//////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
using namespace functions;
using namespace Z3i;
using namespace ccops;

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class ParDirCollapse
/**
 * Description of class 'ParDirCollapse' <p>
 * \brief Aim:
 */
    template < typename KSpace>
    class ParDirCollapse
    {
        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Destructor.
         */
        ~ParDirCollapse();

        // ----------------------- Interface --------------------------------------
    public:





        void selfDisplay ( std::ostream & out ) const;

        /**
         * Checks the validity/consistency of the object.
         * @return 'true' if the object is valid, 'false' otherwise.
         */
        bool isValid() const;

        // ------------------------- Protected Datas ------------------------------
    private:
        // ------------------------- Private Datas --------------------------------

        typedef Ball3D< Space > MyEuclideanShape;
        typedef GaussDigitizer< Space, MyEuclideanShape > MyGaussDigitizer;
        typedef map<Cell, CubicalCellData> Map;
        typedef CubicalComplex< KSpace, Map > CC;
        typedef Viewer3D<Space, KSpace> MyViewer;
        typedef typename CC::CellMapConstIterator CellMapConstIterator;

    public:
        
        void get_orientation(const Cell& F, const Cell& G, const KSpace& K , int& shapeOrientation);
        void get_direction(const Cell& F, const Cell& G, const KSpace& K , int& shapeDirection);
        void assign_values(const KSpace& K , CC::Iterator begin, CC& complex , int d,int orientation, int dim , std::vector<Cell> &SUB , int priority);
        void collapseShape(std::vector<Cell> &SUB , const KSpace& K , CC& complex , int iteration_number );
    private:

        // ------------------------- Hidden services ------------------------------
    protected:

        /**
         * Constructor.
         * Forbidden by default (protected to avoid g++ warnings).
         */
        ParDirCollapse();

    private:

        /**
         * Copy constructor.
         * @param other the object to clone.
         * Forbidden by default.
         */
        ParDirCollapse ( const ParDirCollapse & other );

        /**
         * Assignment.
         * @param other the object to copy.
         * @return a reference on 'this'.
         * Forbidden by default.
         */
        ParDirCollapse & operator= ( const ParDirCollapse & other );

        // ------------------------- Internals ------------------------------------
    private:

    }; // end of class ParDirCollapse


/**
 * Overloads 'operator<<' for displaying objects of class 'ParDirCollapse'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'ParDirCollapse' to write.
 * @return the output stream after the writing.
 */
    std::ostream&
            operator<< ( std::ostream & out, const ParDirCollapse<KSpace> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/topology/ParDirCollapse.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ParDirCollapse_h

#undef ParDirCollapse_RECURSES
#endif // else defined(ParDirCollapse_RECURSES)