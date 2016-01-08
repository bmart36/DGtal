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
#include <cmath>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
// Cellular grid
#include "DGtal/topology/CubicalComplex.h"
#include "DGtal/topology/CubicalComplexFunctions.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class ParDirCollapse
/**
 * Description of class 'ParDirCollapse' <p>
 * \brief Aim:
 */

template < typename CC, typename TSpace>
class ParDirCollapse
{
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Destructor.
     */
    ~ParDirCollapse() {}

    // ----------------------- Interface --------------------------------------
public:

    typedef typename CC::KSpace KSpace;
    typedef typename TSpace::Vector Vector;
    typedef typename CC::CellMapConstIterator CellMapConstIterator;
    typedef typename KSpace::Cell Cell;
    typedef typename KSpace::Cells Cells;
    typedef typename CC::Iterator Iterator;
    ParDirCollapse( const KSpace & k);
    void init ( CC * pComplex ){ complex = pComplex; }
    void exec (std::vector<Cell> &SUB, int iteration_number );

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
private:
    // ------------------------- Private Datas --------------------------------


    int getOrientation(const Cell& F, const Cell& G);
    int getDirection(const Cell& F, const Cell& G);
    void assignValues( Iterator it, int d,int orientation, int dim, std::vector<Cell> &SUB , int priority);
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

    const KSpace& K;
    CC * complex;

}; // end of class ParDirCollapse

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
