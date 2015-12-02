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
* @file CubicalComplex.h
* @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
* Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
*
* @date 2015/08/28
*
* Header file for module CubicalComplex.cpp
*
* This file is part of the DGtal library.
*/

#if defined(CubicalComplex_RECURSES)
#error Recursive header files inclusion detected in CubicalComplex.h
#else // defined(CubicalComplex_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CubicalComplex_RECURSES

#if !defined CubicalComplex_h
/** Prevents repeated inclusion of headers. */
#define CubicalComplex_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/type_traits.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/ContainerTraits.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /**
  * Any cell is stored within a cubical complex with an associated
  * data, which must derive from this class. Its basic usage is to
  * store flags associated to the cells, but it may store other
  * values.
  *
  * @note Predefined flags are CubicalComplex::REMOVED,
  * CubicalComplex::COLLAPSIBLE, CubicalComplex::FIXED,
  * CubicalComplex::USER1. Other bits can be used to associate an
  * integer to the cell. The corresponding mask is
  * CubicalComplex::VALUE. 
  *
  * @note Such data is notably used in collapse operation
  * (CubicalComplex::collapse).
  */
  struct CubicalCellData {
    inline CubicalCellData() : data( 0 ) {}
    CubicalCellData( uint32_t d ) : data( d ) {}
    uint32_t data;
  };


  /////////////////////////////////////////////////////////////////////////////
  // template class CubicalComplex
  /**
  * Description of template class 'CubicalComplex' <p> \brief Aim:
  * This class represents an arbitrary cubical complex living in some
  * Khalimsky space. Cubical complexes are sets of cells of different
  * dimensions related together with incidence relations. Two cells
  * in a cubical complex are incident if and only if they are
  * incident in the surrounding Khalimsky space. In other words,
  * cubical complexes are defined here as subsets of Khalimsky spaces. 
  *
  * @tparam TKSpace any model of concepts::CCellularGridSpaceND, i.e. a type
  * that models a Khalimsky space.
  *
  * @tparam TCellContainer any model of associative container, mapping
  * a KSpace::Cell to a CubicalCellData or any type deriving from
  * it. It could be for instance a std::map or a
  * std::unordered_map. Note that unfortunately, unordered_map are
  * (strangely) not models of boost::AssociativeContainer, hence we
  * cannot check concepts here.
  *
  * @tparam TData any type deriving from CubicalCellData that is
  * boost::DefaultConstructible, boost::Assignable,
  * boost::CopyConstructible.
  */
  template < typename TKSpace, 
             typename TCellContainer = std::map< typename TKSpace::Cell, CubicalCellData > >
  class CubicalComplex
  {
    // ----------------------- associated types ------------------------------
  public:
    
    BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND< TKSpace > ));
    // BOOST_CONCEPT_ASSERT(( boost::AssociativeContainer< TCellContainer > ));
    // BOOST_CONCEPT_ASSERT(( boost::PairAssociativeContainer< TCellContainer > ));

    typedef TKSpace KSpace;
    typedef TCellContainer CellContainer;
    typedef typename CellContainer::mapped_type Data;

    BOOST_STATIC_ASSERT (( boost::is_base_of< CubicalCellData, Data >::value ));
    BOOST_STATIC_ASSERT (( boost::is_same< typename TKSpace::Cell, typename CellContainer::key_type >::value ));

    /// The dimension of the embedding space.
    static const Dimension dimension = KSpace::dimension;
    typedef typename KSpace::Integer     Integer;
    typedef typename KSpace::Cell        Cell;
    typedef typename KSpace::Cells       Cells;
    typedef typename KSpace::Space       Space;
    typedef typename KSpace::Size        Size;
    typedef typename KSpace::Point       Point;
    typedef typename KSpace::DirIterator DirIterator;
    typedef CellContainer                CellMap;
    typedef typename CellMap::const_iterator CellMapConstIterator;
    typedef typename CellMap::iterator   CellMapIterator;

    /// Flag Used to indicate in a cell data that this cell has been (virtually) removed.
    static const uint32_t REMOVED     = 0x10000000;
    /// Flag Used to indicate in a cell data that this cell is collapsible.
    static const uint32_t COLLAPSIBLE = 0x20000000;
    /// Flag Used to indicate in a cell data that this cell is fixed.
    static const uint32_t FIXED       = 0x40000000;
    /// User flag for a cell.
    static const uint32_t USER1       = 0x80000000;
    /// Value for a cell.
    static const uint32_t VALUE       = 0x0fffffff;

    // ----------------------- inner types ------------------------------------

    /**
    * This structure is used to order cells (represented by their
    * iterator in the CubicalComplex container) in a priority queue.
    * By default, it compares the values associated to each cell (the
    * value is the cell data masked by VALUE).
    */
    struct DefaultCellMapIteratorPriority {

      /**
      * This operator compares two cells specified by their iterators
      * \a it1 and \a it2 and returns true whenever \a it1 has
      * smallest data & VALUE than \a it2.
      *
      * @param it1 any iterator on a cell of this cubical complex.
      * @param it2 any iterator on a cell of this cubical complex.
      * @return 'true' iff the value part of the data of \a it1 is smaller than the one of \a it2.
      */
      bool operator()( const CellMapIterator& it1, const CellMapIterator& it2 ) const
      {
        uint32_t v1 = it1->second.data & VALUE; 
        uint32_t v2 = it2->second.data & VALUE; 
        return ( v1 < v2 ) 
          || ( ( v1 == v2 ) && ( it1->first < it2->first ) );
      }
    };
    // ----------------------- Standard services ------------------------------
  public:

    /**
    * Destructor.
    */
    ~CubicalComplex();

    /**
    * Constructor of empty complex. Needs a space to represents
    * cubical cells.
    *
    * @param aK a Khalimsky space.
    */
    CubicalComplex ( ConstAlias<KSpace> aK );

    /**
    * Copy constructor.
    * @param other the object to clone.
    * Forbidden by default.
    */
    CubicalComplex ( const CubicalComplex & other );
    
    /**
    * Constructor a complex from a digital set.
    * @param set - a digital set from which to create a complex. 
    * Set has to be of the same dimension as a Khalimsky space.
    */
    template < typename TDigitalSet >
    void construct ( const TDigitalSet & set );

    /**
    * Assignment.
    * @param other the object to copy.
    * @return a reference on 'this'.
    * Forbidden by default.
    */
    CubicalComplex & operator= ( const CubicalComplex & other );

    /**
    * Clears the cubical complex, which becomes empty.
    */
    void clear();

    /**
    * Clears all cells of dimension \a d of the cubical complex.
    * @param d the dimension of cell \a aCell.
    */
    void clear( Dimension d );

    /**
    * Fills the data of every cell of this cubical complex, which
    * becomes \a data. Default value resets the data to zero.
    *
    * @param data any data.
    */
    void fillData( Data data = Data() );

    /**
    * Fills the data of every cell of dimension \a d this cubical
    * complex, which becomes \a data. Default value resets the data to
    * zero.
    *
    * @param d the dimension of cell \a aCell.
    * @param data any data.
    */
    void fillData( Dimension d, Data data = Data() );

    /**
    * @return the maximal dimension of a cell in the complex, 0 if
    * the complex is empty.
    */
    Dimension dim() const;

    /**
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return the dimension of the cell \a aCell.
    */
    Dimension dim( const Cell& aCell ) const;

    /**
    * @param d the dimension of cells.
    * @return the number of cells of dimension \a d in this complex.
    */
    Size nbCells( Dimension d ) const;

    /**
    * @note For instance, all Platonician solids have euler number
    * equal to one, while their surface have euler number equal to
    * two.
    *
    * @return the Euler number of this complex which equals nbCells( 0 ) - nbCells( 1 ) + nbCells( 2 ) - ...
    */
    Integer euler() const;

    /**
    * @return a reference to the Khalimsky space associated to this complex.
    */
    const KSpace& space() const;

    /**
    * Insert cell \a aCell into CubicalComplex and assign to it the value \a data.
    *
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @param data any value.
    */
    void insertCell( const Cell& aCell, const Data& data = Data() );

    /**
    * Insert cell \a aCell into CubicalComplex and assign to it the
    * value \a data. Faster than the other insertCell method.
    *
    * @param d the dimension of cell \a aCell.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @param data any value.
    */
    void insertCell( Dimension d, const Cell& aCell, const Data& data = Data() );

    /**
    * Insert the cells within range [it,itE) into the
    * CubicalComplex. The value associated to each cell is the
    * default.
    *
    * @param it an iterator pointing at the beginning of a range of (arbitrary) cells.
    * @param itE an iterator pointing after the end of a range of (arbitrary) cells.
    * @param data any value.
    * @tparam CellConstIterator any model of a forward const iterator on Cell.
    */
    template <typename CellConstIterator>
    void insertCells( CellConstIterator it, CellConstIterator itE, const Data& data = Data() );

    /**
    * Insert the cells within range [it,itE) into the
    * CubicalComplex. The value associated to each cell is the
    * default.
    *
    * @param d the dimension of all cells in the range [it,itE).
    * @param it an iterator pointing at the beginning of a range of (arbitrary) cells.
    * @param itE an iterator pointing after the end of a range of (arbitrary) cells.
    * @param data any value.
    * @tparam CellConstIterator any model of a forward const iterator on Cell.
    */
    template <typename CellConstIterator>
    void insertCells( Dimension d, CellConstIterator it, CellConstIterator itE, const Data& data = Data() );
    
    /**
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return 'true' if and only if \a aCell belongs to this complex.
    */
    bool belongs( const Cell& aCell ) const;

    /**
    * @param d the dimension of cell \a aCell.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return 'true' if and only if \a aCell belongs to this complex.
    */
    bool belongs( Dimension d, const Cell& aCell ) const;

    /**
    * Erases cell \a aCell from the complex.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return the number of cells effectively removed from the cubical complex.
    */
    Size eraseCell( const Cell& aCell );

    /**
    * Erases cell \a aCell from the complex.
    * @param d the dimension of cell \a aCell.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return the number of cells effectively removed from the cubical complex.
    */
    Size eraseCell( Dimension d, const Cell& aCell );

    /**
    * Erases cell pointed by iterator \a it.
    * @param it any valid iterator on a cell of this complex.
    */
    void eraseCell( CellMapIterator it );

    /**
    * Erases cells in range [\a it, \a itE).
    * @param it any valid iterator on the first element of range of cells of this complex.
    * @param itE any valid iterator after the last element of range of cells of this complex.
    */
    void eraseCells( CellMapIterator it, CellMapIterator itE );

    /**
    * Erases the cells stored in range [it,itE) from the
    * CubicalComplex. 
    *
    * @param it an iterator pointing at the beginning of a range of (arbitrary) cells.
    * @param itE an iterator pointing after the end of a range of (arbitrary) cells.
    * @return the number of cells effectively removed from the cubical complex.
    * @tparam CellConstIterator any model of a forward const iterator on Cell.
    */
    template <typename CellConstIterator>
    Size eraseCells( CellConstIterator it, CellConstIterator itE );

    /**
    * Erases the cells of dimension \a d stored in range [it,itE) from the
    * CubicalComplex. 
    *
    * @param d the dimension of every cell in range [it,itE).
    * @param it an iterator pointing at the beginning of a range of (arbitrary) cells.
    * @param itE an iterator pointing after the end of a range of (arbitrary) cells.
    * @return the number of cells effectively removed from the cubical complex.
    * @tparam CellConstIterator any model of a forward const iterator on Cell.
    */
    template <typename CellConstIterator>
    Size eraseCells( Dimension d, CellConstIterator it, CellConstIterator itE );

    /**
    * Makes CubicalComplex a functor Cell -> boolean, which represents
    * the characteristic cell function.
    *
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return 'true' if and only if \a aCell belongs to this complex.
    */
    bool operator()( const Cell& aCell ) const;

    /**
    * Outputs all the cells that are proper faces of \a aCell with output iterator \a it.
    *
    * @param outIt the output iterator on Cell that is used for outputing faces.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintClosed when 'true', this hint tells that the complex
    * is closed, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @tparam CellOutputIterator any model of boost::OutputIterator, with value_type Cell.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellOutputIterator>
    void faces( CellOutputIterator& outIt, const Cell& aCell, 
                bool hintClosed = false ) const;

    /**
    * Outputs all the cells that are direct faces of \a aCell with
    * output iterator \a it (direct faces are lower incident cells
    * with a dimension just one below).
    *
    * @param outIt the output iterator on Cell that is used for outputing faces.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintClosed when 'true', this hint tells that the complex
    * is closed, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @tparam CellOutputIterator any model of boost::OutputIterator, with value_type Cell.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellOutputIterator>
    void directFaces( CellOutputIterator& outIt, const Cell& aCell,
                      bool hintClosed = false ) const;

    /**
    * Outputs all the iterators on cells that are direct faces of \a aCell with
    * output iterator \a it (direct faces are lower incident cells
    * with a dimension just one below).
    *
    * @param outIt the output iterator on CellMapIterator that is used for outputing face iterators.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    *
    * @tparam CellMapIteratorOutputIterator any model of boost::OutputIterator, with value_type CellMapIterator.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellMapIteratorOutputIterator>
    void directFacesIterators( CellMapIteratorOutputIterator& outIt, const Cell& aCell );

    /**
    * Outputs all the cells that are proper co-faces of \a aCell with
    * output iterator \a it.
    *
    * @param outIt the output iterator on Cell that is used for outputing faces.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintOpen when 'true', this hint tells that the complex
    * is open, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @tparam CellOutputIterator any model of boost::OutputIterator, with value_type Cell.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellOutputIterator>
    void coFaces( CellOutputIterator& outIt, const Cell& aCell,
                  bool hintOpen = false ) const;

    /**
    * Outputs all the cells that are direct co-faces of \a aCell with
    * output iterator \a it (direct faces are upper incident cells
    * with a dimension just one above).
    *
    * @param outIt the output iterator on Cell that is used for outputing faces.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintOpen when 'true', this hint tells that the complex
    * is open, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @tparam CellOutputIterator any model of boost::OutputIterator, with value_type Cell.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellOutputIterator>
    void directCoFaces( CellOutputIterator& outIt, const Cell& aCell,
                        bool hintOpen = false ) const;

    /**
    * Outputs all the iterators on cells that are direct co-faces of \a aCell with
    * output iterator \a it (direct faces are upper incident cells
    * with a dimension just one above).
    *
    * @param outIt the output iterator on CellMapIterator that is used for outputing face iterators.
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    *
    * @tparam CellMapIteratorOutputIterator any model of boost::OutputIterator, with value_type CellMapIterator.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    template <typename CellMapIteratorOutputIterator>
    void directCoFacesIterators( CellMapIteratorOutputIterator& outIt, const Cell& aCell );

    /**
    * @param d any valid dimension.
    * @return a const iterator pointing on the first cell of dimension \a d of this.
    */
    CellMapConstIterator begin( Dimension d ) const;

    /**
    * @param d any valid dimension.
    * @return a const iterator pointing after the last cell of dimension \a d of this.
    */
    CellMapConstIterator end( Dimension d ) const;

    /**
    * @param d any valid dimension.
    * @return an iterator pointing on the first cell of dimension \a d of this.
    */
    CellMapIterator begin( Dimension d );

    /**
    * @param d any valid dimension.
    * @return an iterator pointing after the last cell of dimension \a d of this.
    */
    CellMapIterator end( Dimension d );

    /**
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return an iterator pointing on the pair (aCell,data) if the cell belongs to the complex, or end( dim( aCell ) ) 
    */
    CellMapConstIterator find( const Cell& aCell ) const;

    /**
    * @param d the dimension of cell \a aCell.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return an iterator pointing on the pair (aCell,data) if the cell belongs to the complex, or end( d ) 
    */
    CellMapConstIterator find( Dimension d, const Cell& aCell ) const;

    /**
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return an iterator pointing on the pair (aCell,data) if the cell belongs to the complex, or end( dim( aCell ) ) 
    */
    CellMapIterator find( const Cell& aCell );

    /**
    * @param d the dimension of cell \a aCell.
    * @param aCell any cell valid in the Khalimsky space associated to the complex.
    * @return an iterator pointing on the pair (aCell,data) if the cell belongs to the complex, or end( d ) 
    */
    CellMapIterator find( Dimension d, const Cell& aCell );

    // ---------- local operations for extracting specific subcomplexes ---------------
  public:

    /**
    * Returns the boundary of the cell \a aCell as a cell collection,
    * i.e. all the cells that are proper faces of \a aCell. Generally
    * faster than method \ref faces, which outputs cells with an
    * output iterator.
    *
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintClosed when 'true', this hint tells that the complex
    * is (locally) closed, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @return the collection of cells that defines the boundary of \a
    * aCell, i.e. its proper faces.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    Cells cellBoundary( const Cell& aCell, bool hintClosed = false ) const;

    /**
    * Returns the co-boundary of the cell \a aCell as a cell collection,
    * i.e. all the cells that are proper co-faces of \a aCell. Generally
    * faster than method \ref coFaces, which outputs cells with an
    * output iterator.
    *
    * @param aCell any cell valid in the Khalimsky space associated to the complex. 
    * @param hintOpen when 'true', this hint tells that the complex
    * is (locally) open, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    *
    * @return the collection of cells that defines the co-boundary of \a aCell, i.e. its proper co-faces.
    *
    * @note all returned cells belong to this complex, while it is
    * not compulsory for \a aCell to belong to it.
    */
    Cells cellCoBoundary( const Cell& aCell, bool hintOpen = false ) const;


    // ----------------------- Standard subcomplexes --------------------------------------
  public:
    
    /**
    * Returns the closure of the cells in \a S within this complex,
    * i.e. the smallest subcomplex that contains each cell in \a S.
    *
    * @tparam CellAssociativeContainer any simple associative container with key Cell.
    * @param S any collection of cells of this complex.
    * @param hintClosed when 'true', this hint tells that the complex
    * is (locally around \a S) closed, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    * @return the closure of \a S within this complex.
    */
    template <typename CellAssociativeContainer>
    CellAssociativeContainer closure( const CellAssociativeContainer& S, bool hintClosed = false ) const;

    /**
    * Returns the star of the cells in \a S within this complex, i.e. the
    * set of all cells of this complex that have any faces in \a S.
    *
    * @tparam CellAssociativeContainer any simple associative container with key Cell.
    * @param S any collection of cells of this complex.
    * @param hintOpen when 'true', this hint tells that the complex
    * is (locally around \a S) open, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    * @return the star of \a S within this complex.
    */
    template <typename CellAssociativeContainer>
    CellAssociativeContainer star( const CellAssociativeContainer& S, bool hintOpen = false ) const;

    /**
    * Returns the link of the cells in \a S within this complex,
    * i.e. the closed star of \a S minus the stars of all faces of \a
    * S.
    *
    * @tparam CellAssociativeContainer any simple associative container with key Cell.
    * @param S any collection of cells of this complex.
    * @param hintClosed when 'true', this hint tells that the complex
    * is (locally around \a S) closed, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    * @param hintOpen when 'true', this hint tells that the complex
    * is (locally around \a S) open, so this speeds up this method, otherwise, the
    * complex may be arbitrary.
    * @return the link of \a S within this complex.
    */
    template <typename CellAssociativeContainer>
    CellAssociativeContainer link( const CellAssociativeContainer& S, bool hintClosed = false, bool hintOpen = false ) const;

    // ----------------------- global operations on complexes --------------------------------------
  public:

    /**
    * Close the whole complex.
    */
    void close();

    /**
    * Close all cells of dimension less or equal to \a k.
    * @param k any strictly positive integer.
    */
    void close( Dimension k );

    /**
    * Collapse a user-specified part of this complex, collapsing cells
    * following priority [priority], in a decreasing sequence until no
    * more collapse is feasible. The range [\a S_itb,\a S_itE)
    * provides the starting cells, generally (but not compulsory)
    * maximal cells. The resulting complex is guaranteed to keep the
    * same homotopy type (a kind of topology equivalence).
    *
    * @note Cells whose data has been marked as FIXED are not removed.
    *
    * @note Only cells that are in the closure of [\a S_itb,\a S_itE)
    * may be removed, and only if they are not marked as FIXED.
    *
    * @advanced If you use a DefaultCellMapIteratorPriority object as
    * \a priority, then the VALUE part of each cell data defines the
    * priority (the highest value the soonest are these cells
    * collapsed). You may thus fill these cell values before calling
    * this method.
    *
    * @tparam CellConstIterator any forward const iterator on Cell.
    *
    * @tparam CellMapIteratorPriority any type defining a method 'bool
    * operator()( const Cell&, const Cell&) const'. Defines the order
    * in which cells are collapsed. @see DefaultCellMapIteratorPriority
    *
    * @param S_itB the start of a range of cells which is included in [K].
    * @param S_itE the end of a range of cells which is included in [K].
    * @param priority the object that assign a priority to each cell.
    * @param hintIsSclosed indicates if [\a S_itb,\a S_ite) is a closed set (faster in this case).
    * @param hintIsKclosed indicates that this complex is closed.
    * @param verbose outputs some information during processing when 'true'.
    */
    template <typename CellConstIterator, typename CellMapIteratorPriority>
    void collapse( CellConstIterator S_itB, CellConstIterator S_itE, 
                   const CellMapIteratorPriority& priority, 
                   bool hintIsSClosed = false, bool hintIsKClosed = false,
                   bool verbose = false );


    // ----------------------- Interface --------------------------------------
  public:

    /**
    * Writes/Displays the object on an output stream.
    * @param out the output stream where the object is written.
    */
    void selfDisplay ( std::ostream & out ) const;

    /**
    * Checks the validity/consistency of the object.
    * @return 'true' if the object is valid, 'false' otherwise.
    */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  protected:

    /// The Khalimsky space in which lives the cubical complex.
    const KSpace* myKSpace;

    // ------------------------- Private Datas --------------------------------
  private:

    /// An array of map Cell -> Data that stores cells dimension per
    /// dimension (i.e. cells of dimension 0 are stored in myCells[0],
    /// cells of dimension 1 in myCells[1] and so on).
    std::vector<CellMap> myCells;
    

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
    * Constructor.
    * Forbidden by default (protected to avoid g++ warnings).
    */
    CubicalComplex();

  private:


    // ------------------------- Internals ------------------------------------
  private:

    /**
    * Given a cell [c], tells if it is a maximal cell in the complex
    * (return 0), or if it is a free face of the cell pointed by
    * [it_cell_up] (return 1) or if it is not a free face.
    *
    * The complex must be closed. In computing the 1-up-incident
    * cells, this method ignores cell marked as REMOVED. Furthermore,
    * if one 1-up-incident cell is not marked as COLLAPSIBLE, the
    * method returns 2.
    *
    * @param c a cubical cell (belonging to 'this')
    *
    * @param it_cell_up (returns) a pointer on a cell d if c is a
    * free face of d.
    *
    * @return 0 if the cell is maximal, 1 if the cell is a free face,
    * 2 otherwise.
    */
    uint32_t computeCellType( Dimension n, const Cell& c, CellMapIterator& it_cell_up );



  }; // end of class CubicalComplex

  namespace detail 
  {
    template <typename AssociativeContainer, bool ordered>
    struct SetOperation {
      /** 
      * Updates the set S1 as S1 - S2. This version does not use the
      * fact that the container is ordered.
      * @param[in,out] S1 an input set, \a S1 - \a S2 as output.
      * @param[in] S2 another input set.
      */
      static void difference( AssociativeContainer& S1, const AssociativeContainer& S2 )
      {
        for ( typename AssociativeContainer::const_iterator it = S2.begin(), 
          itE = S2.end(); it != itE; ++it )
          S1.erase( *it );
      }
    };

    template <typename AssociativeContainer>
    struct SetOperation< AssociativeContainer, true > {
      /** 
      * Updates the set S1 as S1 - S2. This version uses the fact that
      * the container is ordered.
      * @param[in,out] S1 an input set, \a S1 - \a S2 as output.
      * @param[in] S2 another input set.
      */
      static void difference( AssociativeContainer& S1, const AssociativeContainer& S2 )
      {
        AssociativeContainer S;
        std::swap( S,S1 );
        std::set_difference( S.begin(), S.end(), S2.begin(), S2.end(), std::inserter( S1, S1.end() ) );
      }
    };

  }

  /**
  * Overloads 'operator<<' for displaying objects of class 'CubicalComplex'.
  * @param out the output stream where the object is written.
  * @param object the object of class 'CubicalComplex' to write.
  * @return the output stream after the writing.
  */
  template <typename TKSpace, typename TCellContainer>
  std::ostream&
  operator<< ( std::ostream & out, 
               const CubicalComplex<TKSpace, TCellContainer> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/topology/CubicalComplex.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CubicalComplex_h

#undef CubicalComplex_RECURSES
#endif // else defined(CubicalComplex_RECURSES)
