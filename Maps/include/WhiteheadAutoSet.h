
#ifndef _WHITEHEAD_AUTO_SET_H_
#define _WHITEHEAD_AUTO_SET_H_


#include "Map.h"
#include "Word.h"

#include <ext/hash_set>
//#include <ext/stl_hash_fun.h>


// --------------------- structures to create a set of maps2
//! Implements a comparison operator of two maps.
struct compMaps
{
  bool operator()(  Map m1,  Map m2 ) const
  {
    return (m1 == m2);
  }
};

namespace __gnu_cxx {
  struct map_hash
  {
    size_t operator()(Map m) const
    {
      int h = 0;
      const vector<Word>& is =  m.generatingImages( );
      for (int i=0;i<is.size();i++)
	h+=is[i].length();
      
      return h;
    }
  };
}

//! Implements a set fo Map objects
typedef __gnu_cxx::hash_set<Map, __gnu_cxx::map_hash, compMaps> SetOfMaps;


//! Abstract interface for a set of Maps (Automorphisms)
class AutoSet
{
 public:
  //! Returns the corresponding set of maps
  virtual const SetOfMaps& getSet() const = 0;
};


// -------------------------- NielsenAutoSet ---------------------------- //

//! Implements a set of Nielsen automorphisms of a free group.
class NielsenAutoSet : public AutoSet
{

public:

  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Constructors:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  //! Constructor. Generates a set of Nielsen automorphisms
  /*!
    The set of Nielsen automorphisms for a free group of rank \c n  is generated. <br>
    Note, the whole set is generated by enumeration. \f$\O(n^2)\f$ automorphisms.
    \param n - rank of a free group.
   */
  NielsenAutoSet( int  n );

  ~NielsenAutoSet( ) { }

  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Accessors:                                                          //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  //! Returns a random Nielsen automorphism.
  /*!
    \return Nielsen automorphisms selected uniformly randomly form the set.
  */
  const Map& getRandomAuto()const;

  //! Returns the set of Nielsen automorphisms.
  /*!
    \return the set of Nielsen automorphisms.
   */
  const SetOfMaps& getSet() const { return theSet; }
 private:
  SetOfMaps theSet;
  int nGens;
};

// -------------------------- RestrictedWhiteheadAutoSet ---------------------------- //


//! Implements a so-called restricted set of Whitehead automorphisms of a free group.
/*!
  In general the Restricted set of Whitehead automorphisms contains automorphisms
  of the type: <br>
  \f[ x_i \rightarrow x_i x_j \f]
  \f[ x_i \rightarrow x^{-1}_j x_i \f]
  \f[ x_i \rightarrow x^{-1}_j x_i x_j, \f]
where \f$ x_j, x_j \in X\f$ - generating set of the free group.
 */
class RestrictedWhiteheadAutoSet : public AutoSet
{

public:

  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Constructors:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  
  //! Constructor. Generates a set of restricted Whitehead  automorphisms
  /*!
    The set of restricted Whitehead automorphisms for a free group of rank \c n  is generated. <br>
    Note, the whole set is generated by enumeration. 
    \param n - rank of a free group.
    \param use_conj - if \c false the conjugation is excluded from the set. \c true is the default value.
   */
  RestrictedWhiteheadAutoSet( int n, bool use_conj = true );

  ~RestrictedWhiteheadAutoSet( ) { }

  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Accessors:                                                          //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////


  //! Returns a random restricted Whitehead automorphism.
  /*!
    \return restricted Whitehead automorphisms selected uniformly randomly form the set.
  */ 
  const Map& getRandomAuto()const;
  
  //! Returns the set of restricted Whitehead automorphisms.
  /*!
    \return the set of restricted Whitehead automorphisms.
  */
  const SetOfMaps& getSet() const { return theSet; }
 private:
  SetOfMaps theSet;
  int nGens;
};

//
// -------------------------- WhiteheadAutoSetType2 ---------------------------- //
//

//! Implements the  set of Whitehead automorphisms of type II.
/*!
  Whitehead automorphisms of the type II are the Whitehead automorphisms which may 
  alter the length of an input word. Basically it is all the Whitehead automorphisms 
  excluding permutations of the generators.
*/
class WhiteheadAutoSetType2 : public AutoSet
{
  
 public:
  
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Constructors:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  //! Constructor. Generates the set of  Whitehead  automorphisms of type II.
  /*!
    The set of Whitehead automorphisms for a free group of rank \c n  is generated. <br>
    Note, the whole set is generated by enumeration. \f$ \O(2^n)\f$ automorphisms. This is
    a brute force approach and inapplicable for groups with large ranks.
    \param n - rank of a free group.
   */
  WhiteheadAutoSetType2(  int n );

  ~WhiteheadAutoSetType2( );

  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Accessors:                                                          //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////


 //! Returns a random Whitehead automorphism.
  /*!
    \return  Whitehead automorphisms selected uniformly randomly form the set.
  */ 
  const SetOfMaps& getSet() const { return theSet; }
  //! Returns the set of Whitehead automorphisms.
  /*!
    \return the set of Whitehead automorphisms.
  */
  const Map& getRandomAuto() const;
  
private:

  Map getMap(   int n, const vector<int>& tCounts, Word a );
  void computeSet( int n );

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Data Members:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  const static int nElemAutos = 4;

  SetOfMaps theSet;
  int nGens;

};



//! Implements a greedy  procedure of reducing a word to its minimal length.
/*!
  Greedy procedure applies automorphisms from a  Whitehead set of type II  to a word until 
  its length is reduced. If a shorter image is found, procedure is applied to the
  shorten word. If none of the automorphisms can reduce the length, procedure stops and
  the shortest word is returned.
 */
class WhiteheadMinimization
{
 public:
  //! Constructor. 
  /*! A set of Whitehead automorphisms is created in the constructor. 
    May not be applicable when the rank is large.
    \param n - the rank  of a free group.
  */
  WhiteheadMinimization( int n ): wSet( n ) { } 

  //! Test if the word is minimal
  /*!
    Applies automorphisms from the Whitehead set to test if a given word has
    minimal length.
   \param w - word to be tested.
   \return \c true if \c w has the minimal length, \c false otherwise.
   */
  bool isMinimal( const Word& w) const;
  
  
  //! Find a word of the minimal length 
  /*!
    Applies automorphisms from the Whitehead set to find a word of the minimal length
    in the automorphic orbit of a given word.
    \param w - initial word.
    \return a word of the minimal length in the orbit of \c w.
   */
  Word findMinimal( const Word& w, ostream* out = NULL )const;

  //! Get the Whitehead set of type II.
  /*!
    \return the Whitehead set of type II.
   */
  const WhiteheadAutoSetType2& getSet() const { return wSet; }

 private:
  WhiteheadAutoSetType2 wSet;
};


// ---------------------------------------------------------------------------------------- //
// Check if a word is reducible by a given set
// ---------------------------------------------------------------------------------------- //

namespace WhiteheadAutoSet
{
  //! Reduce the length of a given word by automorphisms form a set.
  /*!
    Tries to reduce the length of a given word by applying automorphisms
    from a set. Will return if a shorter word is found <br>
    Note, it does not return the minimal length word.
    \param w - the initial word to be reduced.
    \param theSet - set of automorphisms to reduce by.
    \return a shorter word if found, original word otherwise.
   */ 
  Word reduceBy( const Word& w, const SetOfMaps& theSet );
};

#endif

