/**
 * \file FiberLogic.h
 * \brief Contains class FiberLogic which defines logical operators on VTK fibers.
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef FIBERLOGIC_H_
#define FIBERLOGIC_H_

#include <vector>

/**
 * \class FiberLogic
 * \brief Contains logical operators for fibers
 */
class FiberLogic
{
public:

  /** Constructor */
  FiberLogic()
  {
  }

  /** Destructor */
  ~FiberLogic();

  /**
   * Checks wheather a fiber is contained in a vector of fibers
   * \param fib a single fiber
   * \param fibers a vector of fibers
   * \returns true if fib is contained in fibers, false otherwise.
  */
  inline bool Contains(std::vector<Fiber> & fibers, Fiber & fib)
  {
    for( unsigned int i = 0; i < fibers.size(); i++ )
      {
      if( Equals(fibers[i], fib) )
        {
        return true;
        }
      }
    return false;
  }

  /** Forms the logical OR of two fiber vectors, the result is A plus all elements in B that aren't already in A */
  inline std::vector<Fiber> doOr(std::vector<Fiber> & A, std::vector<Fiber> & B)
  {

    std::vector<Fiber> out;

    out = A;
    for( unsigned int i = 0; i < B.size(); i++ )
      {
      if( !Contains(A, B[i]) ) // add just the ones that aren't already in A
        {
        out.push_back(B[i]);
        }
      }

    return out;
  }

  /** Forms the logical AND of two fibers, returns all fibers that are in the fiber vector A and the fiber vector B */
  inline std::vector<Fiber> doAnd(std::vector<Fiber> & A, std::vector<Fiber> & B)
  {
    std::vector<Fiber> out;

    for( unsigned int i = 0; i < B.size(); i++ )
      {
      if( Contains(A, B[i]) )
        {
        out.push_back(B[i]);
        }
      }
    return out;
  }

  /** The equal operator for two fibers, two fibers are equal if all points are equal */
  inline bool Equals(Fiber & fib1, Fiber & fib2)
  {
    if( fib1.Points.size() == fib2.Points.size() )
      {
      for( unsigned int i = 0; i < fib1.Points.size(); ++i )
        {
        if( (fib1.Points[i][0] != fib2.Points[i][0]) &&
            (fib1.Points[i][1] != fib2.Points[i][1]) &&
            (fib1.Points[i][2] != fib2.Points[i][2]) )
          {
          return false;
          }
        }
      return true;
      }
    else
      {
      return false;
      }
  }

};

#endif  // FIBERCOLLECTION_H_
