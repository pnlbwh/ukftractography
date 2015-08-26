

/**
 * \file filter_model.cc
 * \brief Implements functions defined in filter_model.h
*/

#include "filter_model.h"
#include <iostream>
#include "utilities.h"

ukfPrecisionType FilterModel::CheckZero(const ukfPrecisionType & local_d) const
{
  if( local_d < 0 )
    {
    if( local_d >= -1.0e-4 ) // for small errors just round it to 0
      {
      return ukfZero;
      }
    else   // for errors too big exit with exception
      {
      std::cout << "Error, a variable became negative. Most likely something went wrong in the QP\n";
      exit(1);
      }
    }
  return local_d;
}

void createProtocol(const ukfVectorType& _b_values,
                    ukfVectorType& _gradientStrength, ukfVectorType& _pulseSeparation)
{
  static int count=0;
  static ukfVectorType gradientStrength, pulseSeparation;
  std::vector<double> Bunique, tmpG;
  ukfPrecisionType Bmax = 0;
  ukfPrecisionType tmp, Gmax, GAMMA;

  _gradientStrength.resize(_b_values.size());
  _pulseSeparation.resize(_b_values.size());

  // set maximum G = 40 mT/m
  Gmax = 0.04;
  GAMMA = 267598700;
  if(count ==1)
  {
     // gradient strength and pulseSeparation are once
     // same values are returned
    _gradientStrength = gradientStrength;
    _pulseSeparation = pulseSeparation;
  }
  else{
    for(int i = 0; i < _b_values.size(); ++i )
    {
      int unique = 1;
      for(size_t j = 0; j < Bunique.size(); ++j )
      {
        if (_b_values[i] == Bunique[j])
        {
          unique = 0;
          break;
        }
      }
      if(unique == 1)
      {
        Bunique.push_back(_b_values[i]);
      }
      if(Bmax < _b_values[i])
      {
        Bmax = _b_values[i];
      }
    }

    tmp = cbrt(3*Bmax*1000000/(2*GAMMA*GAMMA*Gmax*Gmax));

    for(int i = 0; i < _b_values.size(); ++i )
    {
      _pulseSeparation[i] = tmp;
    }

    for(size_t i = 0; i < Bunique.size(); ++i )
    {
      tmpG.push_back(std::sqrt(Bunique[i]/Bmax) * Gmax);
    }

    for(size_t i = 0; i < Bunique.size(); ++i )
    {
      for(int j=0; j < _b_values.size(); j++)
      {
        if(_b_values[j] == Bunique[i])
        {
          _gradientStrength[j] = tmpG[i];
        }
      }
    }
    gradientStrength = _gradientStrength;
    pulseSeparation = _pulseSeparation;
    count = 1;
  }
}

