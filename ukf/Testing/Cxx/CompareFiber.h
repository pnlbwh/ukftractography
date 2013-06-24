#ifndef COMPAREFIBER_H_
#define COMPAREFIBER_H_

#include <vector>
#include <cassert>
#include "linalg.h"

struct Fiber
  {

public:

  typedef std::map<std::string, std::vector<float> > FieldMapType;
  std::vector<vec_t> Points;
  FieldMapType Fields;

  };

#endif
