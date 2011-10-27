// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_HDNUM_HH
#define HDNUM_HDNUM_HH

#define HDNUM_DEVEL_MODE 1

#define HDNUM_HAS_GMP 1

#include <math.h>
#include <complex>
#if HDNUM_HAS_GMP
#include <gmpxx.h>
#endif

#include "src/exceptions.hh"
#include "src/vector.hh"
#include "src/densematrix.hh"
#include "src/timer.hh"
#include "src/precision.hh"
#include "src/lr.hh"
#include "src/ode.hh"
#include "src/pde.hh"
#include "src/sgrid.hh"
#include "src/newton.hh"


#endif
