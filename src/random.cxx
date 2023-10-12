// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#include "random.h"


namespace trento3d { namespace random {

// Seed random number generator from hardware device.
Engine engine{std::random_device{}()};

}}  // namespace trento3d::random
