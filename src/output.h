// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "fwd_decl.h"

namespace trento3d {

/// Simple interface for outputting event data.  Determines which output formats
/// to create based on the configuration and writes those formats when called.
class Output {
 public:
  /// Instantiate from the configuration.
  Output(const VarMap& var_map);

  /// \rst
  /// Call the functor to output event data.  Arguments are perfect-forwarded to
  /// each output function.  The required arguments are
  ///
  /// - ``int`` event number
  /// - ``double`` impact parameter
  /// - ``int`` binary collisions
  /// - ``const Event&`` Event object
  ///
  /// \endrst
  template <typename... Args>
  void operator()(Args&&... args) const;

 private:
  /// Internal storage of output functions.
  std::vector<std::function<void(int, double, int, const Event&)>> writers_;

  /// Buffer for composing binary-format output.
  std::vector<float> buf_;
};

template <typename... Args>
void Output::operator()(Args&&... args) const {
  for (const auto& write : writers_)
    write(std::forward<Args>(args)...);
}

}  // namespace trento3d
