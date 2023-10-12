// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#pragma once

#include <array>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include "fast_exp.h"
#include "fwd_decl.h"
#include "random.h"

namespace trento3d {

// TODO explain this class structure

/// \rst
/// Represents a single nucleon.  Stores its transverse position and whether or
/// not it's a participant.  These properties are globally readable, but can
/// only be set through ``Nucleus`` and ``NucleonProfile``.
/// \endrst
class NucleonData {
 public:
  /// Only a default constructor is necessary\---the class is designed to be
  /// constructed once and repeatedly updated.
  NucleonData() = default;

  /// Whether or not this nucleon is a participant.
  bool is_participant() const;

  /// Whether or not its constituents have been sampled.
  bool constituents_exist() const;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// The longitudinal \em z position.
  double z() const;

 private:
  ///
  friend class NucleonCommon;

  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;

  /// Set the transverse position and reset participant status to false.
  void set_position(double x, double y, double z);

  /// Internal storage of the transverse position.
  double x_, y_, z_;

  /// Whether or not nucleon is a participant.
  bool is_participant_ = false;

  /// Whether or not nucleon's constituents are sampled.
  bool constituents_exist_ = false;

  /// Constituent transverse position and fluctuation prefactor.
  struct Constituent {
    double x, y;
    double frac_mid;
    double frac_forward;
  };

  /// Vector of constituent positions and fluctuation prefactors.
  std::vector<Constituent> constituents_;
};

/// \rst
/// Encapsulates properties shared by all nucleons: transverse thickness
/// profile, cross section, fluctuations.  Responsible for sampling
/// nucleon-nucleon participation with given `\sigma_{NN}`.
/// \endrst
class NucleonCommon {
 public:
  /// Instantiate from the configuration.
  explicit NucleonCommon(const VarMap& var_map);

  /// The maximum impact parameter for participation.
  double max_impact() const;

  /// Corners of the tile enclosing the nucleon thickness.
  std::array<double, 4> boundary(const NucleonData& nucleon) const;

  /// Nucleon thickness as a function of transverse position, at mid-rapidity
  double thickness(const NucleonData& nucleon, double x, double y) const;

  /// WK: same as above, but without the Gamma fluctuation,
  /// used in the calculation of binary collision density
  double deterministic_thickness(const NucleonData& nucleon, double x, double y) const;

  /// WK: return Tpp given bpp^2
  double norm_Tpp(double bpp_sqr) const;

  /// Nucleon thickness as a function of transverse position, at forward rapidity
  double fragmentation(const NucleonData& nucleon, double x, double y) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(NucleonData& A, NucleonData& B) const;
  double avg_xloss() const {return avg_xloss_;};

 private:
  double avg_xloss_;

  /// Sample constituent positions inside the nucleon.
  void sample_constituent_positions(NucleonData& nucleon) const;

  /// Flag the nucleon as a participant.
  void set_participant(NucleonData& nucleon) const;

  /// Fast exponential for calculating the thickness profile.
  const FastExp<double> fast_exp_;

  /// Inelastic form factor width:
  const double form_width_;

  /// Gaussian nucleon width.
  const double nucleon_width_;

  /// Gaussian constituent width.
  const double constituent_width_;

  /// Possibly non-integer number of constituents.
  const double constituent_number_raw_;

  /// Minimum number of constituents inside the nucleon.
  const std::size_t constituent_number_;

  /// Probability of having an additional constituent.
  const double extra_constituent_probability_;

  /// Gaussian width of constituent position sampling distribution.
  const double sampling_width_;

  /// Maximum impact parameter for nucleon-nucleon participation.
  const double max_impact_sq_;

  /// Constituent width squared.
  const double constituent_width_sq_;

  /// Calculate thickness out to this distance from constituent center.
  const double constituent_radius_sq_;

  /// Width of Gaussian thickness function.
  const double width_sqr_;

  /// WK 1/4pi
  const double one_div_four_pi_;

  /// Nuclear opacity parameter
  double sigma_partonic_;
  //double sigma_partonic_lonc_, sigma_partonic_hinc_;

  /// Tracks binary collisions if true
  const bool calc_ncoll_;

  /// Gaussian distribution for sampling constituent positions.
  mutable std::normal_distribution<double> constituent_position_dist_;

  mutable boost::random::beta_distribution<> constituent_xloss_dist_;

};

// These functions are short, called very often, and account for a large
// fraction of the total computation time, so request inlining.

// Trivial helper function.
template <typename T>
inline constexpr T sqr(const T& value) {
  return value * value;
}

// NucleonData inline member functions

inline double NucleonData::x() const {
  return x_;
}

inline double NucleonData::y() const {
  return y_;
}

inline double NucleonData::z() const {
  return z_;
}

inline void NucleonData::set_position(double x, double y, double z) {
  x_ = x;
  y_ = y;
  z_ = z;
  is_participant_ = false;
  constituents_exist_ = false;
}

inline bool NucleonData::is_participant() const {
  return is_participant_;
}

inline bool NucleonData::constituents_exist() const {
  return constituents_exist_;
}

// NucleonCommon inline member functions

inline double NucleonCommon::max_impact() const {
  return std::sqrt(max_impact_sq_);
}

inline std::array<double, 4>
NucleonCommon::boundary(const NucleonData& nucleon) const {
  auto constituent = nucleon.constituents_.begin();

  // Initialize the boundary with the position of the first constituent.
  auto xmin = constituent->x, xmax = constituent->x;
  auto ymin = constituent->y, ymax = constituent->y;

  // Check the remaining constituents and update the boundary accordingly.
  // Using this instead of something from std::algorithm because it finds all
  // four quantities {xmin, xmax, ymin, ymax} in a single pass over the constituents.
  for (std::advance(constituent, 1); constituent != nucleon.constituents_.end(); ++constituent) {
    auto x = constituent->x;
    if (x < xmin)
      xmin = x;
    else if (x > xmax)
      xmax = x;

    auto y = constituent->y;
    if (y < ymin)
      ymin = y;
    else if (y > ymax)
      ymax = y;
  }

  // Remember to add and subtract the constituent radius.
  auto r = std::sqrt(constituent_radius_sq_);

  return {xmin - r, xmax + r, ymin - r, ymax + r};
}

inline double NucleonCommon::thickness(
    const NucleonData& nucleon, double x, double y) const {
  auto t = 0.;

  for (const auto& constituent : nucleon.constituents_) {
    auto frac = constituent.frac_mid;
    auto distance_sq = sqr(x - constituent.x) + sqr(y - constituent.y);
    if (distance_sq < constituent_radius_sq_)
      t += frac * fast_exp_(-.5*distance_sq/constituent_width_sq_);
  }

  return math::double_constants::one_div_two_pi
        /constituent_width_sq_/constituent_number_raw_
        * t;
}

inline double NucleonCommon::deterministic_thickness(
    const NucleonData& nucleon, double x, double y) const {
  auto t = 0.;

  for (const auto& constituent : nucleon.constituents_) {
    auto distance_sq = sqr(x - constituent.x) + sqr(y - constituent.y);
    if (distance_sq < constituent_radius_sq_)
      t += fast_exp_(-.5*distance_sq/constituent_width_sq_);
  }

  return math::double_constants::one_div_two_pi
        /constituent_width_sq_/constituent_number_raw_
        * t;
}

inline double NucleonCommon::norm_Tpp(double bpp_sqr) const  {
  return one_div_four_pi_
         / constituent_width_sq_
         * std::exp(-.25 * bpp_sqr / constituent_width_sq_);  // TODO: verify this, and use fast_exp_ (will need to adjust xmin in constructor)
}

inline double NucleonCommon::fragmentation(
    const NucleonData& nucleon, double x, double y) const {
  auto t = 0.;

  for (const auto& constituent : nucleon.constituents_) {
    auto frac = constituent.frac_forward;
    auto distance_sq = sqr(x - constituent.x) + sqr(y - constituent.y);
    if (distance_sq < constituent_radius_sq_)
      t += frac * fast_exp_(-.5*distance_sq/constituent_width_sq_);
  }

  return math::double_constants::one_div_two_pi
        /constituent_width_sq_/constituent_number_raw_
        * t;
}

inline bool NucleonCommon::participate(NucleonData& A, NucleonData& B) const {
  // If both nucleons are already participants, there's nothing to do.
  if (A.is_participant() && B.is_participant() && !calc_ncoll_)
    return true;

  auto distance_sq = sqr(A.x() - B.x()) + sqr(A.y() - B.y());

  // Check if nucleons are out of range.
  if (distance_sq > max_impact_sq_)
    return false;

  sample_constituent_positions(A);
  sample_constituent_positions(B);

  double overlap = std::exp(-distance_sq/4./form_width_/form_width_);
  auto one_minus_prob  = std::exp(-sigma_partonic_*overlap*one_div_four_pi_/form_width_/form_width_);

  // Sample one random number and decide if this pair participates.
  if (one_minus_prob < random::canonical<double>()) {
    set_participant(A);
    set_participant(B);
    return true;
  }

  return false;
}

inline void NucleonCommon::sample_constituent_positions(NucleonData& nucleon) const {
  if (nucleon.constituents_exist())
    return;

  auto nc = constituent_number_ + (((extra_constituent_probability_> 0.) && (random::canonical<double>() < extra_constituent_probability_)) ? 1 : 0);  // for compatibility, don't draw from PRNG if n_c is a whole number
  nucleon.constituents_.resize(nc);

  double xcom = 0.0;
  double ycom = 0.0;

  // Sample nucleon constituent positions
  double overall_fluct = constituent_xloss_dist_(random::engine);
  for (auto&& constituent : nucleon.constituents_) {
    auto xloc = constituent_position_dist_(random::engine);
    auto yloc = constituent_position_dist_(random::engine);

    constituent.x = xloc;
    constituent.y = yloc;
    constituent.frac_mid = overall_fluct;
    constituent.frac_forward = 1. - constituent.frac_mid;

    xcom += xloc;
    ycom += yloc;
  }

  xcom /= nc;
  ycom /= nc;

  // Place nucleon at the desired position
  for (auto&& constituent : nucleon.constituents_) {
    constituent.x += (nucleon.x() - xcom);
    constituent.y += (nucleon.y() - ycom);
  }

  nucleon.constituents_exist_ = true;
}

inline void NucleonCommon::set_participant(NucleonData& nucleon) const {
  if (nucleon.is_participant())
    return;

  nucleon.is_participant_ = true;
}

}  // namespace trento3d
