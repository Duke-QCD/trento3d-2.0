// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#include "nucleon.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include "boost/random.hpp"
using boost::math::quadrature::trapezoidal;

#include "fwd_decl.h"

namespace trento3d {

namespace {

// These constants define distances in terms of the width of the nucleon profile
// Gaussian thickness function.

// Truncation radius of the thickness function.
constexpr double max_radius_widths = 5.;

// Maximum impact parameter for participation.
constexpr double max_impact_widths = 6.;
// Create ctor parameters for unit mean std::gamma_distribution.
//   mean = alpha*beta == 1  ->  beta = 1/alpha
// Used below in NucleonProfile ctor initializer list.

template <typename RealType> using param_type =
  typename std::gamma_distribution<RealType>::param_type;

template <typename RealType>
param_type<RealType> gamma_param_unit_mean(RealType alpha = 1.) {
  return param_type<RealType>{alpha, 1./alpha};
}

// fitted sigma_nn_inel to PDG data within 5 GeV < sqrt{s} < 10^5 GeV
// units: [fm^2]
double SigmaNN(double sqrts){
  double logsqrts = std::log(sqrts);
  return (0.073491826*logsqrts-.19313457)*logsqrts+3.123737545;
}

// Test approximate cross section parameter equality.
bool almost_equal(double a, double b) {
  return fabs(a - b) < 1e-6;
}

// Calculate constituent position sampling width from CL arguments
double calc_sampling_width(const VarMap& var_map) {
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto constituent_width = var_map["constit-width"].as<double>();
  auto constituent_number = var_map["constit-number"].as<double>();

  auto one_constituent = (constituent_number <= 1.);
  auto same_size = (nucleon_width - constituent_width < 1e-6);

  if (one_constituent || same_size)
    return 0.0;
  else {
    auto num = sqr(nucleon_width) - sqr(constituent_width);
    auto denom = 1. - 1./constituent_number;
    return std::sqrt(num / denom);
  }
}

// Determine cross section parameter for sampling nucleon participants.
// This semi-analytic method is only valid for one constituent.
double analytic_partonic_cross_section(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto sigma_nn = SigmaNN(var_map["sqrts"].as<double>());
  //auto width = var_map["nucleon-width"].as<double>();
  auto width = var_map["form-width"].as<double>();

  // Initialize arguments for boost root finding function.

  // Bracket min and max.
  auto a = -10.;
  auto b = 20.;

  // Tolerance function.
  // Require 3/4 of double precision.
  math::tools::eps_tolerance<double> tol{
    (std::numeric_limits<double>::digits * 3) / 4};

  // Maximum iterations.
  // This is overkill -- in testing only 10-20 iterations were required
  // (but no harm in overestimating).
  boost::uintmax_t max_iter = 1000;

  // The right-hand side of the equation.
  auto rhs = sigma_nn / (4 * math::double_constants::pi * sqr(width));

  // This quantity appears a couple times in the equation.
  auto c = sqr(max_impact_widths) / 4;

  try {
    auto result = math::tools::toms748_solve(
      [&rhs, &c](double x) {
        using std::exp;
        using math::expint;
        return c - expint(-exp(x)) + expint(-exp(x-c)) - rhs;
      },
      a, b, tol, max_iter);

    auto cross_section_param = .5*(result.first + result.second);
    return std::exp(cross_section_param) * 4 * math::double_constants::pi * sqr(width);
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon width too small?"};
  }
}

}  // unnamed namespace

NucleonCommon::NucleonCommon(const VarMap& var_map)
    : fast_exp_(-.5*sqr(max_radius_widths), 0., 1000),
      form_width_(var_map["form-width"].as<double>()),
      nucleon_width_(var_map["nucleon-width"].as<double>()),
      constituent_width_(var_map["constit-width"].as<double>()),
      constituent_number_raw_(var_map["constit-number"].as<double>()),
      constituent_number_(std::size_t(std::floor(constituent_number_raw_))),
      extra_constituent_probability_(constituent_number_raw_ - constituent_number_),
      sampling_width_(calc_sampling_width(var_map)),
      max_impact_sq_(sqr(max_impact_widths*nucleon_width_)),
      constituent_width_sq_(sqr(constituent_width_)),
      constituent_radius_sq_(sqr(max_radius_widths*constituent_width_)),
      width_sqr_(sqr(nucleon_width_)),
      one_div_four_pi_(0.5*math::double_constants::one_div_two_pi),
      sigma_partonic_(analytic_partonic_cross_section(var_map)),
      calc_ncoll_(var_map["ncoll"].as<bool>()),
      constituent_position_dist_(0, sampling_width_){
    double Mproton = 0.938;
    double sqrts = var_map["sqrts"].as<double>();
    double eta_max = std::log(sqrts/var_map["kT-min"].as<double>());
    double mid_power = var_map["mid-power"].as<double>();
    double mid_norm = var_map["mid-norm"].as<double>();
    double var = var_map["fluctuation"].as<double>();
    double flatness_ = var_map["flatness"].as<double>();
    double norm_trento = mid_norm * Mproton * std::pow(sqrts/Mproton, mid_power);
    auto f1 = [eta_max,flatness_](double eta){
        //return std::exp(-x*x/2./L_)*std::cosh(x)*std::pow(1.-std::pow(x/eta_max, 4), 2);
        if (std::abs(eta)>eta_max) return 0.;
        double u = eta*eta/2./(eta_max-3.0);
        return std::exp(-std::pow(u, flatness_)) * std::cosh(eta)
             * std::pow(1.-std::pow(eta/eta_max, 4), 4);
    };
    double F1 = norm_trento * trapezoidal(f1, -eta_max, eta_max);
    // Average energy fraction  needed to be depositied from the fragmentation region
    avg_xloss_ = F1/sqrts;
    if (avg_xloss_>1) {
      std::cout << "central fireball too large!" << std::endl;
      exit(-1);
    }
    //std::cout <<  avg_xloss_ << std::endl;
    double sum = 1./var - 1.;
    double a = avg_xloss_ * sum;
    double b = (1. - avg_xloss_) * sum;
    if (a<0 ) {
      std::cout << "Fluctuation too large" << std::endl;
      exit(-1);
    }
    constituent_xloss_dist_ = boost::random::beta_distribution<>(a, b);
}

}  // namespace trento3d
