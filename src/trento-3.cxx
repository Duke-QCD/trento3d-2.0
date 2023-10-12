// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "collider.h"
#include "fwd_decl.h"

// CMake sets this definition.
// Fall back to a sane default.
#ifndef TRENTO3D_VERSION_STRING
#define TRENTO3D_VERSION_STRING "dev"
#endif

namespace trento3d {

namespace {

void print_version() {
  std::cout << "trento-3 " << TRENTO3D_VERSION_STRING << '\n';
}

void print_bibtex() {
  std::cout <<
    "@article{Moreland:2014oya,\n"
    "    author         = \"Moreland, J. Scott and Bernhard, Jonah E. and Bass,\n"
    "                      Steffen A.\",\n"
    "    title          = \"{Alternative ansatz to wounded nucleon and binary\n"
    "                      collision scaling in high-energy nuclear collisions}\",\n"
    "    journal        = \"Phys.Rev.\",\n"
    "    number         = \"1\",\n"
    "    volume         = \"C92\",\n"
    "    pages          = \"011901\",\n"
    "    doi            = \"10.1103/PhysRevC.92.011901\",\n"
    "    year           = \"2015\",\n"
    "    eprint         = \"1412.4708\",\n"
    "    archivePrefix  = \"arXiv\",\n"
    "    primaryClass   = \"nucl-th\",\n"
    "    SLACcitation   = \"%%CITATION = ARXIV:1412.4708;%%\",\n"
    "}\n"
    "\n"
    "@article{Moreland:2018gsh,\n"
    "    author         = \"Moreland, J. Scott and Bernhard, Jonah E. and Bass,\n"
    "                      Steffen A.\",\n"
    "    title          = \"{Bayesian calibration of a hybrid nuclear collision model using p-Pb and Pb-Pb data at energies available at the CERN Large Hadron Collider}\",\n"
    "    eprint         = \"1808.02106\",\n"
    "    archivePrefix  = \"arXiv\",\n"
    "    primaryClass   = \"nucl-th\",\n"
    "    doi            = \"10.1103/PhysRevC.101.024911\",\n"
    "    journal        = \"Phys. Rev. C\",\n"
    "    volume         = \"101\",\n"
    "    number         = \"2\",\n"
    "    pages          = \"024911\",\n"
    "    year           = \"2020\"\n"
    "}\n"
    "\n"
    "@misc{Soeder:2023vdn,\n"
    "    author         = \"Soeder, Derek and Ke, Weiyao and Paquet, J.-F. and Bass,\n"
    "                      Steffen A.\",\n"
    "    title          = \"{Bayesian parameter estimation with a new three-dimensional initial-conditions model for ultrarelativistic heavy-ion collisions}\",\n"
    "    eprint         = \"2306.08665\",\n"
    "    archivePrefix  = \"arXiv\",\n"
    "    primaryClass   = \"nucl-th\",\n"
    "    month          = \"6\",\n"
    "    year           = \"2023\"\n"
    "}\n";
}

// TODO
// void print_default_config() {
//   std::cout << "to do\n";
// }

}  // unnamed namespace

}  // namespace trento3d

int main(int argc, char* argv[]) {
  using namespace trento3d;

  // Parse options with boost::program_options.
  // There are quite a few options, so let's separate them into logical groups.
  using OptDesc = po::options_description;

  using VecStr = std::vector<std::string>;
  OptDesc main_opts{};
  main_opts.add_options()
    ("projectile", po::value<VecStr>()->required()->
     notifier(  // use a lambda to verify there are exactly two projectiles
         [](const VecStr& projectiles) {
           if (projectiles.size() != 2)
            throw po::required_option{"projectile"};
           }),
     "projectile symbols")
    ("number-events", po::value<int>()->default_value(1),
     "number of events");

  // Make all main arguments positional.
  po::positional_options_description positional_opts{};
  positional_opts
    .add("projectile", 2)
    .add("number-events", 1);

  using VecPath = std::vector<fs::path>;
  OptDesc general_opts{"general options"};
  general_opts.add_options()
    ("help,h", "show this help message and exit")
    ("version", "print version information and exit")
    ("bibtex", "print bibtex entry and exit")
    // ("default-config", "print a config file with default settings and exit")
    ("config-file,c", po::value<VecPath>()->value_name("FILE"),
     "configuration file\n(can be passed multiple times)");

  OptDesc output_opts{"output options"};
  output_opts.add_options()
    ("quiet,q", po::bool_switch(),
     "do not print event properties to stdout")
    ("output,o", po::value<fs::path>()->value_name("PATH"),
     "HDF5 file or directory for text files")
    ("no-header", po::bool_switch(),
     "do not write headers to text files")
    ("ncoll", po::bool_switch(),
     "calculate binary collisions")
    ("binary", po::bool_switch()->default_value(false, "FALSE"),
     "write condensed binary format to stdout")
    ("binary-sigma", po::bool_switch()->default_value(false, "FALSE"),
     "append weighted sigma-finder radius column to binary output");

  OptDesc phys_opts{"physical options"};
  phys_opts.add_options()
    // Random seed
    ("random-seed",
     po::value<int64_t>()->value_name("INT")->default_value(-1, "auto"),
     "random seed")
    // Nuclear configuration parameters
    ("form-width",
     po::value<double>()->value_name("FLOAT")->default_value(.5, "0.5"),
     "form-factor width of inelastic collisions [fm]")
    ("nucleon-width,w",
     po::value<double>()->value_name("FLOAT")->default_value(0.5, "0.5"),
     "Gaussian nucleon width [fm]")
    ("constit-width,v",
     po::value<double>()->value_name("FLOAT")->default_value(0.5, "same"),
     "Gaussian constituent width [fm]")
    ("constit-number,m",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1.0"),
     "number of constituents in the nucleon")
    ("nucleon-min-dist,d",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum nucleon-nucleon distance [fm]")
    // Cuts
    ("b-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum impact parameter [fm]")
    ("b-max",
     po::value<double>()->value_name("FLOAT")->default_value(-1., "auto"),
     "maximum impact parameter [fm]")
    ("npart-min",
     po::value<int>()->value_name("INT")->default_value(0, "0"),
     "minimum Npart cut")
    ("npart-max",
     po::value<int>()->value_name("INT")->default_value(
     std::numeric_limits<int>::max(), "INT_MAX"), "maximum Npart cut")
    ("mult-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum multiplicity cut")
    ("mult-max",
     po::value<double>()->value_name("FLOAT")->default_value(
     std::numeric_limits<double>::max(), "DOUBLE_MAX"), "maxmimum multiplicity cut")
    ("mult-etas-low",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "etas lower bound for calculating multiplicity")
    ("mult-etas-high",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "etas upper bound for calculating multiplicity")
    // Model parameters
    ("sqrts,s",
     po::value<double>()->value_name("FLOAT")->default_value(2760., "2760."),
     "CoM N-N energy of collision [GeV], determines cross-section, ybeam, etc.")
    ("fluctuation,k",
     po::value<double>()->value_name("FLOAT")->default_value(0.3, "0.3"),
     "fluctuation of energy sharing of central fireball and fragments")
    ("kT-min,t",
     po::value<double>()->value_name("FLOAT")->default_value(0.4, "0.4"),
     "etas_max ~ ln(sqrts/kT-min) [GeV]")
    ("reduced-thickness,p",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "reduced thickness parameter\n!!currently set to 0 internally!!")
    ("shape-alpha",
     po::value<double>()->value_name("FLOAT")->default_value(4.0, "4.0"),
     "fragmentation region a: [-ln(x)]^a x^b")
    ("shape-beta",
     po::value<double>()->value_name("FLOAT")->default_value(0.5, "0.5"),
     "fragmentation region b: [-ln(x)]^a x^b")
    ("mid-power",
     po::value<double>()->value_name("FLOAT")->default_value(0.4, "0.4"),
     "power of ETmid ~ #1*sqrts^#2")
    ("mid-norm",
     po::value<double>()->value_name("FLOAT")->default_value(0.3, "0.3"),
     "norm of ETmid ~ #1*sqrts^#2")
    ("overall-norm",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "multiplier applied to entire grid")
    ("flatness",
     po::value<double>()->value_name("FLOAT")->default_value(1.5, "1.5"),
     "flatness parameter of central plateau");

  OptDesc grid_opts{"grid options"};
  grid_opts.add_options()
    ("grid-max",
     po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
     "xy max [fm]\n(grid extends from -max to +max)")
    ("grid-step",
     po::value<double>()->value_name("FLOAT")->default_value(0.2, "0.2"),
     "step size [fm]")
    ("nsteps-etas",
     po::value<int>()->value_name("INT")->default_value(1, "1"),
     "number of grid steps in etas")
    ("etas-shift",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "etas at grid center,\nuse if Ebeam1 != Ebeam2");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  OptDesc usage_opts{};
  usage_opts
    .add(general_opts)
    .add(output_opts)
    .add(phys_opts)
    .add(grid_opts);

  // Now a meta-group containing _all_ options.
  OptDesc all_opts{};
  all_opts
    .add(usage_opts)
    .add(main_opts);

  // Will be used several times.
  const std::string usage_str{
    "usage: trento [options] projectile projectile [number-events = 1]\n"};

  try {
    // Initialize a VarMap (boost::program_options::variables_map).
    // It will contain all configuration values.
    VarMap var_map{};

    // Parse command line options.
    po::store(po::command_line_parser(argc, argv)
        .options(all_opts).positional(positional_opts).run(), var_map);

    // Handle options that imply immediate exit.
    // Must do this _before_ po::notify() since that can throw exceptions.
    if (var_map.count("help")) {
      std::cout
        << usage_str
        << "\n"
           "projectile = { p | d | Al | Cu | Cu2 | Xe | Xe2 | Au | Au2 | Pb | U | U2 | U3 }\n"
        << usage_opts
        << "\n"
           "see the online documentation for complete usage information\n";
      return 0;
    }
    if (var_map.count("version")) {
      print_version();
      return 0;
    }
    if (var_map.count("bibtex")) {
      print_bibtex();
      return 0;
    }
    // if (var_map.count("default-config")) {
    //   print_default_config();
    //   return 0;
    // }

    // Merge any config files.
    if (var_map.count("config-file")) {
      // Everything except general_opts.
      OptDesc config_file_opts{};
      config_file_opts
        .add(main_opts)
        .add(output_opts)
        .add(phys_opts)
        .add(grid_opts);

      for (const auto& path : var_map["config-file"].as<VecPath>()) {
        if (!fs::exists(path)) {
          throw po::error{
            "configuration file '" + path.string() + "' not found"};
        }
        fs::ifstream ifs{path};
        po::store(po::parse_config_file(ifs, config_file_opts), var_map);
      }
    }

    // Set default constituent width equal to the nucleon width.
    if (var_map["constit-width"].defaulted()) {
      var_map.at("constit-width").value() = var_map["nucleon-width"].as<double>();
    }

    double nucleon_width = var_map["nucleon-width"].as<double>();
    double constituent_width = var_map["constit-width"].as<double>();
    double constituent_number = var_map["constit-number"].as<double>();

    if (constituent_number < 1.)
      throw po::error{"must have at least one constituent"};

    // Constituent and nucleon widths must be non-negative.
    if ((nucleon_width <= 0.) || (constituent_width <= 0.))
      throw po::error{"nucleon and constituent widths must be positive"};

    // Constituent width cannot be larger than nucleon width.
    if (constituent_width > nucleon_width)
      throw po::error{"constituent width cannot be larger than nucleon width"};

    // Cannot fit nucleon width using single constituent if different sizes.
    if ((constituent_width < nucleon_width) && constituent_number <= 1.)
      throw po::error{"cannot fit nucleon width using single constituent if different sizes"};

    if (var_map["nsteps-etas"].as<int>() <= 0)
      throw po::error{"nsteps-etas must be positive"};

    // Save all the final values into var_map.
    // Exceptions may occur here.
    po::notify(var_map);

    // Go!
    Collider collider{var_map};
    collider.run_events();
  }
  catch (const po::required_option&) {
    // Handle this exception separately from others.
    // This occurs e.g. when the program is executed with no arguments.
    std::cerr << usage_str << "run 'trento-3 --help' for more information\n";
    return 1;
  }
  catch (const std::exception& e) {
    // For all other exceptions just output the error message.
    std::cerr << e.what() << '\n';
    return 1;
  }

  return 0;
}
