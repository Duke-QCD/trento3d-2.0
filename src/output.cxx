// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#include "output.h"

#include <stdio.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>

#include "event.h"
#include "hdf5_utils.h"

namespace trento3d {

namespace {

// These output functions are invoked by the Output class.

void write_stream(std::ostream& os, int width, int num, double impact_param,
    int ncoll, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  // Write a nicely-formatted line of event properties.
  os << setprecision(10)
     << setw(width)              << num
     << ' ' << setw(14) << fixed << impact_param
     << ' ' << setw(4)           << event.npart()
     << ' ' << setw(4)           << event.npartA()
     << ' ' << setw(4)           << event.npartB();

  // Output ncoll if calculated
  if (ncoll > 0) os << ' ' << setw(5) << ncoll;

  os << ' ' << setw(17) << scientific << event.multiplicity();

  const auto & detdetas = event.dET_detas();
  const auto & eccang   = event.ecc_ang();

  for (size_t i = 0; i < detdetas.size(); ++i) {
    os << scientific << ' ' << setw(17) << detdetas[i] << fixed;
    for (const auto & itecc : event.ecc_mag()) {
      os << ' ' << setw(13) << itecc.second[i] << ' ' << setw(13) << eccang.at(itecc.first)[i];
    }
  }

  os << '\n';
}

void write_stream_binary(const Event& event, float * buf) {
  const auto & detdetas = event.dET_detas();

  for (size_t i = 0; i < detdetas.size(); i++)
    buf[i] = static_cast<float>(detdetas[i]);

  fwrite(buf, sizeof(buf[0]), detdetas.size(), stdout);
}

void write_stream_binary_sigma(const Event& event) {
  // note: write_stream_binary will always be called before write_stream_binary_sigma (see Output::Output)

  const auto & detdetas = event.dET_detas();
  float sigma = 0., totwt = 0.;

  for (size_t iz = 0; iz < detdetas.size(); ++iz) {
    auto wt = detdetas.at(iz);
    sigma += wt * std::get<0>(event.get_energy_frac_radii(iz)) * event.dxy();
    totwt += wt;
  }
  sigma /= totwt;

  fwrite(&sigma, sizeof(sigma), 1, stdout);
}

void write_text_file(const fs::path& output_dir, int width, int num,
    double impact_param, int ncoll, const Event& event, bool header) {
  // Open a numbered file in the output directory.
  // Pad the filename with zeros.
  std::ostringstream padded_fname{};
  padded_fname << std::setw(width) << std::setfill('0') << num << ".dat";
  fs::ofstream ofs{output_dir / padded_fname.str()};

  if (header) {
    // Write a commented header of event properties as key = value pairs.
    ofs << std::setprecision(10)
        << "# event "   << num                  << '\n'
        << "# b     = " << impact_param         << '\n'
        << "# npart = " << event.npart()        << '\n';

    // Output ncoll if calculated
    if (ncoll > 0) ofs << "# ncoll = " << ncoll << '\n';

    ofs << "# mult  = " << event.multiplicity() << '\n';

    //for (const auto& ecc : event.eccentricity())
    //  ofs << "# e" << ecc.first << "    = " << ecc.second << '\n';
  }

  // Write IC profile as a block grid.  Use C++ default float format (not
  // fixed-width) so that trailing zeros are omitted.  This significantly
  // increases output speed and saves disk space since many grid elements are
  // zero.
  for (const auto & slice : event.Density()) {
    for (const auto & row : slice){
      auto&& iter = row.begin();
      // Write all row elements except the last with a space delimiter afterwards.
      do {
        ofs << *iter << ' ';
      } while (++iter != --row.end());

      // Write the last element and a linebreak.
      ofs << *iter << '\n';
    }
  }
}

#ifdef TRENTO3D_HDF5

/// Simple functor to write many events to an HDF5 file.
class HDF5Writer {
 public:
  /// Prepare an HDF5 file for writing.
  HDF5Writer(const fs::path& filename);

  /// Write an event.
  void operator()(int num, double impact_param,
      int ncoll, const Event& event) const;

 private:
  /// Internal storage of the file object.
  H5::H5File file_;
};

// Add a simple scalar attribute to an HDF5 dataset.
template <typename T>
void hdf5_add_scalar_attr(
    const H5::DataSet& dataset, const std::string& name, const T& value) {
  const auto& datatype = hdf5::type<T>();
  auto attr = dataset.createAttribute(name, datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

HDF5Writer::HDF5Writer(const fs::path& filename)
    : file_(filename.string(), H5F_ACC_TRUNC)
{}

void HDF5Writer::operator()(int num, double impact_param,
    int ncoll, const Event& event) const {
  // Prepare arguments for new HDF5 dataset.

  // The dataset name is a prefix plus the event number.
  const std::string name{"event_" + std::to_string(num)};

  // Cache a reference to the event grid -- will need it several times.
  const auto& grid = event.Density();

  // Define HDF5 datatype and dataspace to match the grid.
  const auto& datatype = hdf5::type<Event::Grid::element>();
  std::array<hsize_t, Event::Grid3D::dimensionality> shape;
  std::copy(grid.shape(), grid.shape() + shape.size(), shape.begin());
  auto dataspace = hdf5::make_dataspace(shape);

  // Set dataset storage properties.
  H5::DSetCreatPropList proplist{};
  // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
  // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
  // sense to chunk this way, since chunks must be read contiguously and there's
  // no reason to read a partial grid.
  proplist.setChunk(shape.size(), shape.data());
  // Set gzip compression level.  4 is the default in h5py.
  proplist.setDeflate(4);

  // Create the new dataset and write the grid.
  auto dataset = file_.createDataSet(name, datatype, dataspace, proplist);
  dataset.write(grid.data(), datatype);

  // Write event attributes.
  hdf5_add_scalar_attr(dataset, "b", impact_param);
  hdf5_add_scalar_attr(dataset, "npart", event.npart());

  // Write ncoll if calculated
  if (ncoll > 0) hdf5_add_scalar_attr(dataset, "ncoll", ncoll);

  hdf5_add_scalar_attr(dataset, "mult", event.multiplicity());
  //for (const auto& ecc : event.eccentricity())
  //  hdf5_add_scalar_attr(dataset, "e" + std::to_string(ecc.first), ecc.second);
}

#endif  // TRENTO3D_HDF5

}  // unnamed namespace

Output::Output(const VarMap& var_map) {
  // Determine the required width (padding) of the event number.  For example if
  // there are 10 events, the numbers are 0-9 and so no padding is necessary.
  // However given 11 events, the numbers are 00-10 with padded 00, 01, ...
  auto nevents = var_map["number-events"].as<int>();
  auto width = static_cast<int>(std::ceil(std::log10(nevents)));

  // Write to stdout unless the quiet option was specified.
  if (!var_map["quiet"].as<bool>()) {
    if (var_map["binary"].as<bool>()) {
      buf_.resize(static_cast<size_t>(var_map["nsteps-etas"].as<int>()));  // important: write_stream_binary accesses buf_'s storage directly and without safety checks, so it must be properly allocated

      writers_.emplace_back(
        [this](int num, double impact_param, int ncoll, const Event& event) {
          write_stream_binary(event, buf_.data());
        }
      );
      if (var_map["binary-sigma"].as<bool>()) {
        writers_.emplace_back(
          [this](int num, double impact_param, int ncoll, const Event& event) {
            write_stream_binary_sigma(event);
          }
        );
      }
    } else {
      writers_.emplace_back(
        [width](int num, double impact_param, int ncoll, const Event& event) {
          write_stream(std::cout, width, num, impact_param, ncoll, event);
        }
      );
    }
  }

  // Possibly write to text or HDF5 files.
  if (var_map.count("output")) {
    const auto& output_path = var_map["output"].as<fs::path>();
    if (hdf5::filename_is_hdf5(output_path)) {
#ifdef TRENTO3D_HDF5
      if (fs::exists(output_path) && !fs::is_empty(output_path))
        throw std::runtime_error{"file '" + output_path.string() +
                                 "' exists, will not overwrite"};
      writers_.emplace_back(HDF5Writer{output_path});
#else
      throw std::runtime_error{"HDF5 output was not compiled"};
#endif  // TRENTO3D_HDF5
    } else {
      // Text files are all written into a single directory.  Require the
      // directory to be initially empty to avoid accidental overwriting and/or
      // spewing files into an already-used location.  If the directory does not
      // exist, create it.
      if (fs::exists(output_path)) {
        if (!fs::is_empty(output_path)) {
          throw std::runtime_error{"output directory '" + output_path.string() +
                                   "' must be empty"};
        }
      } else {
        fs::create_directories(output_path);
      }
      auto header = !var_map["no-header"].as<bool>();
      writers_.emplace_back(
        [output_path, width, header](int num, double impact_param,
          int ncoll, const Event& event) {
          write_text_file(output_path, width, num,
              impact_param, ncoll, event, header);
        }
      );
    }
  }
}

}  // namespace trento3d
