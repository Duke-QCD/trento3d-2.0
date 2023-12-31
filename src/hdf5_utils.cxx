// TRENTO-3D: Reduced Thickness Event-by-event Nuclear Topology
// Copyright (c) 2015-2023 Jonah E. Bernhard, J. Scott Moreland, Steffen A. Bass, Weiyao Ke, Derek Soeder
// MIT License

#include "hdf5_utils.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>  // include this only after explicitly including the std functionality we need

namespace trento3d {

namespace hdf5 {

bool filename_is_hdf5(const fs::path& path) {
  if (!path.has_extension())
    return false;

  auto hdf5_exts = {".hdf5", ".hdf", ".hd5", ".h5"};
  auto result = std::find(hdf5_exts.begin(), hdf5_exts.end(), path.extension());

  return result != hdf5_exts.end();
}

bool filename_is_hdf5(const std::string& path) {
  return filename_is_hdf5(fs::path{path});
}

#ifdef TRENTO3D_HDF5

H5::H5File try_open_file(const std::string& path, unsigned int flags) {
  if (!fs::exists(path))
    throw std::invalid_argument{"file '" + path + "' does not exist"};

  if (!H5::H5File::isHdf5(path))
    throw std::invalid_argument{"'" + path + "' is not a valid HDF5 file"};

  return H5::H5File{path, flags};
}

#endif  // TRENTO3D_HDF5

}  // namespace hdf5

}  // namespace trento3d
