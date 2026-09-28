# PMTop

PMTop is a C++ research codebase for finite element analysis, damped modal analysis, and eigenvalue/eigenvector sensitivity calculations. It reads ANSYS-exported model data and includes numerical routines for topology optimization.

The current executable compares sensitivity methods for a damped system. It is a research driver with experiment settings in the source, rather than a command-line application with selectable workflows.

## Capabilities

- Import ANSYS node, element, element-type, boundary-condition, and element-matrix files.
- Assemble stiffness, mass, and damping matrices.
- Solve modal problems using ARPACK-NG and Intel oneMKL/PARDISO.
- Compare Nelson, augmented-system, and adjoint sensitivity calculations.
- Evaluate modal indicators including MAC, MIR, and MPC.
- Export model geometry to VTK and record diagnostic output.
- Provide MMA/GCMMA optimization routines and standalone test programs.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/Algorithm/` | Modal analysis and sensitivity calculations |
| `src/FEM/` | ANSYS readers, finite element types, assembly, and boundary conditions |
| `src/NumericalAlgebra/` | Matrix/vector types, linear solvers, eigensolvers, and optimizers |
| `src/General/` | Configuration, logging, containers, and utilities |
| `src/Mesh/` | VTK-related utilities |
| `src/main/` | Current experiment driver |
| `test/` | Standalone development and verification programs |
| `data/` | Configuration and selected TopShellDamp datasets |
| `docs/USAGE.md` | Build, configuration, dataset, and execution details |

## Build

The current build targets a Linux environment with Intel's C++ compiler. Install and expose the following dependencies before configuring:

- CMake 3.10 or newer, as declared by the project.
- Intel oneAPI C/C++ compilers (`icx` and `icpx`), with C++17 support.
- Intel oneMKL, including PARDISO and its CMake package.
- ARPACK-NG, including C++ headers and the `arpackng` CMake package.
- OpenMP, nlohmann/json 3.11.2 or newer, and spdlog.

After loading your oneAPI environment:

```sh
git clone https://github.com/Huang-Kai98/PMTop.git
cd PMTop
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx
cmake --build build --parallel
```

The expected executable location is `build/src/main/PMTop`. If dependencies are installed outside standard search paths, set `CMAKE_PREFIX_PATH` or the relevant package's CMake directory. The release flags include `-xHost`, so the resulting binary targets the build machine's CPU capabilities.

## Run a supplied dataset

The committed `data/Config.json` retains the original experiment configuration. It references `Shell_12060` and `element_library.txt`, which are not included in this repository. Use the supplied full-element example for the uploaded TopShellDamp data:

```sh
mkdir -p build/run-topshelldamp1010
cp data/TopShellDamp1010*.lis data/TopShellDamp1010_elem_*.txt \
  build/run-topshelldamp1010/
cp docs/examples/Config.TopShellDamp1010.json \
  build/run-topshelldamp1010/Config.json
(cd build/run-topshelldamp1010 && ../src/main/PMTop)
```

The program reads `Config.json` from its working directory. It writes `log.txt` and a model VTK file there and prints sensitivity comparisons and CPU timings to the terminal. See [the usage guide](docs/USAGE.md) for configuration fields and the other datasets.

These instructions are derived from the current source and CMake configuration. A successful build and numerical run have not been verified as part of this documentation update; `icpx` was unavailable on the documentation environment's PATH.

## Included data

The repository contains 28 files for `TopShellDamp1010`, `TopShellDamp4040`, and `TopShellDamp6040`, totaling approximately 256 MB. Other local datasets, build outputs, IDE settings, and the original development history are excluded. The published history starts with a snapshot of the current source.

## License and attribution

PMTop's original contributions are copyright (c) 2026 Huang Kai and are distributed under the **GNU Lesser General Public License, version 2.1 only (LGPL-2.1-only)**. See [LICENSE](LICENSE) for the license text and [COPYING](COPYING) for the GNU GPL version 2 referenced by it.

Existing third-party notices and terms remain applicable to their respective code. In particular, the MFEM notice in `src/NumericalAlgebra/gmres.cpp` and the VTK notice in `src/Mesh/vtk.h` are preserved. See [NOTICE.md](NOTICE.md) for attribution details. External dependencies retain their own licenses.

Maintainer: [Huang Kai](https://github.com/Huang-Kai98).
