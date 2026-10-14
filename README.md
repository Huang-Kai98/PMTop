# PMTop

PMTop is a C++ research codebase for finite element analysis, damped modal analysis, and eigenvalue/eigenvector sensitivity calculations. It reads ANSYS-exported model data for structural dynamics and sensitivity studies.

The current executable compares sensitivity methods for a damped system. It is a research driver with experiment settings in the source, rather than a command-line application with selectable workflows.

## Capabilities

- Import ANSYS node, element, element-type, boundary-condition, and element-matrix files.
- Assemble stiffness, mass, and damping matrices.
- Solve modal problems using ARPACK-NG and Intel oneMKL/PARDISO.
- Compare Nelson, augmented-system, and adjoint sensitivity calculations.
- Evaluate modal indicators including MAC, MIR, and MPC.
- Export model geometry to VTK and record diagnostic output.
- Include standalone development and verification programs.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/Algorithm/` | Modal analysis and sensitivity calculations |
| `src/FEM/` | ANSYS readers, finite element types, assembly, and boundary conditions |
| `src/NumericalAlgebra/` | Matrix/vector types, linear solvers, and eigensolvers |
| `src/General/` | Configuration, logging, containers, and utilities |
| `src/Mesh/` | VTK-related utilities |
| `src/main/` | Current experiment driver |
| `test/` | Standalone development and verification programs |
| `data/` | Configuration and selected TopShellDamp datasets |
| `examples/minimal-shell/` | Complete, verified small shell-model inputs and launcher |
| `docs/USAGE.md` | Build, configuration, dataset, and execution details |

## Minimal verified example

The [minimal shell example](examples/minimal-shell/README.md) includes all inputs for a 231-node, 200-element model. It runs the existing C++ executable and completes all five sensitivity methods on a single CPU thread. After installing the build dependencies below:

```sh
source /opt/intel/oneapi/setvars.sh
bash examples/minimal-shell/run.sh
```

The script builds PMTop, creates a separate output directory, runs the model and checks its output. A successful run prints `PASS` and saves `validation.json`, `stdout.log` and `Shell_2010.vtk`. The example README records the settings and numerical comparison; the C++ algorithms are unchanged.

## Build

The current build targets a Linux environment with Intel's C++ compiler. Install and expose the following dependencies before configuring:

- CMake 3.10 or newer, as declared by the project.
- Intel oneAPI C/C++ compilers (`icx` and `icpx`), with C++17 support.
- Intel oneMKL, including PARDISO and its CMake package.
- ARPACK-NG, including C++ headers and the `arpackng` CMake package.
- OpenMP, nlohmann/json 3.11.2 or newer, and spdlog.

Load the Intel oneAPI environment before configuring and keep it active when running PMTop:

```sh
source /opt/intel/oneapi/setvars.sh
git clone https://github.com/Huang-Kai98/PMTop.git
cd PMTop
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx
cmake --build build --parallel
```

The expected executable location is `build/src/main/PMTop`. If dependencies are installed outside standard search paths, set `CMAKE_PREFIX_PATH` or the relevant package's CMake directory. The release flags include `-xHost`, so the resulting binary targets the build machine's CPU capabilities.

## Configure a model run

PMTop selects its model input workflow from the `Config.json` file in the **run directory** (the current working directory). The executable supports three distinct `read_mode` values:

- `all_element`: reads individual element stiffness, mass, and damping matrices from `<filename>_elem_k.txt`, `<filename>_elem_m.txt`, and `<filename>_elem_c.txt`, together with the ANSYS model lists.
- `single_element`: uses an entry selected by `element_tag` in a local `element_library.txt` for element matrices.
- `single_element_two_coeff`: uses the element library for stiffness and mass and generates two-region Rayleigh damping from the coefficients and `Interface` specified in `Config.json`.

Select the input configuration that matches the model files in the run directory. The supplied [minimal shell example](examples/minimal-shell/README.md) includes a compatible configuration and all required inputs. For the `TopShellDamp*` datasets, use `read_mode: "all_element"` and the example configurations in `docs/examples/`. The default `data/Config.json` targets a separate `Shell_12060` study whose model files are not part of this repository.

### Example: supplied full-element dataset

For the supplied `TopShellDamp1010` model, `docs/examples/Config.TopShellDamp1010.json` selects `all_element`. The commands below show the file layout and invocation:

```sh
mkdir -p build/run-topshelldamp1010
cp data/TopShellDamp1010*.lis data/TopShellDamp1010_elem_*.txt \
  build/run-topshelldamp1010/
cp docs/examples/Config.TopShellDamp1010.json \
  build/run-topshelldamp1010/Config.json
(cd build/run-topshelldamp1010 && ../src/main/PMTop)
```

The program reads `Config.json` and writes `log.txt` in its working directory. After successful model initialization, it is designed to export a model VTK file and print sensitivity comparisons and CPU timings to the terminal. See [the usage guide](docs/USAGE.md) for configuration fields and the other datasets.

Release configuration and compilation were verified using Intel oneAPI C++ Compiler 2025.3.0 and oneMKL 2025.3.0. Run PMTop with a matching model dataset and `Config.json`; see the [usage guide](docs/USAGE.md) and [minimal shell example](examples/minimal-shell/README.md).

## Included data

The repository contains 28 files for `TopShellDamp1010`, `TopShellDamp4040`, and `TopShellDamp6040`, totaling approximately 256 MB. Other local datasets, build outputs, IDE settings, and the original development history are excluded. The published history starts with a snapshot of the current source.

## License and attribution

PMTop's original contributions are copyright (c) 2026 Huang Kai and are distributed under the **GNU Lesser General Public License, version 2.1 only (LGPL-2.1-only)**. See [LICENSE](LICENSE) for the license text and [COPYING](COPYING) for the GNU GPL version 2 referenced by it.

Existing third-party notices and terms remain applicable to their respective code. In particular, the MFEM notice in `src/NumericalAlgebra/gmres.cpp` and the VTK notice in `src/Mesh/vtk.h` are preserved. See [NOTICE.md](NOTICE.md) for attribution details. External dependencies retain their own licenses.

Maintainer: [Huang Kai](https://github.com/Huang-Kai98).
