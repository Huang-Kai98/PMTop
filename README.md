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

**Matching the configuration to its input files is essential.** The committed `data/Config.json` describes a `Shell_12060` two-region study, but its model lists are not included in the public repository. It is not the configuration for the included `TopShellDamp*` files or the minimal example. For a complete supplied single-element workflow, use `examples/minimal-shell/run.sh`, which selects `Shell_2010` and its compatible element library. The previous failed `all_element` example test must not be taken to mean that all PMTop configurations are unable to run.

### Example: supplied full-element dataset

For the uploaded `TopShellDamp1010` model, `docs/examples/Config.TopShellDamp1010.json` selects `all_element`. The commands below show the required file layout and invocation; this particular example has not been independently confirmed to finish with the current published source. See [configuration-specific runtime notes](docs/USAGE.md#6-troubleshooting-and-validation-scope).

```sh
mkdir -p build/run-topshelldamp1010
cp data/TopShellDamp1010*.lis data/TopShellDamp1010_elem_*.txt \
  build/run-topshelldamp1010/
cp docs/examples/Config.TopShellDamp1010.json \
  build/run-topshelldamp1010/Config.json
(cd build/run-topshelldamp1010 && ../src/main/PMTop)
```

The program reads `Config.json` and writes `log.txt` in its working directory. After successful model initialization, it is designed to export a model VTK file and print sensitivity comparisons and CPU timings to the terminal. See [the usage guide](docs/USAGE.md) for configuration fields and the other datasets.

Release configuration and compilation were verified on the development machine after sourcing `/opt/intel/oneapi/setvars.sh`, using Intel oneAPI C++ Compiler 2025.3.0 and oneMKL 2025.3.0. The build completed with warnings. Execution requires a dataset and `Config.json` compatible with the selected reader mode. No independently verified end-to-end reproduction log for the paper's numerical results is currently included in this repository.

## Included data

The repository contains 28 files for `TopShellDamp1010`, `TopShellDamp4040`, and `TopShellDamp6040`, totaling approximately 256 MB. Other local datasets, build outputs, IDE settings, and the original development history are excluded. The published history starts with a snapshot of the current source.

## License and attribution

PMTop's original contributions are copyright (c) 2026 Huang Kai and are distributed under the **GNU Lesser General Public License, version 2.1 only (LGPL-2.1-only)**. See [LICENSE](LICENSE) for the license text and [COPYING](COPYING) for the GNU GPL version 2 referenced by it.

Existing third-party notices and terms remain applicable to their respective code. In particular, the MFEM notice in `src/NumericalAlgebra/gmres.cpp` and the VTK notice in `src/Mesh/vtk.h` are preserved. See [NOTICE.md](NOTICE.md) for attribution details. External dependencies retain their own licenses.

Maintainer: [Huang Kai](https://github.com/Huang-Kai98).
