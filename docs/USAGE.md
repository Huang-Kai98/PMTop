# PMTop Usage Guide

## 1. Build environment

The root `CMakeLists.txt` selects C++17 and Intel `icpx`. It requires the CMake packages `MKL`, `arpackng`, `OpenMP`, `nlohmann_json` (at least 3.11.2), and `spdlog`. oneMKL is configured for dynamic linking, Intel threading, and the LP64 interface.

Load the compiler and library environment, then run the build commands in the [README](../README.md):

```sh
source /opt/intel/oneapi/setvars.sh
```

This is the oneAPI installation path on the development machine. Adjust it if your installation uses a different location. Keep the environment active when running the executable so its shared libraries can be found.

The `test/` directory contains standalone programs. The current top-level CMake configuration does not add that directory or register a CTest suite. A successful application build therefore does not imply that those programs have been compiled or executed.

## 2. Configuration and working directory

`src/main/main.cpp` loads the literal filename `Config.json`. Dataset filenames are resolved relative to the process working directory. There is currently no command-line argument for selecting a configuration file.

The minimal configuration for the supplied full-element example is:

```json
{
  "filename": "TopShellDamp1010",
  "read_mode": "all_element",
  "description": "TopShellDamp1010 damped modal sensitivity example"
}
```

A copy is available at [examples/Config.TopShellDamp1010.json](examples/Config.TopShellDamp1010.json). Copy it to `Config.json` in a separate run directory, together with the corresponding model files. This keeps generated output separate from the published data.

### Reader modes

| `read_mode` | Input behavior | Additional configuration |
| --- | --- | --- |
| `all_element` | Read each element's matrices from files sharing the `filename` prefix | None beyond `filename` and the driver's `description` |
| `single_element` | Reuse a matrix entry from `element_library.txt` across the mesh | `element_tag` |
| `single_element_two_coeff` | Reuse stiffness/mass matrices and construct damping from two Rayleigh coefficient pairs | `element_tag`, `RayleighDampingCoef.Coef1`, `RayleighDampingCoef.Coef2`, and `Interface` |

Both single-element modes require `element_library.txt`. A compatible library and complete `Shell_2010` inputs are included in the [minimal shell example](../examples/minimal-shell/README.md). Use `all_element` for the supplied TopShellDamp datasets.

In `single_element_two_coeff` mode, damping is assembled as `C = alpha * M + beta * K`. The implemented interface type is `OneWidth`: elements with centroid `y <= Interface.Interface` use `Coef1`; the remaining elements use `Coef2`.

The original configuration includes `IndexType`, but the active driver selects `DMAC` directly in the source. Changing that field alone does not switch the active sensitivity calculation. Alternative `DMIR` and `DMPC` calls are currently commented out in `src/main/main.cpp`.

## 3. TopShellDamp files

| Prefix | File count | Included content |
| --- | --- | --- |
| `TopShellDamp1010` | 10 | Model lists, stiffness/mass/damping matrices, combined matrix export, eigenvalue/eigenvector files |
| `TopShellDamp4040` | 10 | Model lists, stiffness/mass/damping matrices, combined matrix export, eigenvalue/eigenvector files |
| `TopShellDamp6040` | 8 | Model lists, stiffness/mass/damping matrices, combined matrix export |

For a model prefix `<name>`, the full-element reader uses:

| Suffix | Purpose |
| --- | --- |
| `_NLIST.lis` | Node coordinates and related node data |
| `_ELIST.lis` | Element connectivity and attributes |
| `_ETLIST.lis` | Element type definitions |
| `_DLIST.lis` | Prescribed degrees of freedom |
| `_elem_k.txt` | Element stiffness matrices |
| `_elem_m.txt` | Element mass matrices |
| `_elem_c.txt` | Element damping matrices |

The reader also checks for `_CELIST.lis` and `_CPLIST.lis` constraint files. Their absence disables those optional constraint mechanisms. They are not part of the uploaded TopShellDamp sets.

The `_ElementMatrix.txt`, `.Eigenvalue.txt`, and `.Eigenvector.txt` files are retained as accompanying data; the active `all_element` initialization does not read them. The driver computes its modal solution from the supplied matrices.

To select another uploaded model, change `filename` in your run directory's `Config.json` to `TopShellDamp4040` or `TopShellDamp6040` and copy that model's `.lis` and `_elem_*.txt` files into the same directory. Keep `read_mode` set to `all_element`.

## 4. Current experiment

The current driver:

1. Forces single-thread settings for its timing comparisons.
2. Loads the model and exports its geometry to VTK.
3. Initializes element design variables to 1.0.
4. Requests 10 modes for the right quadratic eigenproblem.
5. Selects mode index 0 and constructs perturbed modal vectors.
6. Compares direct and adjoint sensitivity formulations, reporting CPU time, difference norms, and representative element results.

Mode count, selected mode, repetition count, interpolation functions, and sensitivity target are source-level settings. The noise generator currently uses seed 1024 and a default noise level of 0.2. These settings describe the checked-in experiment; they are not a claim of validation for every dataset.

## 5. Output

- `log.txt`: diagnostic logging in the run directory.
- `<filename>.vtk`: model geometry exported during initialization.
- Standard output: model sizes, modal indicators, sensitivity comparisons, and CPU timings.

The `saveElementValue(...)` calls at the end of the driver are commented out. Enable the relevant calls and rebuild if you want those additional sensitivity exports.

## 6. Troubleshooting and validation scope

| Symptom | Check |
| --- | --- |
| `icpx` cannot be found | Load the Intel oneAPI compiler environment before configuring |
| A CMake package cannot be found | Install its development package and expose its package configuration directory |
| `Config.json` cannot be opened | Run from the directory containing the intended configuration |
| Model files cannot be opened | Match `filename` exactly, including case, and copy the required files into the run directory |
| `element_library.txt` is missing | Use the complete minimal shell example, or provide a compatible library for your single-element model |
| A shared library cannot be loaded | Restore the compiler/library runtime environment used for the build |

### Verified build

After sourcing `/opt/intel/oneapi/setvars.sh`, Release configuration and compilation completed successfully with Intel oneAPI C++ Compiler 2025.3.0, oneMKL 2025.3.0, and CMake 3.28.3. The resulting executable is `build/src/main/PMTop`.

Compiler warnings remain, including a reference to a temporary in `src/NumericalAlgebra/sparsemat.cpp`, intentional null-pointer dereferences in `src/General/error.cpp`, and ambiguous `else` statements. Successful compilation does not establish runtime correctness.

### Configuration-specific runtime notes

The [minimal shell example](../examples/minimal-shell/README.md) was successfully compiled and run on 2026-10-13 with the published C++ driver. Its `single_element_two_coeff` configuration loaded 231 nodes and 200 elements (1,386 DOFs), completed all five sensitivity methods and produced a geometry export. The full-vector ADNE/ADBM L2 difference was approximately `1.51e-13`. The supplied launcher checks the input hashes and actual output and writes `validation.json`. See the example README and its recorded result for the unchanged ASPM stopping settings and accuracy limits.

PMTop has multiple input paths. `AnsysLoad` selects `all_element`, `single_element`, or `single_element_two_coeff` according to the `read_mode` in the working directory's `Config.json`; these modes are **not interchangeable**. Before running, ensure that `filename` points to the correct ANSYS model lists and that the required matrix files or `element_library.txt` are present. The repository's `data/Config.json` targets the locally supplied `Shell_12060`/`element_library.txt` workflow, whose source files are not committed; the `docs/examples/Config.TopShellDamp1010.json` file instead selects the included `TopShellDamp1010` full-element workflow.

An earlier check of **only** the `TopShellDamp1010` `all_element` example read 121 nodes and 300 elements (726 DOFs) but stopped during element-matrix loading. In the published `src/FEM/AnsysLoad.cpp`, `ReadAnasysElementStifAndElementMass()` declares `elem` and `Index` inside each line-reading loop; this can leave the matrix-value lines without a valid retained element identifier or index. That observation is specific to this reader path, **not** a conclusion that other correctly configured PMTop workflows cannot run. Selecting a `single_element*` workflow avoids this function, but requires a compatible element library and model data.

The successful minimal run applies to its exact inputs and reader mode. The repository does not currently include an independently verified end-to-end run of all paper examples. No C++ source code was changed for the minimal example.

## 7. Attribution

See [NOTICE.md](../NOTICE.md) for project and third-party attribution and [LICENSE](../LICENSE) for the LGPL-2.1 license text.
