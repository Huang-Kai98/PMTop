# PMTop Usage Guide

## 1. Build environment

The root `CMakeLists.txt` selects C++17 and Intel `icpx`. It requires the CMake packages `MKL`, `arpackng`, `OpenMP`, `nlohmann_json` (at least 3.11.2), and `spdlog`. oneMKL is configured for dynamic linking, Intel threading, and the LP64 interface.

Load the compiler and library environment for your installation, then run the build commands in the [README](../README.md). Keep that environment active when running the executable so its shared libraries can be found.

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

Both single-element modes require `element_library.txt`, which is not included. Use `all_element` for the supplied TopShellDamp datasets.

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
| `element_library.txt` is missing | Use the supplied `all_element` example, or provide your own element library for a single-element mode |
| A shared library cannot be loaded | Restore the compiler/library runtime environment used for the build |

This documentation update does not change the solver or certify numerical results. Build and runtime verification remain outstanding because Intel `icpx` was not available on PATH during preparation. Check results against an appropriate reference solution for your model.

## 7. Attribution

See [NOTICE.md](../NOTICE.md) for project and third-party attribution and [LICENSE](../LICENSE) for the LGPL-2.1 license text.
