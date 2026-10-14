# Minimal PMTop shell example

This is a complete, small input set for the existing C++ `PMTop` executable: a 20-by-10 shell mesh, 231 nodes, 200 SHELL181 elements and 1,386 retained degrees of freedom. It runs finite element assembly, the quadratic eigenproblem and all five sensitivity methods in `src/main/main.cpp`:

- Traditional Nelson and traditional bordered methods.
- ADNE and ADBM.
- ASPM using the published driver's existing settings.

The input files total approximately 222 KiB. No files from the author's private data directory are required. The example uses the existing C++ algorithms; Python only checks the resulting log.

## Requirements

Linux, Bash, Python 3 (standard library only), `taskset`, and the PMTop build dependencies listed in the root README: Intel oneAPI `icx`/`icpx`, oneMKL, CMake, ARPACK-NG with its C++ headers/CMake package, OpenMP, nlohmann/json and spdlog. MATLAB is not required.

## Build and run

From a fresh clone:

```bash
git clone https://github.com/Huang-Kai98/PMTop.git
cd PMTop
source /opt/intel/oneapi/setvars.sh
bash examples/minimal-shell/run.sh
```

The script builds `build/src/main/PMTop`, copies only this example's six inputs to a new run directory, binds execution to one available logical CPU, and sets MKL/OpenMP thread counts to one. Every invocation uses a new run directory. An optional output directory must not already exist:

```bash
bash examples/minimal-shell/run.sh /tmp/PMTop_minimal_review_new
```

Success ends with:

```text
PASS: actual PMTop run; 231 nodes, 200 elements, 1386 DOFs; all five methods completed.
```

The output directory contains `stdout.log`, `log.txt`, `Shell_2010.vtk`, `validation.json` and the copied inputs. `validation.json` records the input checks, representative sensitivities and source hashes. Build failures or a nonzero solver exit stop the script.

To check an existing run without solving again:

```bash
python3 examples/minimal-shell/validate.py /path/to/run-directory
```

## Model and settings

`Config.json` selects `single_element_two_coeff`, `Shell_2010` and element-library tag 2. `element_library.txt` supplies the element stiffness and mass matrices. Rayleigh coefficients are `(alpha, beta) = (15, 1.5e-5)` and `(400, 3.5e-5)` in the two regions; the `OneWidth` interface parameter is 0.5. The four `.lis` files define nodes, connectivity, element types and constraints.

The existing driver selects ARPACK column 0, adds its existing deterministic observation noise and differentiates MAC. Its stiffness/mass/damping interpolation rules and sensitivity parameter are unchanged; this is not the later regional damage-identification workflow.

ASPM retains the published driver's absolute iterate-change tolerance `1e-5` and maximum of five iterations. These are runtime-example settings, not the later manuscript rerun settings. A successful run does not imply a relative equation residual or sensitivity error of `1e-8`.

## Verified results

The example was compiled from a clean checkout and run on 2026-10-13 using Intel oneAPI C++ 2025.3 and oneMKL 2025.3 on Linux. See [verified-run.json](verified-run.json) for the actual recorded result.

Representative zero-based element index: 199. Typical printed sensitivities:

| Method | MAC sensitivity |
|---|---:|
| Nelson | -0.00377632 |
| Bordered | -0.00377632 |
| ADBM | -0.00377632 |
| ADNE | -0.00377632 |
| ASPM | -0.00377071 |

The full-vector ADNE/ADBM L2 difference was approximately `1.51e-13`. The ASPM representative scalar differs by approximately `1.49e-3` relative under the unchanged loose stopping settings. The checker records this difference rather than claiming all methods have identical accuracy. Printed digits can vary across compatible compiler/library versions.

This example demonstrates the published driver using the supplied `single_element_two_coeff` configuration. Other model datasets use the input modes documented in [docs/USAGE.md](../../docs/USAGE.md).

## License

Copyright (c) 2026 Huang Kai. LGPL-2.1-only; see the repository's `LICENSE`, `COPYING` and `NOTICE.md`. Existing third-party notices are unchanged.
