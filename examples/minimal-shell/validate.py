#!/usr/bin/env python3
"""Check this small model's actual PMTop output; uses only Python's standard library."""
# Copyright (c) 2026 Huang Kai. SPDX-License-Identifier: LGPL-2.1-only
from pathlib import Path
import argparse
import hashlib
import json
import math
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    manifest = json.loads((HERE / 'input_manifest.json').read_text())
    for row in manifest:
        path = run / row['name']
        require(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
                'Input differs: ' + row['name'])
    text = (run / 'stdout.log').read_text()
    for pattern, label in [(r'Number of Nodes: 231\b', 'nodes'),
                           (r'Number of Dof: 1386\b', 'DOFs'),
                           (r'Number of Elements: 200\b', 'elements')]:
        require(re.search(pattern, text), 'Missing or incorrect ' + label)
    require(not re.search(r'(?<![a-zA-Z])(?:nan|[+-]?inf(?:inity)?)(?![a-zA-Z])', text, re.I),
            'Non-finite value in program output')
    indicators = [float(v) for v in re.findall(r'^(' + NUMBER + r')\s*$', text, re.M)]
    require(len(indicators) == 3 and all(0 <= v <= 1+1e-6 for v in indicators),
            'Expected finite MIR, MPC and MAC values')
    labels = ['DC_NelsonSingle', 'DC_NelsonAdd', 'DC_AdjointAdd', 'DC_AdjointNelson', 'DC_Symm']
    values = {}
    for label in labels:
        match = re.search(re.escape(label) + r':\s*(' + NUMBER + ')', text)
        require(match is not None, 'Missing method result: ' + label)
        values[label] = float(match.group(1))
    require(all(math.isfinite(v) for v in values.values()), 'Non-finite sensitivity')
    reference = values['DC_AdjointNelson']
    require(abs(reference) > 1e-10, 'Degenerate representative sensitivity')
    # These scalar comparisons use the precision printed by the existing driver.
    for label in labels[:3]:
        require(abs(values[label]-reference) <= 1e-5*abs(reference)+1e-10,
                'Direct methods disagree at the representative element: ' + label)
    match = re.search(r'Diff AdjointAdd vs AdjointNelson Norm2:\s*(' + NUMBER + ')', text)
    require(match is not None and float(match.group(1)) < 1e-8,
            'ADNE/ADBM full-vector difference is too large')
    difference = float(match.group(1))
    match = re.search(r'Representative element index:\s*(\d+)', text)
    require(match is not None and 0 <= int(match.group(1)) < 200, 'Invalid element index')
    element = int(match.group(1))
    for method in ['NelsonSingle', 'NelsonAdd', 'AdjointNelosnAdd', 'AdjointNelson', 'Adjointsymm']:
        require(re.search(r'^' + method + r'.*CPU time:', text, re.M), 'Method did not finish: ' + method)
    vtk = run / 'Shell_2010.vtk'
    require(vtk.is_file() and vtk.stat().st_size > 0, 'Missing geometry export')
    aspm_delta = abs(values['DC_Symm']-reference) / abs(reference)
    record = dict(status='PASS', nodes=231, elements=200, dofs=1386,
                  read_mode='single_element_two_coeff', element_tag=2,
                  indicators=dict(zip(['MIR', 'MPC', 'MAC'], indicators)),
                  representative_element_index_zero_based=element,
                  representative_sensitivities=values,
                  adne_adbm_full_vector_difference_l2=difference,
                  aspm_representative_relative_difference=aspm_delta,
                  aspm_settings=dict(update_tolerance=1e-5, maximum_iterations=5,
                                     source='existing src/main/main.cpp defaults'),
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                                 [ROOT/'src/main/main.cpp', ROOT/'src/Algorithm/TOP.cpp']},
                  scope='End-to-end execution of the published FEM driver; not a paper-data reproduction or an ASPM 1e-8 accuracy certificate.')
    (run / 'validation.json').write_text(json.dumps(record, indent=2) + '\n')
    print('PASS: actual PMTop run; 231 nodes, 200 elements, 1386 DOFs; all five methods completed.')
    print('ASPM representative relative difference:', format(aspm_delta, '.6g'))


if __name__ == '__main__':
    main()
