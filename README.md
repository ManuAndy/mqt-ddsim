[![PyPI](https://img.shields.io/pypi/v/mqt.ddsim?logo=pypi&style=flat-square)](https://pypi.org/project/mqt.ddsim/)
![OS](https://img.shields.io/badge/os-linux%20%7C%20macos%20%7C%20windows-blue?style=flat-square)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/actions/workflow/status/cda-tum/mqt-ddsim/cmake.yml?branch=main&style=flat-square&logo=github&label=c%2B%2B)](https://github.com/cda-tum/mqt-ddsim/actions/workflows/cmake.yml)
[![Bindings](https://img.shields.io/github/actions/workflow/status/cda-tum/mqt-ddsim/bindings.yml?branch=main&style=flat-square&logo=github&label=python)](https://github.com/cda-tum/mqt-ddsim/actions/workflows/bindings.yml)
[![Documentation](https://img.shields.io/readthedocs/ddsim?logo=readthedocs&style=flat-square)](https://ddsim.readthedocs.io/en/latest/)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/mqt-ddsim?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/mqt-ddsim)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/mqt-ddsim/main/docs/source/_static/mqt_light.png" width="60%">
    <img src="https://raw.githubusercontent.com/cda-tum/mqt-ddsim/main/docs/source/_static/mqt_dark.png" width="60%">
  </picture>
</p>

# MQT DDSIM - A quantum circuit simulator based on decision diagrams written in C++

A tool for classical quantum circuit simulation by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/).
DDSIM is part of the Munich Quantum Toolkit (MQT; formerly known as JKQ and developed by the [Institute for Integrated Circuits](https://iic.jku.at/) at the [Johannes Kepler University Linz](https://jku.at)).
It builds upon [MQT Core](https://github.com/cda-tum/mqt-core).

**Detailed documentation on all available formats, options, and algorithms is available at [ReadTheDocs](https://ddsim.readthedocs.io/en/latest/).**

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by creating an [issue](https://github.com/cda-tum/mqt-ddsim/issues) on GitHub.

## Getting Started

DDSIM bundled with the provider and backends for Qiskit is available via [PyPI](https://pypi.org/project/mqt.ddsim/) for Linux, macOS, and Windows.

```console
(venv) $ pip install mqt.ddsim
```

The following code gives an example on the usage:

```python3
from qiskit import *
from mqt import ddsim

circ = QuantumCircuit(3)
circ.h(0)
circ.cx(0, 1)
circ.cx(0, 2)

print(circ.draw(fold=-1))

backend = ddsim.DDSIMProvider().get_backend('qasm_simulator')

job = execute(circ, backend, shots=10000)
counts = job.result().get_counts(circ)
print(counts)
```

## System Requirements and Building

The implementation is compatible with any C++17 compiler and a minimum CMake version of 3.14.
Please refer to the [documentation](https://ddsim.readthedocs.io/en/latest/) on how to build the project.

Building (and running) is continuously tested under Linux, macOS, and Windows using the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).

## References

DDSIM has been developed based on methods proposed in the following papers:

- A. Zulehner and R. Wille, "[Advanced Simulation of Quantum Computations](https://www.cda.cit.tum.de/files/eda/2018_tcad_advanced_simulation_quantum_computation.pdf)," Transactions on CAD of Integrated Circuits and Systems (TCAD), 2019
- S. Hillmich, I.L. Markov, and R. Wille, "[Just Like the Real Thing: Fast Weak Simulation of Quantum Computation](https://www.cda.cit.tum.de/files/eda/2020_dac_weak_simulation_quantum_computation.pdf)," in Design Automation Conference (DAC), 2020
- S. Hillmich, R. Kueng, I. L. Markov, and R. Wille, "[As Accurate as Needed, as Efficient as Possible: Approximations in DD-based Quantum Circuit Simulation](https://www.cda.cit.tum.de/files/eda/2021_date_approximations_dd_baed_quantum_circuit_simulation.pdf)," in Design, Automation and Test in Europe (DATE), 2021
- L. Burgholzer, H. Bauer, and R. Wille, "[Hybrid Schrödinger-Feynman Simulation of Quantum Circuits With Decision Diagrams](https://www.cda.cit.tum.de/files/eda/2021_qce_hybrid_schrodinger_feynman_simulation_with_decision_diagrams.pdf)," Conference on Quantum Computing and Engineering (QCE), 2021
- L. Burgholzer, A.Ploier, and R. Wille, "[Exploiting Arbitrary Paths for the Simulation of Quantum Circuits with Decision Diagrams](https://www.cda.cit.tum.de/files/eda/2022_date_exploiting_arbitrary_paths_simulation_quantum_circuits_decision_diagrams.pdf)," in Design, Automation and Test in Europe (DATE), 2022
- T. Grurl, J. Fuß, and R. Wille, "[Noise-aware Quantum Circuit Simulation With Decision Diagrams](https://www.cda.cit.tum.de/files/eda/2022_tcad_noise-aware_quantum_circuit_simulation_with_decision_diagrams.pdf)," Transactions on CAD of Integrated Circuits and Systems (TCAD), 2022.
