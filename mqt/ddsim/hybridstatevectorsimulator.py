"""Backend for DDSIM Hybrid Schrodinger-Feynman Simulator."""

import logging
from math import log2

from qiskit.providers.models import BackendConfiguration
from qiskit.utils.multiprocessing import local_hardware_info

from mqt.ddsim import __version__
from mqt.ddsim.hybridqasmsimulator import HybridQasmSimulatorBackend

logger = logging.getLogger(__name__)


class HybridStatevectorSimulatorBackend(HybridQasmSimulatorBackend):
    """Python interface to MQT DDSIM Hybrid Schrodinger-Feynman Simulator"""

    SHOW_STATE_VECTOR = True

    def __init__(self, configuration=None, provider=None):
        conf = {
            "backend_name": "hybrid_statevector_simulator",
            "backend_version": __version__,
            "url": "https://github.com/cda-tum/mqt-ddsim",
            "simulator": True,
            "local": True,
            "description": "MQT DDSIM C++ simulator",
            "basis_gates": [
                "gphase",
                "id",
                "u0",
                "u1",
                "u2",
                "u3",
                "cu3",
                "x",
                "cx",
                "y",
                "cy",
                "z",
                "cz",
                "h",
                "ch",
                "s",
                "sdg",
                "t",
                "tdg",
                "rx",
                "crx",
                "ry",
                "cry",
                "rz",
                "crz",
                "p",
                "cp",
                "cu1",
                "sx",
                "csx",
                "sxdg",
                # 'swap', 'cswap', 'iswap',
                "snapshot",
            ],
            "memory": False,
            "n_qubits": int(log2(local_hardware_info()["memory"] * (1024**3) / 16)),
            "coupling_map": None,
            "conditional": False,
            "max_shots": 1000000000,
            "open_pulse": False,
            "gates": [],
        }
        super().__init__(configuration=configuration or BackendConfiguration.from_dict(conf), provider=provider)

    def _validate(self, _quantum_circuit):
        return
