"""Backend for DDSIM."""

import logging
import time
import uuid
import warnings
from typing import Dict, List, Union

from qiskit import QiskitError, QuantumCircuit
from qiskit.compiler import assemble
from qiskit.providers import BackendV1, Options
from qiskit.providers.models import BackendConfiguration, BackendStatus
from qiskit.qobj import PulseQobj, QasmQobj, QasmQobjExperiment, Qobj
from qiskit.result import Result

from mqt.ddsim import CircuitSimulator, __version__
from mqt.ddsim.job import DDSIMJob

logger = logging.getLogger(__name__)


class QasmSimulatorBackend(BackendV1):
    """Python interface to MQT DDSIM"""

    SHOW_STATE_VECTOR = False

    @classmethod
    def _default_options(cls) -> Options:
        return Options(
            shots=None,
            parameter_binds=None,
            simulator_seed=None,
            approximation_step_fidelity=1.0,
            approximation_steps=0,
            approximation_strategy="fidelity",
        )

    def __init__(self, configuration=None, provider=None):
        conf = {
            "backend_name": "qasm_simulator",
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
                "ccx",
                "mcx_gray",
                "mcx_recursive",
                "mcx_vchain",
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
                "mcrx",
                "ry",
                "cry",
                "mcry",
                "rz",
                "crz",
                "mcrz",
                "p",
                "cp",
                "cu1",
                "mcphase",
                "sx",
                "csx",
                "sxdg",
                "swap",
                "cswap",
                "iswap",
                "dcx",
                "ecr",
                "rxx",
                "ryy",
                "rzz",
                "rzx",
                "xx_minus_yy",
                "xx_plus_yy",
                "snapshot",
            ],
            "memory": False,
            "n_qubits": 64,
            "coupling_map": None,
            "conditional": False,
            "max_shots": 1000000000,
            "open_pulse": False,
            "gates": [],
        }
        super().__init__(configuration=configuration or BackendConfiguration.from_dict(conf), provider=provider)

    def run(self, quantum_circuits: Union[QuantumCircuit, List[QuantumCircuit]], **options) -> DDSIMJob:
        if isinstance(quantum_circuits, (QasmQobj, PulseQobj)):
            msg = "QasmQobj and PulseQobj are not supported."
            raise QiskitError(msg)

        if not isinstance(quantum_circuits, list):
            quantum_circuits = [quantum_circuits]

        out_options = {}
        for key in options:
            if not hasattr(self.options, key):
                warnings.warn("Option %s is not used by this backend" % key, UserWarning, stacklevel=2)
            else:
                out_options[key] = options[key]
        circuit_qobj = assemble(quantum_circuits, self, **out_options)

        job_id = str(uuid.uuid4())
        local_job = DDSIMJob(self, job_id, self._run_job, circuit_qobj, **options)
        local_job.submit()
        return local_job

    def _run_job(self, job_id, qobj_instance: Qobj, **options) -> Result:
        self._validate(qobj_instance)

        start = time.time()
        result_list = [self.run_experiment(qobj_exp, **options) for qobj_exp in qobj_instance.experiments]
        end = time.time()

        result = {
            "backend_name": self.configuration().backend_name,
            "backend_version": self.configuration().backend_version,
            "qobj_id": qobj_instance.qobj_id,
            "job_id": job_id,
            "results": result_list,
            "status": "COMPLETED",
            "success": True,
            "time_taken": (end - start),
            "header": qobj_instance.header.to_dict(),
        }
        return Result.from_dict(result)

    def run_experiment(self, qobj_experiment: QasmQobjExperiment, **options) -> Dict:
        start_time = time.time()
        approximation_step_fidelity = options.get("approximation_step_fidelity", 1.0)
        approximation_steps = options.get("approximation_steps", 1)
        approximation_strategy = options.get("approximation_strategy", "fidelity")
        seed = options.get("seed", -1)

        sim = CircuitSimulator(
            qobj_experiment,
            approximation_step_fidelity=approximation_step_fidelity,
            approximation_steps=approximation_steps,
            approximation_strategy=approximation_strategy,
            seed=seed,
        )
        counts = sim.simulate(options.get("shots", 1024))
        end_time = time.time()
        counts_hex = {hex(int(result, 2)): count for result, count in counts.items()}

        result = {
            "header": qobj_experiment.header.to_dict(),
            "name": qobj_experiment.header.name,
            "status": "DONE",
            "time_taken": end_time - start_time,
            "seed": options.get("seed", -1),
            "shots": options.get("shots", 1024),
            "data": {"counts": counts_hex},
            "success": True,
        }
        if self.SHOW_STATE_VECTOR:
            result["data"]["statevector"] = sim.get_vector()
        return result

    def _validate(self, _quantum_circuit):
        return

    def status(self) -> BackendStatus:
        """Return backend status.
        Returns:
            BackendStatus: the status of the backend.
        """
        return BackendStatus(
            backend_name=self.name(),
            backend_version=self.configuration().backend_version,
            operational=True,
            pending_jobs=0,
            status_msg="",
        )
