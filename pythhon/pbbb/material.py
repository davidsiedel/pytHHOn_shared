from typing import Dict

import mgis.behaviour as mgis_bv
from mgis.behaviour import Hypothesis
from mgis.behaviour import Behaviour
from mgis.behaviour import MaterialDataManager
from mgis.behaviour import IntegrationType
from mgis.behaviour import MaterialStateManagerStorageMode
from pythhon.pbbb.field import Field
from pythhon.parameters import *


class Material:
    behaviour: Behaviour
    mat_data: MaterialDataManager
    stabilization_parameter: float
    lagrange_parameter: float
    integration_type: IntegrationType
    storage_mode: MaterialStateManagerStorageMode
    temperature: float

    def __init__(
        self,
        nq: int,
        library_path: str,
        library_name: str,
        hypothesis: Hypothesis,
        stabilization_parameter: float,
        lagrange_parameter: float,
        field: Field,
        parameters: Dict[str, float] = None,
        integration_type: IntegrationType = IntegrationType.IntegrationWithConsistentTangentOperator,
        storage_mode: MaterialStateManagerStorageMode = MaterialStateManagerStorageMode.ExternalStorage,
        temperature: float = 293.15,
    ):
        """

        Args:
            nq:
            library_path:
            library_name:
            hypothesis:
            stabilization_parameter:
            field:
            parameters:
            integration_type:
            storage_mode:
            temperature:
        """
        self.nq = nq
        if field.strain_type in [StrainType.DISPLACEMENT_SYMMETRIC_GRADIENT]:
            if field.stress_type == StressType.CAUCHY:
                self.behaviour = mgis_bv.load(library_path, library_name, hypothesis)
        elif field.strain_type == StrainType.DISPLACEMENT_TRANSFORMATION_GRADIENT:
            if field.stress_type == StressType.PIOLA_KIRCHOFF_1:
                opt = mgis_bv.FiniteStrainBehaviourOptions()
                opt.stress_measure = mgis_bv.FiniteStrainBehaviourOptionsStressMeasure.PK1
                opt.tangent_operator = mgis_bv.FiniteStrainBehaviourOptionsTangentOperator.DPK1_DF
                self.behaviour = mgis_bv.load(opt, library_path, library_name, hypothesis)
        self.mat_data = mgis_bv.MaterialDataManager(self.behaviour, self.nq)
        self.stabilization_parameter = stabilization_parameter
        self.lagrange_parameter = lagrange_parameter
        if not parameters is None:
            for key, val in parameters.items():
                mgis_bv.setParameter(self.behaviour, key, val)
        self.integration_type = integration_type
        self.storage_mode = storage_mode
        self.temperature = temperature

    def set_temperature(self):
        """

        """
        T = self.temperature * np.ones(self.nq)
        mgis_bv.setExternalStateVariable(self.mat_data.s0, "Temperature", T, self.storage_mode)
        mgis_bv.setExternalStateVariable(self.mat_data.s1, "Temperature", T, self.storage_mode)
