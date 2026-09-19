from dataclasses import dataclass
from pathlib import Path
from ase import Atoms
import numpy as np


@dataclass
class Potential:
    data: np.ndarray
    atoms: object
    averaged: np.ndarray
    z: np.ndarray


@dataclass
class PotentialAveraged:
    coordinates: np.ndarray
    planar: np.ndarray


@dataclass
class ChargeDensity:
    data: np.ndarray
    atoms: object
    averaged: np.ndarray
    z: np.ndarray

@dataclass
class ChargeDensityAveraged:
    coordinates: np.ndarray
    planar: np.ndarray


@dataclass
class LDOS:
    averaged_ldos: np.ndarray
    z: np.ndarray
    energies: np.ndarray
    atoms: Atoms

@dataclass
class LDOSAveraged:
    energies: np.ndarray
    coordinates: np.ndarray
    planar: np.ndarray


@dataclass
class DOS:
    energy: np.ndarray
    dos: np.ndarray
    idos: np.ndarray

@dataclass
class System:
    name: str
    atoms: Atoms
    path: Path
    fermi_energy: float
    
    potential_raw: Potential
    charge_density_raw: ChargeDensity
    dos: DOS
    ldos: LDOS


@dataclass
class Results:
    coupled: System
    bottom: System
    top: System
