from dataclasses import dataclass
from pathlib import Path
from ase import Atoms
import numpy as np


@dataclass
class Potential:
    data: np.ndarray
    atoms: object

@dataclass
class PotentialAveraged:
    coordinates: np.ndarray
    planar: np.ndarray
    macroscopic: np.ndarray


@dataclass
class ChargeDensity:
    data: np.ndarray
    atoms: object

@dataclass
class ChargeDensityAveraged:
    coordinates: np.ndarray
    planar: np.ndarray
    macroscopic: np.ndarray


@dataclass
class DOS:
    energy: np.ndarray
    dos: np.ndarray
    idos: np.ndarray


@dataclass
class LDOSAveraged:
    energies: np.ndarray
    coordinates: np.ndarray
    planar: np.ndarray


@dataclass
class System:
    name: str
    atoms: Atoms
    path: Path
    fermi_energy: float
    
    potential_raw: Potential
    potential_averaged: PotentialAveraged
    charge_density_raw: ChargeDensity
    charge_density_averaged: ChargeDensityAveraged
    dos: DOS
    ldos: LDOSAveraged


@dataclass
class Results:
    coupled: System
    bottom: System
    top: System
