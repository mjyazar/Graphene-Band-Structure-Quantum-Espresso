from dataclasses import dataclass
from pathlib import Path
from ase import Atoms
import numpy as np


@dataclass
class PotentialRaw:
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
class LDOS:
    data: np.ndarray
    atoms: object


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
    
    potential_raw: PotentialRaw
    potential_averaged: PotentialAveraged
    dos: DOS


@dataclass
class Results:
    coupled: System
    bottom: System
    top: System
