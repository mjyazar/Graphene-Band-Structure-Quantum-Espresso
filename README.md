# Graphene Band Structure (Quantum ESPRESSO + ASE)

This repository constructs monolayer and bilayer graphene using Density Functional Theory (DFT). The graphene structures undergo relaxation, computation of the self-consistent-field (SCF), convergence testing, with a subsequent computation and plotting of the electronic band structure using Quantum Espresso (QE) and the Atomic Simulation Environment (ASE).


---

## Structure

```
.
|-- main.py         # main script for execution
|-- graphene.py     # Construction of graphene structure with GrapheneStructure class
|-- qe.py           # QE relaxation, scf, and band calculations
|-- convergence.py  # convergence testing
|-- pseudo/
|   |-- C.pbe-n-kjpaw_psl.1.0.0.UPF  # psuedopotential of Carbon (https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency)
|-- outputs/
    |-- monolayer/
    |   |-- relax.pwi / relax.pwo  # Relaxation input/output
    |   |-- scf.pwi / scf.pwo.     # Self Consistent Field input/output
    |   |-- bands.pwi / bands.pwo  # Bands input/output
    |   |-- data/                  # Machine-level QE data from computations
    |   |-- convergence/           # Convergence test results
    |       |-- kgrid/ / ecutwfc/  # per-parameter subrun calculation data
    |       |-- data/              # Machine-level QE data from computations
    |       |-- kgrid Convergence.txt
    |       |-- ecutwfc Convergence.txt
    |
    |-- bilayer/  # same structure as monolayer
    |-- figures/  # Band structure and convergence testing figures
        |-- Monolayer Graphene Band Structure.png
        |-- Bilayer Graphene Band Structure.png
        |-- Monolayer kgrid Convergence.png
        |-- Monolayer ecutwfc Convergence.png
        |-- Bilayer kgrid Convergence.png
        |-- Bilayer ecutwfc Convergence.png
```

---

## Results

### Monolayer Graphene
#### Band Structure
![Band Structure](<outputs/figures/Monolayer Graphene Band Structure.png>)

#### kgrid Convergence
![kgrid Convergence](<outputs/figures/Monolayer kgrid Convergence.png>)

[kgrid Measurements](outputs/monolayer/convergence/kgrid%20Convergence.txt)

### ecutwfc Convergence
![ecutwfc Convergence](<outputs/figures//Monolayer ecutwfc Convergence.png>)

[ecutwfc Measurements](outputs/monolayer/convergence/ecutwfc%20Convergence.txt)


### Bilayer Graphene

#### Band Structure
![Band Structure](<outputs/figures/Bilayer Graphene Band Structure.png>)

#### kgrid Convergence
![kgrid Convergence](<outputs/figures/Bilayer kgrid Convergence.png>)

[kgrid Measurements](outputs/bilayer/convergence/kgrid%20Convergence.txt)

### ecutwfc Convergence
![ecutwfc Convergence](<outputs/figures/Bilayer ecutwfc Convergence.png>)

[ecutwfc Measurements](outputs/bilayer/convergence/ecutwfc%20Convergence.txt)
