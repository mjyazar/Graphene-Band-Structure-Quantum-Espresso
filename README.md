# Graphene Band Structure (QE + ASE)

This repository constructs Monolayer using Density Functional Theory (DFT). The grraphene undergoes relaxation and convergence testing, with a subsequent computation and plotting of the electronic band structure using Quantum Espresso (QE) and the Atomic Simulation Environment (ASE).

---

## Structure

```
.
|-- main.py    # main script for execution
|-- outputs/
|   |-- conv/  # Convergence test results (.txt)
|   |   |-- ECUTWFC Convergence.txt
|   |   |-- KGRIDS Convergence
|   |
|   |-- figs/  # Figures
|   |   |-- Graphene Band Structure.png
|   |   |-- KGRIDS Convergence.png
|   |   |-- ECUTWFC Convergence.png
|   |
|   |-- temp/      # Machine-level data from computations
|   |-- relax.pwi / relax.pwo  # Relaxation input/output
|   |-- bands.pwi / bands.pwo  # Bands input/output
|   |-- scf.pwi / scf.pwo.     # Self Consistent Field input/output
|
|-- pseudo/
|    |-- C.pbe-n-kjpaw_psl.1.0.0.UPF  # psuedopotential of Carbon (https://pseudopotentials.quantum-espresso.org/)


```

---


## Plots

### Graphene Band Structure
![Band Structure](<outputs/figs/Graphene Band Structure.png>)

### KGRIDS Convergence
![KGRIDS Convergence](<outputs/figs/KGRIDS Convergence.png>)
[KGRIDS Measurements](outputs/conv/KGRIDS%20Convergence.txt)

### ECUTWFC Convergence
![ECUTWFC Convergence](<outputs/figs/ECUTWFC Convergence.png>)
[ECUTWFC Measurements](outputs/conv/ECUTWFC%20Convergence.txt)
