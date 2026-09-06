# Graphene Band Structure (Quantum ESPRESSO + ASE)

This repository constructs monolayer and bilayer graphene using Density Functional Theory (DFT). The graphene structures undergo relaxation, computation of the self-consistent-field (SCF), convergence testing, with a subsequent computation and plotting of the electronic band structure using Quantum Espresso (QE) and the Atomic Simulation Environment (ASE).


---

## Structure

```
.
├── config.py
├── convergence.py
├── graphene.py
├── main.py
├── outputs
│   └── figures
│       ├── Bottom Layer DOS.png
│       ├── Coupled Layer DOS.png
│       └── Top Layer DOS.png
├── plotting.py
├── pseudo
│   └── C.pbe-n-kjpaw_psl.1.0.0.UPF
├── qe
│   ├── dos.py
│   ├── pp.py
│   ├── pw.py
│   └── runner.py
└── README.md
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
