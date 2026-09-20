# Dielectric Embedding of Bilayer Graphene using Quantum ESPRESSO

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
│   ├── field_0
│   │   ├── bottom      # relax/scf/nscf, charge, potential, dos, ldos
│   │   ├── coupled     # same, plus relax.pwi/relax.pwo
│   │   └── top         # same as bottom
│   └── figures
│       ├── charge_density
│       │   └── Subtracted_Charge_Density_0au.png
│       ├── dos
│       │   └── DOS_Subtracted_Fermi-aligned_0au.png
│       ├── ldos
│       │   └── LDOS_Subtracted_0au.png
│       └── potential
│           ├── Potential_Bottom_Layer_0au.png
│           ├── Potential_Coupled_Layer_0au.png
│           ├── Potential_Top_Layer_0au.png
│           └── Subtracted_Potential_0au.png
├── plotter.py
├── pseudo
│   ├── C.pbe-n-kjpaw_psl.1.0.0.UPF
│   └── C.upf
├── qe
│   ├── average.py
│   ├── dos.py
│   ├── pp.py
│   ├── projwfc.py
│   ├── pw.py
│   └── runner.py
├── README.md
├── results.py
├── run.sh
└── timer.py
```


---
## Results (Bilayer Graphene)

### LDOS
#### Subtracted LDOS (E-field=0au)
![Subtracted DOS](<outputs/figures/ldos/LDOS_Subtracted_0au.png>)

#### Subtracted LDOS (E-field=0.05au)
![Subtracted DOS](<outputs/figures/ldos/LDOS_Subtracted_0.05au.png>)

#### Subtracted LDOS (E-field=0.1au)
![Subtracted DOS](<outputs/figures/ldos/LDOS_Subtracted_0.1au.png>)


### DOS
#### Subtracted DOS
![Subtracted DOS](<outputs/figures/dos/DOS_Subtracted_Fermi-aligned_0au.png>)


### Subtracted Charge Density
#### Subtracted Charge Density
![Subtracted Charge Density](<outputs/figures/charge_density/Subtracted_Charge_Density_0au.png>)


### Potential
#### Coupled Layers
![Coupled Layers](<outputs/figures/potential/Potential_Coupled_Layer_0au.png>)

#### Subtracted Potential
![Subtracted Potential](<outputs/figures/potential/Subtracted_Potential_0au.png>)
