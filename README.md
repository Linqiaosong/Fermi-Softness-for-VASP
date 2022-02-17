# Fermi softness calculation

**Warning: In current version, FSCAR is wrong, that the Bader division was according to the local Fermi softness. The strict method, the Bader division should be according to the total charge density. Therefore, the correct method should be to calculate an total charge density (core + valance) grid (assuming it is called chgsum.cube), and then use ```bader LFS.cube -ref chgsum.cube```.**


## Update
* 1.2.0: **Big update:** 
    1. Can be installed by pip.
    2. New input file.
    3. Support QE and CP2K. (experiment)
* 1.1.3: Bug fixed. Fixed the bug of ISPIN error.
* 1.1.2: Bug fixed. Fixed the bug of intermidiate files path.
* 1.1.1: Bug fixed. Fixed the bug of ISPIN error.
* 1.1.0: **Big update:** 
    1. Rewrote the code.
    2. Use Bader atomic division instead of W-S division. 
    3. No longer supported W-S division. 
    4. Changed the output format of FSCAR. 
    5. Support Fermi-Softness calculation of insulators and semiconductors. (experiment) 
* 1.0.1: Bug fixed. Change intergral to sum with weight.

## Introduction
In 2016, Huang and Zhuang proposed the theory and applications of a concept dubbed "Fermi softness", which distinguishes itself by enabling prediction of surface reactivity with spatial as well as atomic resolution. Herein, we provide a script for calculating Fermi softness.

## Dependency
* Anaconda3 (Python >= 3.8.5): [Installation](https://www.anaconda.com/products/individual#Downloads)
* ASE: [Installation](https://wiki.fysik.dtu.dk/ase/install.html)
* vaspkit (for VASP Users): [Installation](https://vaspkit.com/installation.html)
* Bader: [Installation](http://theory.cm.utexas.edu/henkelman/code/bader/)

## Installation

You can install ```FermiSoftness``` by ```pip```:
```bash
pip install FermiSoftness
```
or use ```git``` to download source code:
```bash
git clone https://github.com/Linqiaosong/Fermi-Softness-for-VASP.git
cd Fermi-Softness-for-VASP
pip install -e .
```

## Download tutorials
[Tutorial for VASP](https://github.com/Linqiaosong/fermi-softness/releases/download/1.2.0/tutorial-vasp.pdf)


[Tutorial for QE](https://github.com/Linqiaosong/fermi-softness/releases/download/1.2.0/tutorial-qe.pdf)


[Tutorial for CP2K](https://github.com/Linqiaosong/fermi-softness/releases/download/1.2.0/tutorial-cp2k.pdf)

## Reference
* [B. Huang, L. Xiao, J. Lu, L. Zhuang, Angew. Chem. Int. Ed. 2016, 55, 6239–6243](https://onlinelibrary.wiley.com/doi/abs/10.1002/ange.201601824)

## What is Fermi-Softness

![image](https://github.com/Linqiaosong/fermi-softness/blob/main/img/img.jpg)


![image](https://github.com/Linqiaosong/Fermi-Softness-for-VASP/blob/main/img/img2.png)
