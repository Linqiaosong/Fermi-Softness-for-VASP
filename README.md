# Fermi softness calculation for Vienna Ab initio Simulation Package (VASP)

## Update

* 1.1.0: **Big update:** 
    1. Rewrote the code.
    2. Use Bader atomic division instead of W-S division. 
    3. No longer supported W-S division. 
    4. Changed the output format of FSCAR. 
    5. Support Fermi-Softness calculation of insulators and semiconductors. (experiment) 
* 1.0.1: Bug fixed. Change intergral to sum with weight.

## Introduction
In 2016, Huang and Zhuang proposed the theory and applications of a concept dubbed "Fermi softness", which distinguishes itself by enabling prediction of surface reactivity with spatial as well as atomic resolution. Herein, we provide a script depending on Vienna Ab initio Simulation Package (VASP) output for calculating Fermi softness.

## Dependency
* Anaconda3: [Installation](https://www.anaconda.com/products/individual#Downloads)
* ASE: [Installation](https://wiki.fysik.dtu.dk/ase/install.html)
* vaspkit: [Installation](https://vaspkit.com/installation.html)
* Bader: [Installation](http://theory.cm.utexas.edu/henkelman/code/bader/)

## Download script
[Click to get script](https://github.com/Linqiaosong/Fermi-Softness-for-VASP/releases/download/1.1.0/runfs.py)

## Download tutorials
[Click to get tutorials](https://github.com/Linqiaosong/Fermi-Softness-for-VASP/releases/download/1.1.0/How-to-calculate-Fermi-Softness.pdf)

## Reference
* B. Huang, L. Xiao, J. Lu, L. Zhuang, Angew. Chem. Int. Ed. 2016, 55, 6239–6243

## What is Fermi-Softness

![image](https://github.com/Linqiaosong/Fermi-Softness-for-VASP/blob/main/img/img.jpg)


![image](https://github.com/Linqiaosong/Fermi-Softness-for-VASP/blob/main/img/img2.png)
