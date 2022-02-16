
####################################################
#
#        Fermi-Softness Calculation v1.2           
#
#             Author: Qiaosong Lin
#            Wuhan University, China
#
#  Notice:
#  1) You have already finished Non-SCF calculation
#  2) Make sure vaspkit and bader are in your $PATH
#  3) Make sure INCAR OUTCAR WAVECAR POSCAR vasprun.xml exist
#  4) Make sure ASE was installed correctly
#
#  Website:
#  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
#
#####################################################


#-------parameters----------
kbT=0.4                            # Electron temperature (eV): recommended 0.4 by B. Huang

dfdd_threshold=0.001               # Derivation of Fermi-Dirac distribution threshold: recommended 0.001 by B. Huang

intermediate_file_options=False     # Save intermediate files? False or True (default: False)

bader_dir='bader'                  # Path of bader, if bader is in your $PATH, you don't need to change it

vaspkit_dir='vaspkit'              # Path of vaspkit, if vaspkit is in your $PATH, you don't need to change it

band_gap={'VBM':[0.0],             # If band gap exists (You might need to confirm the occupation of VBM and CBM):
          'CBM':[0.0]}             # non-spin polarization: set as 'VBM':[E_VBM],'CBM':[E_CBM] (Do not minus E_fermi)
                                   # spin polarization: set as 'VBM':[E_VBM_UP,E_VBM_DW],'CBM':[E_CBM_UP,E_CBM_DW]
                                   # Otherwise: set as 'VBM':[0.0],'CBM':[0.0]
#----------------------------    

from FermiSoftness.vasp import run_fs
run_fs(kbT,dfdd_threshold,band_gap,intermediate_file_options,bader_dir,vaspkit_dir)
        