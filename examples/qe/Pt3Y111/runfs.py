
####################################################
#
#        Fermi-Softness Calculation v1.2           
#
#             Author: Qiaosong Lin
#            Wuhan University, China
#
#  Notice:
#  1) You have already finished Non-SCF calculation
#  2) Make sure pp.x and bader are in your $PATH
#  3) Make sure prefix and outdir were set correctly
#  4) Make sure ASE was installed correctly
#
#  Website:
#  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
#
#####################################################


#-------parameters----------

prefix='pwscf'

outdir='./tmp'

kbT=0.4                            # Electron temperature (eV): recommended 0.4 by B. Huang

dfdd_threshold=0.001               # Derivation of Fermi-Dirac distribution threshold: recommended 0.001 by B. Huang

intermediate_file_options=False     # Save intermediate files? False or True (default: False)

bader_dir='bader'                  # Path of bader, if bader is in your $PATH, you don't need to change it

pp_laucher='mpirun -np 4 pp.x'                 # Laucher of pp.x, e.g.: 'pp.x' or 'mpirun -np 4 pp.x'

band_gap={'VBM':[0.0],             # If band gap exists (You might need to confirm the occupation of VBM and CBM):
          'CBM':[0.0]}             # non-spin polarization: set as 'VBM':[E_VBM],'CBM':[E_CBM] (Do not minus E_fermi)
                                   # spin polarization: set as 'VBM':[E_VBM_UP,E_VBM_DW],'CBM':[E_CBM_UP,E_CBM_DW]
                                   # Otherwise: set as 'VBM':[0.0],'CBM':[0.0]
#----------------------------              

from FermiSoftness.qe import run_fs
run_fs(prefix,outdir,kbT,dfdd_threshold,band_gap,intermediate_file_options,bader_dir,pp_laucher)
        