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



#-------import------------
import subprocess
import numpy as np
import copy
from xml.etree.ElementTree import parse
#from pathlib import Path
from ase.io.cube import read_cube_data
from ase.atoms import Atoms
from ase.io import read
from ase.units import Bohr




#------uniform wavefunction----------
def uniform(atoms, data=None, origin=None):
# return data (np.array[i,j,k])

    dx=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])

    if data is None:
        data = np.ones((2, 2, 2))
    data = np.asarray(data)

    if data.dtype == complex:
        data = np.abs(data)

    if origin is None:
        origin = np.zeros(3)
    else:
        origin = np.asarray(origin) / Bohr

    for i in range(3):
        n = data.shape[i]
        d = atoms.cell[i] / n / Bohr
        dx[i] = d    

    s=np.linalg.det(dx)*np.sum(data)

    data=data/s

    return data



#-------write cube file function----------
def write_cube(fileobj, atoms, data=None, origin=None, comment=None):

    if data is None:
        data = np.ones((2, 2, 2))
    data = np.asarray(data)

    if data.dtype == complex:
        data = np.abs(data)

    if comment is None:
        comment = 'Cube file from ASE, written on ' + time.strftime('%c')
    else:
        comment = comment.strip()
    fileobj.write(comment)

    fileobj.write('\nOUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n')

    if origin is None:
        origin = np.zeros(3)
    else:
        origin = np.asarray(origin) / Bohr

    fileobj.write('{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n'
                  .format(len(atoms), *origin))

    for i in range(3):
        n = data.shape[i]
        d = atoms.cell[i] / n / Bohr
        fileobj.write('{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n'.format(n, *d))



    positions = atoms.positions / Bohr
    numbers = atoms.numbers
    for Z, (x, y, z) in zip(numbers, positions):
        fileobj.write('{0:5}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n'
                      .format(Z, 0.0, x, y, z))

    k=len(data[0,0,:])
    data=np.reshape(data,(-1,k))
    for i in range(0,len(data[:,0])):
        for j in range(0,k):
            fileobj.write(f"  {data[i,j]:e}")
            if (j+1) % 6 == 0:
                fileobj.write("\n")
        if k % 6 != 0:
            fileobj.write("\n")
    fileobj.write("\n")



#-------get VASP parameters--------------
def get_paraments(file_dir):
# return para (dict{})

    fileobj=open(file_dir)
    et = parse(fileobj)
    root = et.getroot()

    # get eigenvalues
    eigen=np.array([])
    for eigenvalues in root.findall('./calculation/eigenvalues/array/set/set/set/r'):
        eigen=np.append(eigen,float(eigenvalues.text.split()[0]))

    # get efermi
    ef=float(root.findall('./calculation/dos/i[@name="efermi"]')[0].text)

    # get weight
    kweight=np.array([])
    w=root.findall('./kpoints/varray[@name="weights"]')[0]
    for weight in w.findall('./v'):
        kweight=np.append(kweight,float(weight.text))

    # get kpoint
    kpoint=len(kweight)

    # get ion number
    ion=int(root.findall('./atominfo/atoms')[0].text)

    # get ispin
    ispin=int(root.findall('./parameters/separator[@name="electronic"]/separator[@name="electronic spin"]/i[@name="ISPIN"]')[0].text)

    # band number
    band=int(len(eigen)/kpoint/ispin)

    fileobj.close()

    paraments={'NIONS':ion,
               'NKPTS':kpoint,
               'NBANDS':band,
               'Ef':[ef,ef],
               'ISPIN':ispin,
               'EIGENVAL':eigen,
               'WEIGHT':kweight}
    return paraments


#-------generate wavefunction.cube-------
def run_vaspkit_wfn(para,k_index,band_index,vaspkit_dir):
# band_index from 1 to band_number
# k_index from 1 to kpoint_number

    ispin=para['ISPIN']

    if ispin==1:
        # p1=Path(f"WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube")
        # p2=p1
        tmp=subprocess.getstatusoutput(f"ls WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube")
    else:
        # p1=Path(f"WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_UP.vasp.cube")
        # p2=Path(f"WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube")
        tmp=subprocess.getstatusoutput(f"ls WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_UP.vasp.cube WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube")
    
    # if p1.exists() and p2.exists():
    #     return 0
    if "No such file or directory" not in tmp[1]:
        return 0
    else:
        vaspkit_ini=open('vaspkit.ini','w')
        vaspkit_ini.write(f'51\n516\n{k_index:d}\n{band_index:d}\n')
        vaspkit_ini.close()
        tmp=subprocess.getstatusoutput(vaspkit_dir+" < vaspkit.ini > vaspkit.log")
        if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
            print("\n\t**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
            exit()


#------get eigenvalue-----------
def get_eigenvalue(para,k_index,band_index,spin):
# band_index from 1 to band_number
# k_index from 1 to kpoint_number    
# spin: 1 or 2
# return eigenvalue_k_n (float)
    kpoint_number=para['NKPTS']
    band_number=para['NBANDS']
    eigen=para['EIGENVAL']
    return eigen[(spin-1)*kpoint_number*band_number+(k_index-1)*band_number+band_index-1]


#-------calculate LFS----------
def calc_lfs(para,kbT,dfdd_threshold,intermediate_file_options,vaspkit_dir):
# return fs (np.array[i,j,k]) , atoms (ase.Atoms)

    kpoint_number=para['NKPTS']
    band_number=para['NBANDS']
    ispin=para['ISPIN']
    ef=para['Ef']
    ion_number=para['NIONS']
    kweight=para['WEIGHT']
    eigen=para['EIGENVAL']

    if ispin==2:
        tagspin=['_UP','_DW']
    else:
        tagspin=['']
    
    fs=[[],[]]

    for s in range(ispin):
        i=0
        spin=s+1
        print(f'\n\tStart calculating intergral of spin={spin:d}:\n\tKpoint\tBand\tE-Ef/eV\t\t-dFDD\t\t\tweight')
        for k in range(kpoint_number):
            for b in range(band_number):
                k_index=k+1
                band_index=b+1
                e_ef=get_eigenvalue(para,k_index,band_index,spin)-ef[s]
                dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                if dfdd >= dfdd_threshold:
                    run_vaspkit_wfn(para,k_index,band_index,vaspkit_dir)
                    i=i+1
                else:
                    continue
                
                data, atoms = read_cube_data(f'WFN_SQUARED_B{band_index:04d}_K{k_index:04d}'+tagspin[s]+'.vasp.cube')
                data=uniform(atoms,data)*dfdd

                if intermediate_file_options==False:
                    # p=Path('.')
                    # wfn=list(p.glob(f'WFN_SQUARED_B{band_index:04d}_K{k_index:04d}*'))
                    # for q in wfn:
                    #     q.unlink(True)
                    tmp=subprocess.getstatusoutput(f"rm WFN_SQUARED_B{band_index:04d}_K{k_index:04d}*.vasp.cube")
                
                print(f"\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k]:.6f}")

                if i==1:
                    fs[s]=data*kweight[k]
                else:
                    fs[s]=fs[s]+data*kweight[k]            
    
    return fs,atoms


#-------write LFScube-----------
def write_lfs(para,fs,atoms,tag=''):

    ispin=para['ISPIN']

    if ispin==2:
        tagspin=['_UP','_DW']
    else:
        tagspin=['']
    for s in range(len(tagspin)):
        fs_file=open("LFS"+tagspin[s]+tag+".cube",'w')
        write_cube(fs_file,atoms,fs[s],[0.0,0.0,0.0],"Fermi_Softness"+tagspin[s]+tag)


#-------write FSCAR----------
def write_fscar(para,bader_dir,tag=''):

    ispin=para['ISPIN']

    if ispin==1:
        subprocess.getstatusoutput(bader_dir+' LFS'+tag+'.cube')
        # p=Path('ACF.dat')
        # target=Path('FSCAR'+tag)
        # p.rename(target)
        subprocess.getstatusoutput('mv ACF.dat FSCAR'+tag)
    else:
        subprocess.getstatusoutput(bader_dir+' LFS_UP'+tag+'.cube')
        # p=Path('ACF.dat')
        # target=Path('FSCAR_UP'+tag)
        # p.rename(target)
        subprocess.getstatusoutput('mv ACF.dat FSCAR_UP'+tag)
        subprocess.getstatusoutput(bader_dir+' LFS_DW'+tag+'.cube')
        # p=Path('ACF.dat')
        # target=Path('FSCAR_DW'+tag)
        # p.rename(target)
        subprocess.getstatusoutput('mv ACF.dat FSCAR_DW'+tag)


#-------FS modudle-----------
def run_fs(kbT,dfdd_threshold,band_gap,intermediate_file_options,bader_dir,vaspkit_dir):
    #----------Initialization---------
    print('''
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
    ''')

    para=get_paraments('vasprun.xml')

    kpoint_number=para['NKPTS']
    band_number=para['NBANDS']
    ispin=para['ISPIN']
    ef=para['Ef']
    ion_number=para['NIONS']
    kweight=para['WEIGHT']
    eigen=para['EIGENVAL']

    if ispin != 1 and ispin != 2:
        print('\n\t**** !!!! ISPIN error !!!! ****')
        exit()



    print(f'''
    Parameters:
    Electron temperature    =    {kbT:.6f} eV
    dFDD threshold          =    {dfdd_threshold:f}
    Fermi Energy            =    {ef[0]:.6f} eV
    ISPIN                   =    {ispin:d}
    Kpoint Numbers          =    {kpoint_number:d}
    Band Numbers            =    {band_number:d}
    Ion Numbers             =    {ion_number:d}
    CBM Energy              =    {band_gap['CBM']} eV
    VBM Energy              =    {band_gap['VBM']} eV
    Save Intermediate Files =    {intermediate_file_options}
    Bader PATH              =    {bader_dir:s}
    vaspkit PATH            =    {vaspkit_dir:s}

    Initialization is complete, start calculating:
    ''')

    print('\tKpoint\tWeight')

    for i in range(len(kweight)):
        print(f'\t{i+1}\t{kweight[i]}')


    if intermediate_file_options==True:
        # p=Path('WFNSQR')
        # wfn=list(p.glob('WFN_SQUARED_*'))
        # for q in wfn:
        #     target=q.name
        #     q.link_to(target)
        #     q.unlink(True)
            
        tmp=subprocess.getstatusoutput(f"ls ./WFNSQR/WFN_SQUARED_*.vasp.cube")
        if "No such file or directory" not in tmp[1]:
            subprocess.getstatusoutput(f"mv ./WFNSQR/WFN_SQUARED_*.vasp.cube .")
    #----------End:Initialization--------------


    if band_gap['CBM'] == [0.0] and band_gap['VBM'] == [0.0]:
        # no gap, calculate FS
        fs,atoms=calc_lfs(para,kbT,dfdd_threshold,intermediate_file_options,vaspkit_dir)
        write_lfs(para,fs,atoms)
        write_fscar(para,bader_dir)

    else:
        #----------calculate CB--------
        para_cbm=copy.deepcopy(para)
        # remove band under E_CBM
        for i in range(len(para_cbm['EIGENVAL'])):
            if para_cbm['EIGENVAL'][i] < min(band_gap['CBM']):
                para_cbm['EIGENVAL'][i] = 99.0
        # change Ef to E_CBM
        para_cbm['Ef']=band_gap['CBM']
        # calculate FS
        fs,atoms=calc_lfs(para_cbm,kbT,dfdd_threshold,intermediate_file_options,vaspkit_dir)
        write_lfs(para_cbm,fs,atoms,'_CB')
        write_fscar(para_cbm,bader_dir,'_CB')    

        #---------calculate VB-----------
        para_vbm=copy.deepcopy(para)
        # remove band above E_VBM
        for i in range(len(para_vbm['EIGENVAL'])):
            if para_vbm['EIGENVAL'][i] > max(band_gap['VBM']):
                para_vbm['EIGENVAL'][i] = -99.0
        # change Ef to E_CBM
        para_vbm['Ef']=band_gap['VBM']
        # calculate FS
        fs,atoms=calc_lfs(para_vbm,kbT,dfdd_threshold,intermediate_file_options,vaspkit_dir)
        write_lfs(para_vbm,fs,atoms,'_VB')
        write_fscar(para_vbm,bader_dir,'_VB')      

    #-----------save intermediate files-----------
    if intermediate_file_options==True:
        # p=Path('WFNSQR')
        # p.mkdir(exist_ok=True)
        # q=Path('.')
        # wfn=list(q.glob('WFN_SQUARED_*'))
        # for f in wfn:
        #     target= p / f
        #     f.link_to(target)
        #     f.unlink(True)
        subprocess.getstatusoutput('mkdir WFNSQR')
        subprocess.getstatusoutput('mv WFN_SQUARED* WFNSQR')


    #-----------remove temp----------
    # Path('vaspkit.ini').unlink(True)
    # Path('vaspkit.log').unlink(True)
    # Path('AVF.dat').unlink(True)
    # Path('BCF.dat').unlink(True)
    subprocess.getstatusoutput('rm vaspkit.ini vaspkit.log AVF.dat BCF.dat')

    #-----------print success--------
    print('\nThe calculation ends normally.')



#----------main-----------------
if __name__ == "__main__":
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
    run_fs(kbT,dfdd_threshold,band_gap,intermediate_file_options,bader_dir,vaspkit_dir)
