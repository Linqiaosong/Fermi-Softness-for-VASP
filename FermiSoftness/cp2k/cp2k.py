####################################################
#
#        Fermi-Softness Calculation v1.2          
#
#             Author: Qiaosong Lin
#            Wuhan University, China
#
#  Notice:
#  1) You have already finished single point calculation with &MO_CUBES
#  2) Make sure bader is in your $PATH
#  3) Make sure ASE was installed correctly
#
#  Website:
#  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
#
#####################################################

#-------import------------
import subprocess
import numpy as np
import copy
import re
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

    data=np.power(data,2)
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



def get_paraments(filename,ispin):
    result_homo_spin=[[],[]]
    result_lumo_spin=[[],[]]
    eigen=[[],[]]
    with open(filename,'r') as file_obj:
        str=file_obj.read()
    pattern=re.compile('Fermi Energy .*')
    result_ef=re.search(pattern,str).group()
    ef=float(result_ef.split()[-1])
    pattern=re.compile('Eigenvalues of the occupied subspace spin.+Fermi Energy',re.S)
    result_homo=re.search(pattern,str).group()
    result_homo=result_homo.splitlines()
    pattern=re.compile('Lowest eigenvalues of the unoccupied subspace spin.+ENERGY\|',re.S)
    result_lumo=re.search(pattern,str).group()
    result_lumo=result_lumo.splitlines()
    if ispin==1:
        result_homo_spin[0]='\n'.join(result_homo[2:len(result_homo)//ispin-1]).split()
        result_lumo_spin[0]='\n'.join(result_lumo[2:len(result_lumo)//ispin-1]).split()
        eigen[0]=result_homo_spin[0]+result_lumo_spin[0]
        eigen[0]=np.array(list(map(float,eigen[0])))*27.2114
    elif ispin==2:
        result_homo_spin[0]='\n'.join(result_homo[2:len(result_homo)//ispin-1]).split()
        result_homo_spin[1]='\n'.join(result_homo[len(result_homo)//ispin+3:len(result_homo)-1]).split()
        result_lumo_spin[0]='\n'.join(result_lumo[2:len(result_lumo)//ispin-1]).split()
        result_lumo_spin[1]='\n'.join(result_lumo[len(result_lumo)//ispin+2:len(result_lumo)-2]).split()
        eigen[0]=result_homo_spin[0]+result_lumo_spin[0]
        eigen[1]=result_homo_spin[1]+result_lumo_spin[1]
        eigen[0]=np.array(list(map(float,eigen[0])))*27.2114
        eigen[1]=np.array(list(map(float,eigen[1])))*27.2114
    else:
        print('ispin is error!')
        exit(1)

    paraments={'NBANDS':len(eigen[0]),
               'Ef':[ef,ef],
               'ISPIN':ispin,
               'EIGENVAL':eigen}
    return paraments


#-------calculate LFS----------
def calc_lfs(para,kbT,dfdd_threshold,project_name):
# return fs (np.array[i,j,k]) , atoms (ase.Atoms)

    band_number=para['NBANDS']
    ispin=para['ISPIN']
    ef=para['Ef']
    eigen=para['EIGENVAL']
    
    fs=[[],[]]

    for s in range(ispin):
        i=0
        spin=s+1
        print(f'\n\tStart calculating intergral of spin={spin:d}:\n\tBand\tE-Ef/eV\t\t-dFDD')
        for b in range(band_number):
            e_ef=eigen[s][b]-ef[s]
            dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
            if dfdd >= dfdd_threshold:
                data, atoms = read_cube_data(project_name+f'-WFN_{b+1:05d}_{s+1:d}-1_0.cube')
                data=uniform(atoms,data)*dfdd
                i=i+1
            else:
                continue
            
            print(f"\t{b+1:d}\t{e_ef:.6f}\t{dfdd:.8e}")

            if i==1:
                fs[s]=data
            else:
                fs[s]=fs[s]+data         
    
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
def run_fs(kbT,dfdd_threshold,band_gap,bader_dir,filename,project_name,ispin):
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
    #  1) You have already finished single point calculation with &MO_CUBES
    #  2) Make sure bader is in your $PATH
    #  3) Make sure ASE was installed correctly
    #
    #  Website:
    #  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
    #
    #####################################################
    ''')

    para=get_paraments(filename,ispin)

    band_number=para['NBANDS']
    ispin=para['ISPIN']
    ef=para['Ef']
    eigen=para['EIGENVAL']



    print(f'''
    Parameters:
    Electron temperature    =    {kbT:.6f} eV
    dFDD threshold          =    {dfdd_threshold:f}
    Fermi Energy            =    {ef[0]:.6f} eV
    ISPIN                   =    {ispin:d}
    Band Numbers            =    {band_number:d}
    CBM Energy              =    {band_gap['CBM']} eV
    VBM Energy              =    {band_gap['VBM']} eV
    Bader PATH              =    {bader_dir:s}

    Initialization is complete, start calculating:
    ''')


            


    if band_gap['CBM'] == [0.0] and band_gap['VBM'] == [0.0]:
        # no gap, calculate FS
        fs,atoms=calc_lfs(para,kbT,dfdd_threshold,project_name)
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
        fs,atoms=calc_lfs(para,kbT,dfdd_threshold,project_name)
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
        fs,atoms=calc_lfs(para,kbT,dfdd_threshold,project_name)
        write_lfs(para_vbm,fs,atoms,'_VB')
        write_fscar(para_vbm,bader_dir,'_VB')      



    #-----------remove temp----------
    # Path('vaspkit.ini').unlink(True)
    # Path('vaspkit.log').unlink(True)
    # Path('AVF.dat').unlink(True)
    # Path('BCF.dat').unlink(True)
    subprocess.getstatusoutput('rm AVF.dat BCF.dat')

    #-----------print success--------
    print('\nThe calculation ends normally.')



#----------main-----------------
if __name__ == "__main__":
    #-------parameters----------
    filename='pt3y.out'
    project_name='pt3y'
    ispin=1
    kbT=0.4                            # Electron temperature (eV): recommended 0.4 by B. Huang
    dfdd_threshold=0.001               # Derivation of Fermi-Dirac distribution threshold: recommended 0.001 by B. Huang
    bader_dir='bader'                  # Path of bader, if bader is in your $PATH, you don't need to change it
    band_gap={'VBM':[0.0],             # If band gap exists (You might need to confirm the occupation of VBM and CBM):
            'CBM':[0.0]}             # non-spin polarization: set as 'VBM':[E_VBM],'CBM':[E_CBM] (Do not minus E_fermi)
                                    # spin polarization: set as 'VBM':[E_VBM_UP,E_VBM_DW],'CBM':[E_CBM_UP,E_CBM_DW]
                                    # Otherwise: set as 'VBM':[0.0],'CBM':[0.0]
    #----------------------------  
    run_fs(kbT,dfdd_threshold,band_gap,bader_dir,filename,project_name,ispin)

