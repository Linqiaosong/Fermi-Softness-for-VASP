####################################################
#
#        Fermi-Softness Calculation v1.0           
#
#             Author: Qiaosong Lin
#            Wuhan University, China
#
#  Notice:
#  1) You have already finished Non-SCF calculation
#  2) Make sure vaspkit is in your $PATH
#  3) Make sure INCAR OUTCAR WAVECAR POSCAR EIGENVAL 
#     PROCAR DOSCAR exist
#  4) Make sure ASE was installed correctly
#
#  Website:
#  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
#
#####################################################


#-------parameters----------
kbT=0.4                            # Electron temperature (eV): recommended 0.4 by B. Huang
dfdd_threshold=0.001               # Derivation of Fermi-Dirac distribution threshold: recommended 0.001 by B. Huang
intermediate_file_options=True     # Save intermediate files? False or True (default: False) 
#-------End:parameters------




#-------import------------
import subprocess
import numpy as np
from scipy.integrate import simps
from ase.io.cube import read_cube_data
from ase.atoms import Atoms
from ase.io import read
from ase.units import Bohr
#--------End:import-----------



#------uniform wavefunction----------
def uniform(atoms, data=None, origin=None):
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
#-------End uniform wavefunction-----


#-------write cube file function----------
def write_cube(fileobj, atoms, data=None, origin=None, comment=None):

    dx=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])

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
#----------End:write cube file function-------






#----------main-----------------
if __name__ == "__main__":





    #----------Initialization---------
    print('''
    ####################################################
    #
    #        Fermi-Softness Calculation v1.0           
    #
    #             Author: Qiaosong Lin
    #            Wuhan University, China
    #
    #  Notice:
    #  1) You have already finished Non-SCF calculation
    #  2) Make sure vaspkit is in your $PATH
    #  3) Make sure INCAR OUTCAR WAVECAR POSCAR EIGENVAL 
    #     PROCAR DOSCAR exist
    #  4) Make sure ASE was installed correctly
    #
    #  Website:
    #  https://github.com/Linqiaosong/Fermi-Softness-for-VASP
    #
    #####################################################
    ''')
    
    outcar_line=subprocess.getstatusoutput('grep NKPTS OUTCAR')
    outcar=outcar_line[1].split()

    kpoint_number=int(outcar[3])
    band_number=int(outcar[-1])

    outcar_line=subprocess.getstatusoutput('grep ISPIN OUTCAR')
    outcar=outcar_line[1].split()
    ispin=int(outcar[2])

    outcar_line=subprocess.getstatusoutput('grep E-fermi OUTCAR')
    outcar=outcar_line[1].split()
    ef=float(outcar[2])

    outcar_line=subprocess.getstatusoutput('grep NIONS OUTCAR')
    outcar=outcar_line[1].split()
    ion_number=int(outcar[-1])

    kweight=[1.0 for i in range(kpoint_number)]
    outcar_line=subprocess.getstatusoutput(f'grep -A{kpoint_number:d} weights: OUTCAR')
    outcar=outcar_line[1].split()
    for i in range(kpoint_number):
        kweight[i]=float(outcar[11+4*(i+1)])

    print(f'''
    Parameters:
    Electron temperature    =    {kbT:.6f} eV
    dFDD threshold          =    {dfdd_threshold:f}
    Fermi Energy            =    {ef:.6f} eV
    ISPIN                   =    {ispin:d}
    Kpoint Numbers          =    {kpoint_number:d}
    Band Numbers            =    {band_number:d}
    Ion Numbers             =    {ion_number:d}
    Save Intermediate Files =    {intermediate_file_options:d}

    Initialization is complete, start calculating:
    ''')

    if intermediate_file_options==True:
        tmp=subprocess.getstatusoutput(f"ls ./WFNSQR/WFN_SQUARED_*.vasp.cube")
        if "No such file or directory" not in tmp[1]:
            subprocess.getstatusoutput(f"mv ./WFNSQR/WFN_SQUARED_*.vasp.cube .")
        tmp=subprocess.getstatusoutput(f"ls ./DOS/PDOS_A*.dat")
        if "No such file or directory" not in tmp[1]:
            subprocess.getstatusoutput(f"mv ./DOS/PDOS_A*.dat .")
    #----------End:Initialization--------------



    #---------- Generate and read Total DOS---------
    print("\nStart calculating total and projected Fermi-Softness.")
    print("\nGenerating TDOS file:")
    print("Running vaspkit for generating TDOS file.")
    vaspkit_ini=open('vaspkit.ini','w')
    vaspkit_ini.write('11\n111\n')
    vaspkit_ini.close()
    tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
    if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
        print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
        exit()
    #else:
        #print("---> Successed.")

    dos=np.loadtxt('TDOS.dat',dtype=np.float64)
    col_number=len(dos[0,:])
    row_number=len(dos[:,0])
    intelgral=[0 for i in range(col_number)]
    for j in range(col_number):
        if j==0:
            continue
        else:
            for i in range(row_number):
                e_ef=dos[i,0]
                dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                if dfdd < dfdd_threshold:
                    dos[i,j]=0
                else:
                    dos[i,j]=dos[i,j]*dfdd



    #-----------non spin polarization---------
    if ispin==1:

        fscar=open('FSCAR','w')



        #-----------Total FS-----------
        tmp=subprocess.getstatusoutput('head -1 TDOS.dat')
        tmp=tmp[1].split()
        tmp[0]=f'{0.0:f}'
        tmp='\t\t'.join(tmp)
        fscar.write('----Total Fermi-Softness.----\n')
        fscar.write(f'{tmp:s}\n')
        for j in [0,1]:
            if j==0:
                pass
            else:
                intelgral[j]=simps(dos[:,j],dos[:,0])
            fscar.write(f'{intelgral[j]:f}\t')
        fscar.write('\n\n')
        #---------End:Total FS---------------



        #-----------Projected FS of element---------

        # Generate Projected DOS of element
        print("\nGenerating PDOS file of each element:")
        print("Running vaspkit for generating PDOS file for each element.")
        vaspkit_ini=open('vaspkit.ini','w')
        vaspkit_ini.write('11\n113\nall\n')
        vaspkit_ini.close()
        tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
        if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
            print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
            exit()
        #else:
            #print("---> Successed.")

        datfile=subprocess.getstatusoutput('ls PDOS*.dat | grep "PDOS_[A-Z][(a-z|.)]"')
        datfile=datfile[1].split()

        tmp=subprocess.getstatusoutput(f'head -1 {datfile[0]:s}')
        tmp=tmp[1].split()
        tmp[0]='Element'
        tmp='\t\t'.join(tmp)
        fscar.write('----Projected Fermi-Softness of each element.----\n')
        fscar.write(f'{tmp:s}\n')


        # read PDOS and intergral
        for f in datfile:
            dos=np.loadtxt(f,dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            # write datas
            for j in range(col_number):
                if j==0:
                    fscar.write(f[5:-4]+'\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar.write(f'{intelgral[j]:f}\t')
            fscar.write('\n')
            
        fscar.write('\n')
        #--------End:Projected FS of element----------------



        #-----------Projected FS of atom----------
        # Generate Projected DOS of atom
        print("\nGenerating PDOS file of each atom:")
        for i in range(ion_number):
            tmp=subprocess.getstatusoutput(f"ls PDOS_A{i+1:d}.dat")
            if "No such file or directory" in tmp[1]:
                print(f"Running vaspkit for generating PDOS file of atom {i+1:d}.")
                vaspkit_ini=open('vaspkit.ini','w')
                vaspkit_ini.write(f'11\n112\n{i+1:d}\nall\n')
                vaspkit_ini.close()
                tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
                if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
                    print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
                    exit()
                #else:
                    #print("---> Successed.")
            #else:
                #print(f"---> PDOS file of Atom {i+1:d} exists, abord running vaspkit.")

        tmp=subprocess.getstatusoutput(f'head -1 PDOS_A1.dat')
        tmp=tmp[1].split()
        tmp[0]='Atom'
        tmp='\t\t'.join(tmp)
        fscar.write('----Projected Fermi-Softness of each atom, the atom list is the same as POSCAR.----\n')
        fscar.write(f'{tmp:s}\n')


        # read PDOS and intergral
        for f in range(ion_number):
            dos=np.loadtxt(f'PDOS_A{f+1:d}.dat',dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            # write datas
            for j in range(col_number):
                if j==0:
                    fscar.write(f'{f+1:d}\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar.write(f'{intelgral[j]:f}\t')
            fscar.write('\n')
            
        fscar.write('\n')
        #----------End:Projected FS of atom---------------


        fscar.close()

        #--------remove intermediate files------------
        if intermediate_file_options==False:
            subprocess.getstatusoutput('rm *DOS*.*')
            subprocess.getstatusoutput('rm SELECTED_ORBITALS_LIST')


        #-------------Local FS------------
        print("\nStart calculating local Fermi-Softness.")
        eigen_file=open('EIGENVAL',"r")
        eigen_text=eigen_file.readlines()
        eigen_file.close()

        eigen_up={}

        # read energy
        for i in range(6,(band_number+2)*kpoint_number+7):
            k_index=int((i-6)/(band_number+2))
            band_index=(i-6)%(band_number+2)-2
            if band_index < 0:
                continue
            else:
                print(f"\n[{k_index*band_number+(band_index+1):d}/{band_number*kpoint_number:d}] Reading kpoint={k_index+1:d} band={band_index+1:d} eigenvalue from EIGENVAL.")
                eigen_text[i]=eigen_text[i].split( )

                eigen_val_up=float(eigen_text[i][1])
                e_ef_up=eigen_val_up-ef
                dfdd_up=(1.0/kbT)*np.exp(e_ef_up/kbT)/(np.exp(e_ef_up/kbT)+1)/(np.exp(e_ef_up/kbT)+1)
                if dfdd_up < dfdd_threshold:
                    print(f"kpoint={k_index+1:d} band={band_index+1:d}: -dFDD={dfdd_up:.8e} is too low, ignore.")
                else:
                    eigen_up[eigen_val_up]=[k_index+1,band_index+1]
                    
        # sort energy
        eigen_val_up=sorted(eigen_up.keys())

        # generate and read wavefunction, intergral
        print('''\nStart calculating intergral:\nIndex\tKpoint\tBand\tE-Ef/eV\t\t-dFDD\t\t\tweight''')
        for i in range(len(eigen_val_up)):
            band_index=eigen_up[eigen_val_up[i]][1]
            k_index=eigen_up[eigen_val_up[i]][0]

            # generate wavefunction
            tmp=subprocess.getstatusoutput(f"ls WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube")
            if "No such file or directory" in tmp[1]:
                #print("Running vaspkit for generating wavefunction file.")
                vaspkit_ini=open('vaspkit.ini','w')
                vaspkit_ini.write(f'51\n516\n{k_index:d}\n{band_index:d}\n')
                vaspkit_ini.close()
                tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
                if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
                    print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
                    exit()
                #else:
                    #print("---> Successed.")
            #else:
                #print("---> Wavefunction file exists, abord running vaspkit.")
            #print(f"Reading WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube")

            # read wavefunction
            data, atoms = read_cube_data(f'WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube')
            data=uniform(atoms,data)


            # remove intermediate files
            if intermediate_file_options==False:
                tmp=subprocess.getstatusoutput(f"rm WFN_SQUARED_B{band_index:04d}_K{k_index:04d}.vasp.cube")


            # intergral
            e_ef=eigen_val_up[i]-ef
            dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
            data=data*dfdd
            if i==0:
                print(f"{i+1:d}/{len(eigen_val_up):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:f}")
                fs_up=data*kweight[k_index-1]
            else:
                #dE=eigen_val_up[i]-eigen_val_up[i-1]
                print(f"{i+1:d}/{len(eigen_val_up):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:f}")
                fs_up=fs_up+data*kweight[k_index-1]


        # write cube
        fs_up_file=open("LFS.cube",'w')
        write_cube(fs_up_file,atoms,fs_up,[0.0,0.0,0.0],"Fermi_Softness")
        fs_up_file.close()
        #----------End:Local FS--------------
    #--------End:non spin polarization----------------




    #-----------spin polarization--------------
    elif ispin==2:

        fscar=open('FSCAR_UP','w')
        fscar_dw=open('FSCAR_DW','w')
        



        #------------Total FS-------------
        # spin=1
        tmp=subprocess.getstatusoutput('head -1 TDOS.dat')
        tmp=tmp[1].split()
        tmp[0]='Spin'
        tmp[2]=' '
        tmp='\t\t'.join(tmp)
        fscar.write('----Total Fermi-Softness.----\n')
        fscar.write(f'{tmp:s}\n')
        for j in [0,1]:
            if j==0:
                intelgral[j]=1.0
            else:
                intelgral[j]=simps(dos[:,j],dos[:,0])
            fscar.write(f'{intelgral[j]:f}\t')
        fscar.write('\n\n')

        # spin=2
        tmp=subprocess.getstatusoutput('head -1 TDOS.dat')
        tmp=tmp[1].split()
        tmp[0]='Spin'
        tmp[1]=' '
        tmp='\t'.join(tmp)
        fscar_dw.write('----Total Fermi-Softness.----\n')
        fscar_dw.write(f'{tmp:s}\n')
        for j in [0,2]:
            if j==0:
                intelgral[j]=-2.0
            else:
                intelgral[j]=simps(dos[:,j],dos[:,0])
            fscar_dw.write(f'{-intelgral[j]:f}\t')
        fscar_dw.write('\n\n')
        #----------End:Total FS--------------





        #----------Projected FS of element---------
        print("\nGenerating PDOS file of each element:")
        print("Running vaspkit for generating PDOS file for each element.")
        vaspkit_ini=open('vaspkit.ini','w')
        vaspkit_ini.write('11\n113\nall\n')
        vaspkit_ini.close()
        tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
        if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
            print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
            exit()
        #else:
        #    print("---> Successed.")


        # spin=1
        datfile=subprocess.getstatusoutput('ls PDOS*.dat | grep "PDOS_[A-Z][(a-z|_)]*UP"')
        datfile=datfile[1].split()

        tmp=subprocess.getstatusoutput(f'head -1 {datfile[0]:s}')
        tmp=tmp[1].split()
        tmp[0]='Element'
        tmp='\t\t'.join(tmp)
        fscar.write('----Projected Fermi-Softness of each element.----\n')
        fscar.write(f'{tmp:s}\n')

        for f in datfile:
            dos=np.loadtxt(f,dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            for j in range(col_number):
                if j==0:
                    fscar.write(f[5:-7]+'\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar.write(f'{intelgral[j]:f}\t')
            fscar.write('\n')
            
        fscar.write('\n')


        # spin=2
        datfile=subprocess.getstatusoutput('ls PDOS*.dat | grep "PDOS_[A-Z][(a-z|_)]*DW"')
        datfile=datfile[1].split()

        tmp=subprocess.getstatusoutput(f'head -1 {datfile[0]:s}')
        tmp=tmp[1].split()
        tmp[0]='Element'
        tmp='\t\t'.join(tmp)
        fscar_dw.write('----Projected Fermi-Softness of each element.----\n')
        fscar_dw.write(f'{tmp:s}\n')

        for f in datfile:
            dos=np.loadtxt(f,dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            for j in range(col_number):
                if j==0:
                    fscar_dw.write(f[5:-7]+'\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar_dw.write(f'{-intelgral[j]:f}\t')
            fscar_dw.write('\n')
            
        fscar_dw.write('\n')
        #----------End:Projected FS of element--------------





        #---------Projected FS of atom------------
        print("\nGenerating PDOS file of each atom:")
        for i in range(ion_number):
            tmp=subprocess.getstatusoutput(f"ls PDOS_A{i+1:d}_UP.dat PDOS_A{i+1:d}_DW.dat")
            if "No such file or directory" in tmp[1]:
                print(f"Running vaspkit for generating PDOS file of atom {i+1:d}.")
                vaspkit_ini=open('vaspkit.ini','w')
                vaspkit_ini.write(f'11\n112\n{i+1:d}\nall\n')
                vaspkit_ini.close()
                tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
                if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
                    print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
                    exit()
                #else:
                #    print("---> Successed.")
            #else:
            #    print(f"---> PDOS file of Atom {i+1:d} exists, abord running vaspkit.")



        # spin=1
        tmp=subprocess.getstatusoutput(f'head -1 PDOS_A1_UP.dat')
        tmp=tmp[1].split()
        tmp[0]='Atom'
        tmp='\t\t'.join(tmp)
        fscar.write('----Projected Fermi-Softness of each atom, the atom list is the same as POSCAR.----\n')
        fscar.write(f'{tmp:s}\n')

        for f in range(ion_number):
            dos=np.loadtxt(f'PDOS_A{f+1:d}_UP.dat',dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            for j in range(col_number):
                if j==0:
                    fscar.write(f'{f+1:d}\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar.write(f'{intelgral[j]:f}\t')
            fscar.write('\n')
            
        fscar.write('\n')


        # spin=2
        tmp=subprocess.getstatusoutput(f'head -1 PDOS_A1_DW.dat')
        tmp=tmp[1].split()
        tmp[0]='Atom'
        tmp='\t\t'.join(tmp)
        fscar_dw.write('----Projected Fermi-Softness of each atom, the atom list is the same as POSCAR.----\n')
        fscar_dw.write(f'{tmp:s}\n')

        for f in range(ion_number):
            dos=np.loadtxt(f'PDOS_A{f+1:d}_DW.dat',dtype=np.float64)
            col_number=len(dos[0,:])
            row_number=len(dos[:,0])
            intelgral=[0 for i in range(col_number)]
            for j in range(col_number):
                if j==0:
                    continue
                else:
                    for i in range(row_number):
                        e_ef=dos[i,0]
                        dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
                        if dfdd < dfdd_threshold:
                            dos[i,j]=0
                        else:
                            dos[i,j]=dos[i,j]*dfdd

            for j in range(col_number):
                if j==0:
                    fscar_dw.write(f'{f+1:d}\t')
                else:
                    intelgral[j]=simps(dos[:,j],dos[:,0])
                    fscar_dw.write(f'{-intelgral[j]:f}\t')
            fscar_dw.write('\n')
            
        fscar_dw.write('\n')
        #---------End:Projected FS of atom---------------



        fscar.close()
        fscar_dw.close()



        #---------------remove intermediate files---------
        if intermediate_file_options==False:
            subprocess.getstatusoutput('rm *DOS*.*')
            subprocess.getstatusoutput('rm SELECTED_ORBITALS_LIST')




        #----------Local FS--------------
        print("\nStart calculating local Fermi-Softness.")
        eigen_file=open('EIGENVAL',"r")
        eigen_text=eigen_file.readlines()
        eigen_file.close()

        eigen_up={}
        eigen_dw={}

        for i in range(6,(band_number+2)*kpoint_number+7):
            k_index=int((i-6)/(band_number+2))
            band_index=(i-6)%(band_number+2)-2
            if band_index < 0:
                continue
            else:
                print(f"\n[{k_index*band_number+(band_index+1):d}/{band_number*kpoint_number:d}] Reading kpoint={k_index+1:d} band={band_index+1:d} eigenvalue from EIGENVAL.")
                eigen_text[i]=eigen_text[i].split( )

                # spin=1
                eigen_val_up=float(eigen_text[i][1])
                e_ef_up=eigen_val_up-ef
                dfdd_up=(1.0/kbT)*np.exp(e_ef_up/kbT)/(np.exp(e_ef_up/kbT)+1)/(np.exp(e_ef_up/kbT)+1)
                if dfdd_up < dfdd_threshold:
                    print(f"kpoint={k_index+1:d} band={band_index+1:d} spin=1: -dFDD={dfdd_up:.8e} is too low, ignore.")
                else:
                    eigen_up[eigen_val_up]=[k_index+1,band_index+1]

                #spin=2
                eigen_val_dw=float(eigen_text[i][2])
                e_ef_dw=eigen_val_dw-ef
                dfdd_dw=(1.0/kbT)*np.exp(e_ef_dw/kbT)/(np.exp(e_ef_dw/kbT)+1)/(np.exp(e_ef_dw/kbT)+1)
                if dfdd_dw < dfdd_threshold:
                    print(f"kpoint={k_index+1:d} band={band_index+1:d} spin=2: -dFDD={dfdd_dw:.8e} is too low, ignore.")
                else:
                    eigen_dw[eigen_val_dw]=[k_index+1,band_index+1]

                    


        eigen_val_up=sorted(eigen_up.keys())
        eigen_val_dw=sorted(eigen_dw.keys())


        # spin=1
        print('''\nStart calculating intergral of spin=1:\nIndex\tKpoint\tBand\tE-Ef/eV\t\t-dFDD\t\t\tweight''')
        for i in range(len(eigen_val_up)):
            band_index=eigen_up[eigen_val_up[i]][1]
            k_index=eigen_up[eigen_val_up[i]][0]
            tmp=subprocess.getstatusoutput(f"ls WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_UP.vasp.cube WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube")
            if "No such file or directory" in tmp[1]:
                #print("Running vaspkit for generating wavefunction file.")
                vaspkit_ini=open('vaspkit.ini','w')
                vaspkit_ini.write(f'51\n516\n{k_index:d}\n{band_index:d}\n')
                vaspkit_ini.close()
                tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
                if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
                    print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
                    exit()
                #else:
                    #print("---> Successed.")
                    #pass
            #else:
                #print("---> Wavefunction file exists, abord running vaspkit.")
                #pass

            #print(f"Reading WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_UP.vasp.cube")
            data, atoms = read_cube_data(f'WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_UP.vasp.cube')
            data=uniform(atoms,data)

            if intermediate_file_options==False:
                subprocess.getstatusoutput(f"rm WFN_SQUARED_B{band_index:04d}_K{k_index:04d}*.vasp.cube")

            e_ef=eigen_val_up[i]-ef
            dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
            data=data*dfdd
            if i==0:
                print(f"{i+1:d}/{len(eigen_val_up):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:f}")
                fs_up=data*kweight[k_index-1]
            else:
                #dE=eigen_val_up[i]-eigen_val_up[i-1]
                print(f"{i+1:d}/{len(eigen_val_up):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:f}")
                fs_up=fs_up+data*kweight[k_index-1]
        
        # write cube
        fs_up_file=open("LFS_UP.cube",'w')
        write_cube(fs_up_file,atoms,fs_up,[0.0,0.0,0.0],"Fermi_Softness_Spin=1")
        fs_up_file.close()


        # spin=2    
        print('''\nStart calculating intergral of spin=2:\nIndex\tKpoint\tBand\tE-Ef/eV\t\t-dFDD\t\t\tweight''')
        for i in range(len(eigen_val_dw)):
            band_index=eigen_dw[eigen_val_dw[i]][1]
            k_index=eigen_dw[eigen_val_dw[i]][0]
            tmp=subprocess.getstatusoutput(f"ls WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube")
            if "No such file or directory" in tmp[1]:
                #print("Running vaspkit for generating wavefunction file.")
                vaspkit_ini=open('vaspkit.ini','w')
                vaspkit_ini.write(f'51\n516\n{k_index:d}\n{band_index:d}\n')
                vaspkit_ini.close()
                tmp=subprocess.getstatusoutput("vaspkit < vaspkit.ini > vaspkit.log")
                if tmp[0] != 0 or "unknown" in tmp[1].lower() or "error" in tmp[1].lower():
                    print("**** !!!! Running vaspkit error! check the vaspkit.log !!!! ****")
                    exit()
                #else:
                    #print("---> Successed.")
                    #pass
            #else:
                #print("---> Wavefunction file exists, abord running vaspkit.")  
                #pass          
            #print(f"Reading WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube")
            data, atoms = read_cube_data(f'WFN_SQUARED_B{band_index:04d}_K{k_index:04d}_DW.vasp.cube')
            data=uniform(atoms,data)

            if intermediate_file_options==False:
                subprocess.getstatusoutput(f"rm WFN_SQUARED_B{band_index:04d}_K{k_index:04d}*.vasp.cube")

            e_ef=eigen_val_dw[i]-ef
            dfdd=(1.0/kbT)*np.exp(e_ef/kbT)/(np.exp(e_ef/kbT)+1)/(np.exp(e_ef/kbT)+1)
            data=data*dfdd
            if i==0:
                print(f"{i+1:d}/{len(eigen_val_dw):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:.6f}")
                fs_dw=data*kweight[k_index-1]
            else:
                #dE=eigen_val_dw[i]-eigen_val_dw[i-1]
                print(f"{i+1:d}/{len(eigen_val_dw):d}\t{k_index:d}\t{band_index:d}\t{e_ef:.6f}\t{dfdd:.8e}\t\t{kweight[k_index-1]:.6f}")
                fs_dw=fs_dw+data*kweight[k_index-1] 


        # write cube
        fs_dw_file=open("LFS_DW.cube",'w')
        write_cube(fs_dw_file,atoms,fs_dw,[0.0,0.0,0.0],"Fermi_Softness_Spin=2")
        fs_dw_file.close()
        #-------End:Local FS-----------------
    #----------End:spin polarization--------------


    #----------if spin error!---------------------
    else:
        print("**** !!!! Reading ISPIN from OUTCAR error !!!! ****")
        exit()


    #-----------save intermediate files-----------
    if intermediate_file_options==True:
        subprocess.getstatusoutput('mkdir WFNSQR DOS')
        subprocess.getstatusoutput('mv WFN_SQUARED* WFNSQR')
        subprocess.getstatusoutput('mv *DOS*.* DOS')
        subprocess.getstatusoutput('mv SELECTED_ORBITALS_LIST DOS')


    #-----------remove temp----------
    subprocess.getstatusoutput('rm vaspkit.ini vaspkit.log')

    #-----------print success--------
    print('\nThe calculation ends normally.')
    
#-------End:main------------



