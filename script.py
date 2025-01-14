import cmath
import numpy as np
import math


def read_vector(name):
    vector = []
    tmp = [0,0,0]
    vector_file = open(name,"r")
    lines = vector_file.readlines()
    q_point = [float(lines[0].split()[0]), float(lines[0].split()[1]), float(lines[0].split()[2])]
    for i in range(2,len(lines),3):
        for j in range(0,3):    tmp[j] = ( complex(float(lines[i+j].split()[0]),float(lines[i+j].split()[1])) )
        vector.append([tmp[0],tmp[1],tmp[2]])
    vector_file.close()
    return q_point, vector 



def read_poscar(name):
    cell = []
    atomic_positions = []
    poscar_file = open(name,"r")
    poscar_file.readline()
    s = float(poscar_file.readline())
    for i in range(0,3):
        line = poscar_file.readline()
        cell.append([float(line.split()[0])*s, float(line.split()[1])*s, float(line.split()[2])*s])
    atom_types = poscar_file.readline().split()
    atoms_type_numbers_degen =  poscar_file.readline().split()
    total_n_atoms = 0
    for i in range(0,len(atoms_type_numbers_degen)):
        atoms_type_numbers_degen[i] = int(atoms_type_numbers_degen[i])
        total_n_atoms = total_n_atoms + atoms_type_numbers_degen[i]
    coordinate_type = str(poscar_file.readline())
    for i in range(0,total_n_atoms):
        line = poscar_file.readline()
        atomic_positions.append([ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ])
    poscar_file.close()
    return s, coordinate_type, cell, atomic_positions, atoms_type_numbers_degen, atom_types

def write_poscar(name,s, coordinate_type, cell, atomic_positions, atoms_type_numbers_degen, atom_types):
    poscar_file = open(name,"w")
    poscar_file.write("File generate by Francesco\n")
    poscar_file.write(str(s) + "\n")
    for i in cell:
        poscar_file.write(str(i[0]) + "   "  + str(i[1]) +  "   "  + str(i[2]) + "\n")
    for i in atom_types:
        poscar_file.write(str(i) + "   "  )
    poscar_file.write("\n")
    for i in atoms_type_numbers_degen:
        poscar_file.write(str(i) + "   "  )
    poscar_file.write("\n")
    poscar_file.write(coordinate_type)
    for i in atomic_positions:
        poscar_file.write(str(i[0]) + "   "  + str(i[1]) +  "   "  + str(i[2]) + "\n")
    poscar_file.close()

def write_scf(name,s, coordinate_type, cell, atomic_positions, atoms_type_numbers_degen, atom_type):
    scf_file = open(name,"w")
    header_file = open("HEADER",'r')
    header_data = header_file.read()
    header_file.close()
    scf_file.write(header_data)
    scf_file.write("\nCELL_PARAMETERS {angstrom}\n")
    for i in cell:
        scf_file.write(str(i[0]) + "   "  + str(i[1]) +  "   "  + str(i[2]) +  "   "  + "\n")
    scf_file.write("\nATOMIC_POSITIONS {" + str(coordinate_type) + "}\n")
    x = 0
    for j in range(0,len(atoms_type_numbers_degen)):
        for i in range(0, atoms_type_numbers_degen[j]):
            scf_file.write(str(atom_type[j]) + "   "  + str(atomic_positions[x][0]) + "   "  + str(atomic_positions[x][1]) +  "   "  + str(atomic_positions[x][2]) +  "\n")
            x += 1
    scf_file.close()


def crystal_to_cartesian(cell,positions):
    tmp_positions = []
    tmp = [0,0,0]
    for i in range(0,len(positions)):
        tmp = [0,0,0]
        for j in range(0,3):
            for k in range(0,3):
                tmp[j] += cell[k][j]*positions[i][k]
        tmp_positions.append(tmp)
    return tmp_positions




print("Unit cell poscar will be read from POSCAR.")
print("q pint and vectors will be read from VECTOR.")
xx,yy,zz = input("Enter a truple for the supercell size: ").split()
supercell = [int(xx),int(yy),int(zz)]
factor = float(input("Enter scaling factor: "))






s, coordinate_type, cell, positions, atoms_type_numbers_degen, atom_types = read_poscar("POSCAR")
#write_poscar("TEST2",s, coordinate_type, cell, positions, atoms_type_numbers_degen, atom_types)
q,vector = read_vector("VECTOR")


new_cell = [[0 for _ in range(3)] for _ in range(3)]



if coordinate_type.strip() == "Direct": positions = crystal_to_cartesian(cell,positions)



for i in range(0,3):
    for j in range(0,3):
        new_cell[i][j] = cell[i][j]*supercell[i]


tmp = [0,0,0]
new_positions = []
new_atoms_type_numbers_degen = atoms_type_numbers_degen.copy()
new_atom_types = atom_types.copy()
for i in  range(0,len(new_atoms_type_numbers_degen)): new_atoms_type_numbers_degen[i]*=supercell[0]*supercell[1]*supercell[2]

for i in range(-10,10,1):
    step = i/factor
    for j in range(0,len(positions)):
        for x in range(0,supercell[0]):
            for y in range(0,supercell[1]):
                for z in range(0,supercell[2]):
                    for k in range(0,3):
                        vec = vector[j][k]*cmath.exp(2j*math.pi*np.dot(q,[x,y,z]))
                        tmp[k] = positions[j][k] + cell[0][k]*x + cell[1][k]*y + cell[2][k]*z  + step * vec.real/math.sqrt(supercell[0]*supercell[1]*supercell[2])
                    new_positions.append([tmp[0],tmp[1],tmp[2]])
    write_poscar("POSCAR_"+str(step)+".vasp",s,"Cartesian\n", new_cell, new_positions, new_atoms_type_numbers_degen, new_atom_types)
    write_scf("scf_"+str(step)+".in",s, "angstrom", new_cell, new_positions, new_atoms_type_numbers_degen, new_atom_types)
    new_positions = []
print(new_cell)
