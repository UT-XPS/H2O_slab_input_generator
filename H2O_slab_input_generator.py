import math

## Plan.
##
## Denote the positions of the water molecules within the unit cell. Enter this data manually.
## The list of water molecules should be a list of lists, where each sublist is one water molecule.
## Each sublist contains the following elements.
## An index
## A sub-sublist of four elements: Xpos, Ypos, Zpos, identity=O
## A sub-sublist of four elements: Xpos, Ypos, Zpos, identity=H
## A sub-sublist of four elements: Xpos, Ypos, Zpos, identity=H
##
## From the list of water molecules, form lists of atoms, bonds and angles in a suitable format
## for LAMMPS output.

## a - size of repeating unit along x
## b - size of repeating unit along y
## c - size of repeating unit along z
## constant - should be the same for all inputs: ensures that the water molecules are generated with the right bond length (1 Angstrom)
## l - number of repeating units along x
## m - number of repeating units along y
## n - number of repeating units along z
## z_gap - only affects the LAMMPS input file: adjusts the zmax parameter to leave extra space
## O_charge - in units of elementary charge
## H_charge - in units of elementary charge
## O_mass - in atomic mass units
## H_mass - in atomic mass units
def make_water_LAMMPS_input(a,b,c,constant,l,m,n,z_gap,O_charge,H_charge,O_mass,H_mass):
    ice_cell = make_pseudocubic_water_cell(a,b,c,constant)
    supercell = replicate_l_m_n_times(ice_cell,a,b,c,l,m,n)
    atoms = get_atoms(supercell,O_charge,H_charge)
    bonds = get_bonds(supercell)
    angles = get_angles(supercell)
    output(create_LAMMPS_output(atoms,bonds,angles,a,b,c,l,m,n,z_gap,O_mass,H_mass))

def make_pseudocubic_water_cell(a,b,c,constant):
    water_unit_cell = [[1,[0.0,0.0,0.25,"O"],[-0.125,-0.125,0.38683,"H"],[0.125,0.125,0.38683,"H"]],
                       [2,[0.0,0.5,0.5,"O"],[-0.125,0.375,0.63683,"H"],[0.125,0.625,0.63683,"H"]],
                       [3,[0.5,0.0,0.5,"O"],[0.375,-0.125,0.63683,"H"],[0.625,0.125,0.63683,"H"]],
                       [4,[0.5,0.5,0.25,"O"],[0.375,0.375,0.38683,"H"],[0.625,0.625,0.38683,"H"]],
                       [5,[0.75,0.25,0.75,"O"],[0.625,0.125,0.88683,"H"],[0.875,0.375,0.88683,"H"]],
                       [6,[0.75,0.75,0.25,"O"],[0.625,0.625,0.38683,"H"],[0.875,0.875,0.38683,"H"]],
                       [7,[0.25,0.25,0.25,"O"],[0.125,0.125,0.38683,"H"],[0.375,0.375,0.38683,"H"]],
                       [8,[0.25,0.75,0.75,"O"],[0.125,0.625,0.88683,"H"],[0.375,0.875,0.88683,"H"]]]
    cell_in_Cart_lattice_ref = []
    for element in water_unit_cell:
        cell_in_Cart_lattice_ref += [convert_to_Cart_lattice_ref(element,constant)]
    strained_cell_in_Cart_lattice_ref = strain(cell_in_Cart_lattice_ref,a,b,c,constant)
    return strained_cell_in_Cart_lattice_ref

def convert_to_Cart_lattice_ref(element,constant):
    new_element = []
    new_element += [element[0]]
    for sub_element in element[1:]:
        new_sub_0 = sub_element[0]*constant
        new_sub_1 = sub_element[1]*constant
        new_sub_2 = sub_element[2]*constant
        new_sub_3 = sub_element[3]
        new_element += [[new_sub_0,new_sub_1,new_sub_2,new_sub_3]]
    return new_element

def strain(cell_in_Cart_ref,a,b,c,constant):
    a_ratio = a/constant
    b_ratio = b/constant
    c_ratio = c/constant
    strained_cell = []
    for molecule in cell_in_Cart_ref:
        start_O = molecule[1]
        end_O = [molecule[1][0]*a_ratio,molecule[1][1]*b_ratio,molecule[1][2]*c_ratio,"O"]
        translator = [end_O[0]-start_O[0],end_O[1]-start_O[1],end_O[2]-start_O[2]]
        end_H1 = translate(molecule[2],translator)
        end_H2 = translate(molecule[3],translator)
        moved_molecule = [molecule[0],end_O,end_H1,end_H2]
        strained_cell += [moved_molecule]
    return strained_cell

def replicate_l_m_n_times(cell_in_Cart_ref,a,b,c,l,m,n):
    supercell = []
    cell_counter = 0
    for i in range(l):
        for j in range(m):
            for k in range(n):
                translator = [a*i,b*j,c*k]
                for molecule in cell_in_Cart_ref:
                    new_0 = molecule[0]+cell_counter*8
                    new_1 = translate(molecule[1],translator)
                    new_2 = translate(molecule[2],translator)
                    new_3 = translate(molecule[3],translator)
                    supercell += [[new_0,new_1,new_2,new_3]]
                cell_counter += 1
    return supercell

def translate(atom,translator):
    moved_x = atom[0]+translator[0]
    moved_y = atom[1]+translator[1]
    moved_z = atom[2]+translator[2]
    name = atom[3]
    moved_atom = [moved_x,moved_y,moved_z,name]
    return moved_atom

def display(supercell,a,b,c,l,m,n):
    for molecule in supercell:
        H_iter = "a"
        for atom in molecule[1:]:
            if atom[3] == "O":
                print("    O"+str(molecule[0])+"    O    1.0000    "+str(atom[0]/(a*l))+"    "+str(atom[1]/(b*m))+"    "+str(atom[2]/(c*n)))
            if atom[3] == "H":
                print("    H"+str(molecule[0])+H_iter+"    H    1.0000    "+str(atom[0]/(a*l))+"    "+str(atom[1]/(b*m))+"    "+str(atom[2]/(c*n)))
                H_iter = "b"

def get_atoms(supercell,O_charge,H_charge):
    atoms = []
    for molecule in supercell:
        mol_no = molecule[0]
        atom_index_in_molecule = 1
        for atom in molecule[1:]:
            atom_no = atom_index_in_molecule + (mol_no - 1)*3
            x_pos = atom[0]
            y_pos = atom[1]
            z_pos = atom[2]
            identity = atom[3]
            if identity == "O":
                charge = O_charge
                atom_type = 1
            elif identity == "H":
                charge = H_charge
                atom_type = 2
            atoms += [[atom_no,mol_no,atom_type,charge,x_pos,y_pos,z_pos]]
            atom_index_in_molecule += 1
    return atoms

def get_bonds(supercell):
    bonds = []
    for molecule in supercell:
        mol_no = molecule[0]
        O_no = 1 + (mol_no - 1)*3
        H1_no = 2 + (mol_no - 1)*3
        H2_no = 3 + (mol_no - 1)*3
        bond_1_no = 1 + (mol_no - 1)*2
        bond_2_no = 2 + (mol_no - 1)*2
        bonds += [[bond_1_no,1,O_no,H1_no],[bond_2_no,1,O_no,H2_no]]
    return bonds

def get_angles(supercell):
    angles = []
    for molecule in supercell:
        mol_no = molecule[0]
        O_no = 1 + (mol_no - 1)*3
        H1_no = 2 + (mol_no - 1)*3
        H2_no = 3 + (mol_no - 1)*3
        angle_no = mol_no
        angles += [[angle_no,1,H1_no,O_no,H2_no]]
    return angles

def create_LAMMPS_output(atoms,bonds,angles,a,b,c,l,m,n,z_gap,O_mass,H_mass):
    output = []
    output += ["LAMMPS H2O ice input file"+"\n"]
    output += ["\t"+str(len(atoms))+"\t"+"atoms"+"\n"]
    output += ["\t"+str(len(bonds))+"\t"+"bonds"+"\n"]
    output += ["\t"+str(len(angles))+"\t"+"angles"+"\n"]
    output += ["\t"+"0"+"\t"+"dihedrals"+"\n"]
    output += ["\t"+"0"+"\t"+"impropers"+"\n"]
    output += ["\t"+"2"+"\t"+"atom types"+"\n"]
    output += ["\t"+"1"+"\t"+"bond types"+"\n"]
    output += ["\t"+"1"+"\t"+"angle types"+"\n"]
    output += ["\t"+"0"+"\t"+"dihedral types"+"\n"]
    output += ["\t"+"0"+"\t"+"improper types"+"\n"]
    output += ["\n"]
    output += ['{:f}'.format(0)+"\t"+'{:f}'.format(a*l)+"\t"+"xlo xhi"+"\n"]
    output += ['{:f}'.format(0)+"\t"+'{:f}'.format(b*m)+"\t"+"ylo yhi"+"\n"]
    output += ['{:f}'.format(0)+"\t"+'{:f}'.format((c*n)+z_gap)+"\t"+"zlo zhi"+"\n"]
    output += ["\n"]
    output += ["Masses"+"\n"]
    output += ["\n"]
    output += ["\t"+"1"+"\t"+'{:f}'.format(O_mass)+"\n"]
    output += ["\t"+"2"+"\t"+'{:f}'.format(H_mass)+"\n"]
    for entry_type in [["Atoms",atoms],["Bonds",bonds],["Angles",angles]]:
        output += ["\n"]
        output += [entry_type[0]+"\n"]
        output += ["\n"]
        for entry in entry_type[1]:
            line = ""
            for element in entry:
                line += "\t"
                if type(element) == int:
                    line += str(element)
                if type(element) == float:
                    line += '{:f}'.format(element)
            line += "\n"
            output += [line]
    return output

def output(stringlist):
    output = open("water.txt","w")
    output.flush()
    output.writelines(stringlist)
    output.close()

constant = 4.28190841155
#standard# a = 6.2086589795683516
#standard# b = 6.2086589795683516
#standard# c = 6.2086589795683516
a = 6.24
b = 6.24
c = 6.24
l = 10
m = 10
n = 10
O_charge = -1.1128
H_charge = 0.5564
z_gap = 0.0
O_mass = 15.9994
H_mass = 1.008
make_water_LAMMPS_input(a,b,c,constant,l,m,n,z_gap,O_charge,H_charge,O_mass,H_mass)
