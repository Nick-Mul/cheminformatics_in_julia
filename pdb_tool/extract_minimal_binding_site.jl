using PDBTools
using Comonicon

TYR_list = [
    "CB",
    "HB1",
    "HB2",
    "CG",
    "CD1",
    "HD1",
    "CE1",
    "HE1",
    "CD2",
    "HD2",
    "CE2",
    "HE2",
    "CZ",
    "OH",
    "HH"
]

ARG_list = [
    "CD",
    "HD1",
    "HD2",
    "NE",
    "HE",
    "CZ",
    "NH1",
    "1HH1",
    "2HH1",
    "NH2",
    "1HH2",
    "2HH2"
]

TRP_list = [
    "CB",
    "HB1",
    "HB2",
    "CG",
    "CD2",
    "CE2",
    "CE3",
    "HE3",
    "CD1",
    "HD1",
    "NE1",
    "HE1",
    "CZ2",
    "HZ2",
    "CZ3",
    "HZ3",
    "CH2",
    "HH2"
]

LYS_list = ["CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3"]

function write_ligand(ligand)
    open(pwd() * "/ligand.xyz", "w") do f
        atom_count = length(ligand)
        write(f, "$atom_count \n")
        write(f, "ligand \n")
        for j in ligand
            println(index(j))
        end
        for i in ligand
            atom = element(i)
            x, y, z = position(i)
            write(f, "$atom $x $y $z \n")
        end
    end
end

function define_active_site_atoms(protein, ligand)
    active_site_atoms = Atom[]
    for residue in eachresidue(protein)
        if distance(residue, ligand) < 3.5
            append!(active_site_atoms, atom for atom in residue)
        end
    end
    return active_site_atoms
end


function write_binding_site(ligand, protein)
    active_site_atoms = Atom[]
    for residue in eachresidue(protein)
        if distance(residue, ligand) < 8
            append!(active_site_atoms, atom for atom in residue)
        end
    end
    open(pwd() * "/binding_site.xyz", "w") do f
        for i in active_site_atoms
            x, y, z = position(i)
            atom = element(i)
            write(f, "$atom $x $y $z \n")
        end
    end
    writePDB(active_site_atoms,"file.pdb")
for i in active_site_atoms
    println(index(i))
end
end



function trucate_residues(residues, protein, ligand)
    """
    trucates residue and returns atom_list
    """
    atom_list = []
    for atom in protein
        if closest(atom, ligand)[3] < 3.0
            push!(atom_list, atom.index)
        end
    end
    for residue in residues
        for atoms in residue
            if atoms.resname == "TYR" && atoms.name in TYR_list
                #println(atoms.resname ," " , atoms.name, " ", element(atoms))
                push!(atom_list, atoms.index)
            end
            if atoms.resname == "ARG" && atoms.name in ARG_list
                #println(atoms.resname ," " , atoms.name, " ", element(atoms))
                push!(atom_list, atoms.index)
            end
            if atoms.resname == "TRP" && atoms.name in TRP_list
                #println(atoms.resname ," " , atoms.name, " ", element(atoms))
                push!(atom_list, atoms.index)
            end
        end
    end
    return sort(unique(atom_list))
end

function get_close_atoms(protein, ligand)
    atom_list = []
    for atom in protein
        if closest(atom, ligand)[3] < 3.0
            push!(atom_list, atom.index)
        end
    end
    return atom_list
end

"""
get structre from pdb
# Args
- `pdb_file_name`: first argument (e.g. 1BSX)
- `ligand_name` : second argument (e.g. T3)
"""
@cast function get_pdb(pdb_file::String, ligand_name::String)
    #T3
    #1BSX
    model = wget("$pdb_file", "chain A")
    protein = select(model,"protein")
    ligand = select(model,"resname $ligand_name")
    #write_ligand(ligand)
    #write_binding_site(ligand, protein)
    #
    return
end

"""
read_pdb structre from pdb
# Args
- `pdb_file_name`: input file
- `ligand_name` : second argument (e.g. T3)
- `addH <arg>` : add hydrogens if not present in pdb (uses obabel), default is false
- `chain <arg>` : ligand chain (default chain A)
- `waters <arg>` : include water in receptor
- `cofactor <arg>` : additioan cofactor to include (e.g. PLP)
"""
@cast function read_pdb(pdb_file::String, ligand_name::String; addH::Bool=false)
    #T3
    #1BSX
    if addH true
        println("adding Hydrogens with o")
        run(`obabel $pdb_file -h -O output_with_hydrogens.pdb`)
        model = readPDB("output_with_hydrogens.pdb")
        protein = select(model,"protein")
        ligand = select(model,"resname $ligand_name")
else
    #run(`obabel path_to_your_pdb_file.pdb -h -O output_with_hydrogens.pdb`)
    model = readPDB("$pdb_file")
    protein = select(model,"protein")
    ligand = select(model,"resname $ligand_name and chain G")
end
    write_ligand(ligand)
    write_binding_site(ligand, protein)
    active_site_atoms = define_active_site_atoms(protein, ligand)
    residues = (eachresidue(active_site_atoms))
    at_list = trucate_residues(residues, protein, ligand)
    get_close_atoms(protein, ligand)
    print(at_list)
    for at in at_list
        x = select(protein, "index $at")
        for i in x
            println(position(i))
        end
    end
    return
end

"""
My main command.
"""
@main # declare the entry


