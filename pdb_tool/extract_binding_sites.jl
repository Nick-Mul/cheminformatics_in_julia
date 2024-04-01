using PDBTools
#using Molecules
using Comonicon


function write_ligand(ligand)
    open(pwd() * "/ligand.xyz", "w") do f
        atom_count = length(ligand)
        #write(f, "$atom_count \n")
        #write(f, "ligand \n")
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


function write_binding_site(ligand, protein)
    active_site_atoms = Atom[]
    for residue in eachresidue(protein)
        if distance(residue, ligand) < 3.5
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
for i in active_site_atoms
    println(index(i))
end
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
    write_ligand(ligand)
    write_binding_site(ligand, protein)

    return
end

"""
read_pdb structre from pdb
# Args
- `pdb_file_name`: input file
- `ligand_name` : second argument (e.g. T3)
- `addH <arg>` : an option
"""
@cast function read_pdb(pdb_file::String, ligand_name::String; addH::Bool)
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
    ligand = select(model,"resname $ligand_name")
end
    write_ligand(ligand)
    write_binding_site(ligand, protein)
    return
end

"""
My main command.
"""
@main # declare the entry