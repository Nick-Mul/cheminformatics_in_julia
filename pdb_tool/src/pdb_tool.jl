
using PDBTools
using Molecules


ligand_name = "T3"

model = wget("1BSX")
protein = select(model,"protein")
ligand = select(model,"resname $ligand_name and chain A")

##### write ligand as xyz ######
function write_ligand(ligand)
    open(pwd() * "/ligand.xyz", "w") do f
        atom_count = length(ligand)
        #write(f, "$atom_count \n")
        #write(f, "ligand \n")
        for i in ligand
            atom = element(i)
            x, y, z = position(i)
            write(f, "$atom $x $y $z \n")
        end
    end
end

function write_binding_site(ligand, protein)
    binding_site = []
    binding_site_residues = []
    charge = 0
    for atom in ligand
        lig, protein_atom, dist = (closest(atom, protein))
        push!(binding_site, (protein[protein_atom]))
    end
    for i in binding_site
        push!(binding_site_residues, i.resnum)
    end
    binding_site_residues_clean = sort(union(binding_site_residues))
    res = [select(protein, "resnum = $i and chain = A") for i in binding_site_residues_clean]
    open(pwd() * "/binding_site.xyz", "w") do f
        for i in res
            for j in i
                x, y, z = position(j)
                atom = element(j)
                write(f, "$atom $x $y $z \n")
            end
        end
    end
end
write_ligand(ligand)
write_binding_site(ligand, protein)


s = open("ligand.xyz") do file
    read(file, String)
end

ligand_parsed = Molecule(s)

println(ligand_parsed.charge)
println(Molecules.nuclear_dipole(ligand_parsed))
println(Molecules.center_of_mass(ligand_parsed))
# print(Molecules.Symmetry.find_point_group(ligand_parsed))