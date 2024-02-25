
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
    active_site = []
    for r in eachresidue(protein)
        if distance(r, ligand) < 3.5
            append!(active_site, at for at in r)
        end
    end
    open(pwd() * "/binding_site.xyz", "w") do f
        for i in active_site
            x, y, z = position(i)
            atom = element(i)
            write(f, "$atom $x $y $z \n")
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