using MolecularGraph
 

mol = try smilestomol("CCCC") catch err end


LIST_HYBRIDE = [:sp, :sp2, :sp3, :none]

custom_mtd!(mol::SimpleMolGraph, input_data
    ) = set_cache!(mol, :v_custom, input_data)

# make more robust with 
#if item1 ∉ permitted_list
#        id = permitted_list[end]

function encode(list1, list2)
    """
    Encode a list of items from `list1` based on their positions in `list2`.

    This function takes two lists as input, `list1` and `list2`, and returns an array of
    floating-point numbers representing the positions of items from `list1` in `list2`.

    Parameters:
    - list1 (Array): The list of items to be encoded.
    - list2 (Array): The list in which to find the positions of items from `list1`.

    Returns:
    - encoded_positions (Array{Float64, 1}): An array of floating-point numbers representing
      the positions of items from `list1` in `list2`. If an item from `list1` is not found
      in `list2`, its position is set to `NaN`.

    Example:
    ```julia
    list1 = ["apple", "banana", "cherry"]
    list2 = ["banana", "cherry", "apple"]
    encoded_positions = encode(list1, list2)
    println(encoded_positions)
    # Output: [3.0, 2.0, 1.0]
    ```
    """
    positions = []
    #Find the position of item1 in list2
    for item1 in list1
        index = (findfirst(x -> x == item1, list2))
        push!(positions, index)
end
return  convert(Array{Float64,1}, positions)
end

function mol2vec(mol)
    """
    Generate a feature vector for a molecule.

    This function takes a molecule object `mol` as input and returns a feature vector
    representing various properties of its atoms.

    Parameters:
    - mol: A molecule object (e.g., from a chemical informatics library like RDKit).

    Returns:
    - atom_feature_vector (Array{Float64, 1}): A feature vector representing atom properties
      of the given molecule. The vector includes information such as atomic charge,
      aromaticity, apparent valence, hybridization, implicit hydrogens, and more.
"""
    atom_feature_vector = vcat(charge(mol),
        is_aromatic(mol) .* 1,
        apparent_valence(mol),
        encode(hybridization(mol), LIST_HYBRIDE),
        encode(hybridization_delocalized(mol), LIST_HYBRIDE),
        implicit_hydrogens(mol),
        is_aromatic(mol) .* 1,
        is_in_ring(mol) .* 1,
    )
return atom_feature_vector
end

print(mol2vec(mol))

println(mol)


custom_mtd!(mol, [999,999,9999999])
println(mol)

function is_custom(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_custom) && return get_cache(mol, :v_custom)
    return is_custom(mol.graph)
end

println(is_custom(mol))

permitted_list_of_atoms =  ["C","N","O","S","F","Si","P","Cl","Br","Mg","Na","Ca","Fe","As","Al","I", "B","V","K","Tl","Yb","Sb","Sn","Ag","Pd","Co","Se","Ti","Zn", "Li","Ge","Cu","Au","Ni","Cd","In","Mn","Zr","Cr","Pt","Hg","Pb","Unknown"]
degree_list = [0, 1, 2, 3, 4, "MoreThanFour"]

n_node_features = length(permitted_list_of_atoms)+length(degree_list)
print("lenghth of vector ", n_node_features)
