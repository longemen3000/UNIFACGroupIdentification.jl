#=
Julia port of:

Class for fragmenting molecules into molecular subgroups
MIT License
Copyright (C) 2019, Simon Mueller <simon.mueller@tuhh.de>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

#=

=#
function get_num_atoms(mol)
    atoms = MolecularGraph.atomcounter(mol)
    return sum(values(atoms))
end

function get_atoms(mol)
    atoms = MolecularGraph.atomcounter(mol)
    return keys(atoms)
end

function get_heavy_atom_count(mol)
    atoms = MolecularGraph.atomcounter(mol)
    atoms[:H] = 0
    return sum(values(atoms))
end

struct Fragmenter{S,O,F}
    name::Symbol
    algorithm::Symbol
    n_max_fragmentations_to_find::Int
    n_atoms_cuttoff::Int
    match_hydrogens::Bool
    scheme::S
    scheme_order::O
    function_to_choose_fragmentation::F
    scheme_pattern_lookup::Dict{String,Any}
    scheme_group_number_lookup::Dict{String,String}
end

function Base.show(io::IO,::MIME"text/plain",frag::Fragmenter)
    space = "  "
    alg = frag.algorithm
    println(io," Fragmenter(",frag.name,"): ")
    println(io,space,"algorithm: ",alg)

    if alg != :simple
        println(io,space,"max fragmentations to find: ",frag.n_max_fragmentations_to_find)
        println(io,space,"atoms cutoff: ",frag.n_atoms_cuttoff)
    end

    print(io,space,"scheme: ",frag.scheme)
    if frag.scheme_order !== nothing && alg != :simple
        println(io)
        print(io,space,"scheme order: ",frag.scheme_order)
    end
end

function Base.show(io::IO,frag::Fragmenter)
    comma = ", "
    alg = frag.algorithm
    print(io,"Fragmenter(")
    print(io,frag.name,comma)
    print(io,"algorithm = ",frag.algorithm,comma)

    if alg != :simple
        print(io,"n_max_fragmentations_to_find = ",frag.n_max_fragmentations_to_find,comma)
        print(io,"n_atoms_cutoff = ",frag.n_atoms_cuttoff)
    end

    print(io,"scheme = ",frag.scheme)
    if frag.scheme_order !== nothing && alg != :simple
        print(io,comma,"scheme_order = ",frag.scheme_order)
    end
    print(io,")")
end

function mol(smiles)
    return MolecularGraph.smilestomol(smiles)
end

function get_substruct_matches(mol_searched_in,mol_searched_for)
    return MolecularGraph.structmatches(mol_searched_in,mol_searched_for,:substruct)
end

function get_atom_with_idx(i)
    return 1
end

function get_idx(m)
    return 1
end

function get_neighbors(m)
    return 1
end

function get_mol_frags(m)
    mol = MolecularGraph.removehydrogens(m)
    return MolecularGraph.Graph.connectedcomponents(mol)
end

function get_substruct_matches(frag::Fragmenter,mol_searched_for, mol_searched_in, atomIdxs_to_which_new_matches_have_to_be_adjacent)
    valid_matches = []
    if get_num_atoms(mol_searched_in) >= get_num_atoms(mol_searched_for)
        matches = get_substruct_matches(mol_searched_in,mol_searched_for)
        add_this_match = false
        for match in matches #if is empty, it will not iterate
            add_this_match = true
            length(atomIdxs_to_which_new_matches_have_to_be_adjacent) > 0 && (add_this_match = true)
            for i in match
                for neighbor in (mol_searched_in |> get_atom_with_idx |> get_neighbors)
                    if (neighbor |> get_idx) in atomIdxs_to_which_new_matches_have_to_be_adjacent
                        add_this_match = true
                    end
                end
            end
            if add_this_match
                append!(valid_matches,match)
            end
        end
    end
    return valid_matches
end

function Fragmenter(scheme;
    scheme_order = nothing,
    match_hydrogens=true,
    algorithm::Symbol = :combined,
    n_atoms_cutoff = -1,
    function_to_choose_fragmentation = first,
    n_max_fragmentations_to_find = -1,
    name = :auto)

    if algorithm ∉ (:simple,:complete,:combined)
        throw(error("Algorithm must be either `:simple` ,`:complete` or `:combined`"))
    end

    if algorithm == :simple
        n_max_fragmentations_to_find != -1 && throw(
            error("Setting `n_max_fragmentations_to_find` 
            only makes sense with complete or combined algorithm`"))
    end

    if algorithm ∈ (:complete,:combined)
        if n_atoms_cutoff == -1 
            throw(error("n_atoms_cutoff needs to be specified for complete or combined algorithms"))
        end

        if isnothing(function_to_choose_fragmentation)
            throw(error("`function_to_choose_fragmentation` needs to be specified for complete or combined algorithms"))
        end

        if max(0,n_max_fragmentations_to_find) < 1
            throw(error("n_max_fragmentations_to_find has to be 1 or higher."))
        end
    end

    if algorithm ∈ (:complete,:combined)
        if isnothing(scheme_order)
            scheme_order = first.(scheme)
            #im supposing scheme is a vector of 2-tuples or a vector of pairs.
            #python dicts are ordered, we can't count on that on julia.
        end
    end
    scheme_group_number_lookup = Dict{String,String}()
    scheme_pattern_lookup = Dict{String,Any}()

    for (group_number, list_SMARTS) in scheme #works with Vector{Pair{String, Any}}:
        #more than group number is a group text.
        if list_SMARTS isa String
            list_SMARTS = (list_SMARTS,) #seems to accept arbitrary SMARTS per each match, nice
        end
        for smarts in list_SMARTS     
            try
                querymol = MolecularGraph.smartstomol(smarts)
                scheme_pattern_lookup[smarts] = querymol
                scheme_group_number_lookup[smarts] = group_number
            catch
                throw(error("$smarts is not a valid SMARTS"))
            end
        end
    end

    return Fragmenter(:auto,
        algorithm,
        n_max_fragmentations_to_find,
        n_atoms_cutoff,
        match_hydrogens,
        scheme,
        scheme_order,
        function_to_choose_fragmentation,
        scheme_pattern_lookup,
        scheme_group_number_lookup)
end

function fragment(smiles::String,frag::Fragmenter)
    try
        molecule = MolecularGraph.smilestomol(smiles)
    catch
        throw(error("$smiles is not a valid SMILES"))
    end
    if frag.match_hydrogens
        molecule = MolecularGraph.addhydrogens(molecule)
    end
    return fragment(molecule,frag)
end

#the real algorithm starts here
function fragment(molecule::MolecularGraph.GraphMol,fragmenter::Fragmenter)
    sucess = []
    fragmentation = Dict{Int,Any}()
    fragmentation_matches = Dict{Int,Any}()
    frags = get_mol_frags(molecule)
    for frag in frags
        this_mol_fragmentation, this_mol_success = __get_fragmentation(molecule,frag,fragmenter)
        for (SMARTS, matches) in pairs(this_mol_fragmentation)
            group_number = fragmenter.scheme_group_number_lookup[SMARTS]
            if group number ∉ fragmentation
                fragmentation[group_number] = 0
                fragmentation_matches[group_number] = []
            end
            fragmentation[group_number] += length(matches)
            append!(fragmentation_matches[group_number],matches)   
        end
        append!(success,this_mol_success)
    end
    return fragmentation,all(success),fragmentation_matches    
end

function fragment_complete(molecule,frag::Fragmenter)
    try
        molecule = MolecularGraph.smilestomol(smiles)
        if frag.match_hydrogens
            molecule = MolecularGraph.addhydrogens(molecule)
        end
        return fragment_complete(molecule,frag)
    catch
        throw(error("$smiles is not a valid SMILES"))
    end
end

function fragment_complete(molecule::MolecularGraph.GraphMol,frag::Fragmenter)
    frags = get_mol_frags(molecule)
    if !isone(length(frags))
        throw(error("`fragment_complete` does not accept multifragment molecules"))
    end

    temp_fragmentations, success = __complete_fragmentation(molecule)
    fragmentations = []
    fragmentations_matches = []
    for temp_fragmentation in temp_fragmentations
        fragmentation = Dict{Int,Any}()
        fragmentation_matches = Dict{Int,Any}()
        for (SMARTS, matches) in pairs(temp_fragmentation)
            group_number = frag.scheme_group_number_lookup[SMARTS]
            fragmentation[group_number] = length(matches)
            fragmentation_matches[group_number] = matches
        end
        push!(fragmentations,fragmentation)
        push!(fragmentations_matches,fragmentation_matches)
    end
    return fragmentations, success, fragmentations_matches
end

function __get_fragmentation(molecule,frag,fragmenter::Fragmenter)
   success = false
    if fragmenter.algorithm in (:simple,:combined)
        fragmentation, success = __simple_fragmentation(molecule,frag,fragmenter)
        success && (return fragmentation, success)
    end

    if frag.algorithm in (:combined,:complete)
        fragmentation, success = __complete_fragmentation(molecule,frag,fragmenter)
        success && (return frag.function_to_choose_fragmentation(fragmentation), success)
    end
end

function __simple_fragmentation(molecule,frag,fragmenter::Fragmenter)
    if fragmenter.match_hydrogens
        target_atom_count = length(get_atoms(molecule))
    else
        target_atom_count = length(get_heavy_atoms(molecule))
    end
    success = false
    fragmentation_so_far = Dict{String,Any}()
    fragmentation_so_far, atomIdxs_included_in_fragmentation = __search_non_overlapping_solution(fragmenter, molecule, deepcopy(fragmentation_so_far), Set(), Set())
    success = length(atomIdxs_included_in_fragmentation) == target_atom_count
    level = 1
    while !success
        fragmentation_so_far , atomIdxs_included_in_fragmentation_so_far = fragmenter.__clean_molecule_surrounding_unmatched_atoms(molecule, fragmentation, atomIdxs_included_in_fragmentation, level)
        level +=1
        len = length(atomIdxs_included_in_fragmentation_so_far)
        fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far = fragmenter.__search_non_overlapping_solution(molecule, fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far)
        
        success = (len == target_atom_count)
        iszero(len) && break
        success && break
    end
    return fragmentation_so_far, success
end

