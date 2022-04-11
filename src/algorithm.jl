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

uses MolecularGraph.jl instead of RDKit.jl
=#

function get_num_atoms(mol)
    atoms = MolecularGraph.atomcounter(mol)
    return sum(values(atoms))
end

function get_heavy_atom_count(mol)
    atoms = MolecularGraph.atomcounter(mol)
    atoms[:H] = 0
    return sum(values(atoms))
end

struct Fragmenter{S,O,F}
    name::Symbol
    n_max_fragmentations_to_find::Int
    n_atoms_cuttoff::Int
    match_hydrogens::Bool
    scheme::S
    scheme_order::O
    function_to_choose_fragmentation::F
    scheme_pattern_lookup::Dict{String,Any}
    scheme_group_number_lookup::Dict{String,String}
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

function Fragmenter(scheme,
    scheme_order = nothing;
    match_hydrogens=false,
    algorithm::Symbol = :combined,
    n_atoms_cutoff = -1,
    function_to_choose_fragmentation = nothing,
    n_max_fragmentations_to_find = -1
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
        if n_atoms_cuttoff == -1 
            throw(error("n_atoms_cuttoff needs to be specified for complete or combined algorithms"))
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
        if frag.match_hydrogens
            molecule = MolecularGraph.addhydrogens(molecule)
        end
        return fragment(molecule,frag)
    catch
        throw(error("$smiles is not a valid SMILES"))
    end
end

#the real algorithm starts here
function fragment(molecule::MolecularGraph.GraphMol,frag::Fragmenter)
    
end