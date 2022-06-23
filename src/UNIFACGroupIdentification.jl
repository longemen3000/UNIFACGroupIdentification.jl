module UNIFACGroupIdentification

using MolecularGraph
using RDKitMinimalLib

#MolFromSmarts: yes
#MolFromSmiles: yes
#AddHs : yes
#GetMolFrags yes

export Fragmenter

include("dictionaries.jl")
include("algorithm.jl")


end
