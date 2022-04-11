# UNIFACGroupIdentification

A Group Fragmentation Algorithm, julia port of [the python repository](https://github.com/simonmb/fragmentation_algorithm). 

It requires to specify an order, in contrast with the original paper, but in line with the python repository.

## Usage

(TODO, proposed interface)
```
frag = Fragmenter(args...,kwargs...) #default constructor
unifac_frag = AutoFragmenter(:UNIFAC) #Fragmenter with UNIFAC info
gammaMie_frag = AutoFragmenter(:SAFTGammaMie) 
fragment("CCCC",unifac_frag) #["CH3" => 2 "CH2" => 2]
```
This repository provides the same groups as those available in [Clapeyron](https://github.com/ypaul21/Clapeyron.jl).


  