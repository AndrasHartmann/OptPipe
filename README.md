# OptPipe - a Pipeline for Optimizing Metabolic Engineering Targets
---
This software contains up-to date standardized methods for the optimization of metabolic engineering targets in genome-scale metabolic networks. This software was developed in the frame of the
[BaCHBerry] project (European Union FP7- 613793)

The code was written by Ana Vila Santa and András Hartmann under the supervision of Susana Vinga.

# License
---
© 2015 - 2016 [Susana Vinga], [András Hartmann] and [Ana Vila Santa]

OptPipe is a free software: you can redistribute it and/or modify
it under the terms of the GNU Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Public License for more details.

You should have received a copy of the GNU Public License
along with this code. If not, see http://www.gnu.org/licenses/

# Contents
---

| Folder             |          Description                                |
|--------------------|-----------------------------------------------------|
|common_functions    | Folder for the common function shared by the methods|
|examples            | Folder for working examples including models|
|external            | Folder containing external libraries and toolboxes|
|methods             | Folder containing the methods|
|installPipeline.m   | Install script|
|LICENSE.txt         | License file (GPL v3)|
|README.md           | This file|

# Install
---
## Prerequisites
You will need:
- [MATLAB] software for scientific computing 
    - Tested on version 7.14
    - There are known compatibility issues with 2013b or later
    - Parallel toolbox (optional, but highly recommended in order to shorten running time)
- [COBRA toolbox] from the openCOBRA project
- [SBML toolbox]
- [gurobi] optimization software installed for MATLAB (free for academic users)

# Installation
Place the content of the package in a separate folder, start MATLAB.
Within the MATLAB shell go to the init folder and run the init script
```matlab
>> cd '<code folder>/init'
>> initPipeline
```

[BaCHBerry]: <http://www.bachberry.eu/>

[MATLAB]: <http://www.mathworks.com/products/matlab>
[COBRA toolbox]: <http://opencobra.github.io/cobratoolbox/>
[SBML toolbox]: <http://sbml.org/Software/SBMLToolbox>
[gurobi]: <http://www.gurobi.com/>

[Susana Vinga]: <mailto:susanavinga@tecnico.ulisboa.pt>
[András Hartmann]: <mailto:andras.hartmann@gmail.com>
[Ana Vila Santa]: <mailto:vilasanta.ana@gmail.com>
