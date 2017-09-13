# OptPipe - a Pipeline for Optimizing Metabolic Engineering Targets
---
This software contains up-to date standardized methods for the optimization of metabolic engineering targets in genome-scale metabolic networks. This software was developed in the frame of the
[BaCHBerry] project (European Union FP7- 613793)

The code was written by Ana Vila Santa and András Hartmann under the supervision of Susana Vinga.

# License
---
© 2015 - 2017 [Susana Vinga], [András Hartmann] and [Ana Vila Santa]

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
    - Tested on version 2015a and 2016b
    - Compatibility issues might arise between gurobi and MATLAB 2017a
    - Parallel toolbox (optional, but highly recommended in order to shorten running time)
- [COBRA toolbox] from the openCOBRA project (release 2.0)
- [SBML toolbox]
- [gurobi] optimization software installed for MATLAB (free for academic users)

# Installation
Place the content of the package in a separate folder, start MATLAB.
You can use git to download the current version from the repository
```shell
$ git clone --recursive https://github.com/AndrasHartmann/OptPipe.git
```
Within the MATLAB shell go to the code folder and run the install script
```matlab
>> cd '<code folder>'
>> installPipeline 
```

[BaCHBerry]: <http://www.bachberry.eu/>

[MATLAB]: <http://www.mathworks.com/products/matlab>
[COBRA toolbox]: <http://opencobra.github.io/cobratoolbox/>
[SBML toolbox]: <http://sbml.org/Software/SBMLToolbox>
[gurobi]: <http://www.gurobi.com/>

[Susana Vinga]: <mailto:susanavinga@tecnico.ulisboa.pt>
[András Hartmann]: <mailto:andras.hartmann@gmail.com>
[Ana Vila Santa]: <mailto:vilasanta.ana@gmail.com>
