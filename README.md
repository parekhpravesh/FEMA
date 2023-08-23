# FEMA: Fast and efficient mixed-effects algorithm for large sample whole-brain imaging data 

Pravesh Parekh, Chun Chieh Fan, Oleksandr Frei, Clare E. Palmer, Diana M. Smith, Carolina Makowski, John R. Iversen, Diliana Pecheva, Dominic Holland, Robert Loughnan, Pierre Nedelec, Wesley K. Thompson, Donald J. Hagler Jr., Ole A. Andreassen, Terry L. Jernigan, Thomas E. Nichols, and Anders M. Dale 

**Preprint available [on bioRxiv](https://doi.org/10.1101/2021.10.27.466202)**

*This repository contains the code for performing the analyses and creating figures as reported in the FEMA methods paper*


## Requirements
* [MATLAB](http://mathworks.com)

* To run the code(s), you need a copy of FEMA and associated tools/utility functions. [Download the repository here](https://github.com/cmig-research-group/cmig_tools)

* For running the applications on the ACBD Study data, you will need access to the [Adolescent Brain Cognitive Development<sup>SM</sup> data](https://abcdstudy.org/)

## Table of Contentes
* Simulations
  - Simulation 1: effect of binning and parameter recovery
  - Simulation 2: comparison with `fitlme`
  - Simulation 3: comparison of computational time
  - Simulation 4: examining type I error rate
  - Simulation 5: parameter recovery as a function of number of observations
* Applications using the ABCD Study data
  - Application 1: cortical thickness
    - ROI-level analysis and comparison with `fitlmematrix`
    - vertex-wise analysis
  - Application 2: correlation matrix dervied from resting state functional MRI
