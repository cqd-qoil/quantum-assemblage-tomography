# Quantum Assemblage Tomography

## Overview

This project implements a quantum assemblage tomography method using conical optimization techniques. The primary function, `tomography_m3`, analyzes quantum measurement data to estimate the quantum assemblage, the reduced quantum state of the trusted party, detection efficiencies and outcome biases of the untrusted party, and their corresponding likelihoods.
It assumes a general loss model for the detection of the untrusted party, as described in the manuscript: "Quantum Assemblage Tomography."

Two other algorithms named `tomography_m1` and `tomography_m2` are simplifications of the main algorithm, with stronger assumptions about the loss in the system. A description of these is also found in the manuscript.

## Requirements

To run this project, you need the following:

- MATLAB (version compatible with YALMIP and MOSEK)
- YALMIP toolbox
- MOSEK solver

## Installation

1. **Install YALMIP**: Download and add YALMIP to your MATLAB path. You can find YALMIP [here](https://yalmip.github.io/download/).
2. **Install MOSEK**: Download and install the MOSEK solver from [MOSEK's official website](https://www.mosek.com/).
3. **Set up MOSEK**: Follow the instructions provided by MOSEK to configure it within MATLAB.

## Usage

The tomography functions assume that Bob and Alice both perform three measurements each, with two (three) outcomes for Bob (Alice). The number of measurement inputs and outputs can be modified within the algorithm itself, but care must be taken to adjust the optimitization constraints accordingly.

To use the tomography functions, you need to provide them with measurement data in the form of a structured array. The function returns several outputs, including the quantum assemblage data and efficiency/bias arrays. Note that to run `tomography_m3`, you must first call `tomography_m2` on the experimental data and then pass the output assemblage to `tomography_m3` alongside measurement data. This is because `tomography_m3` starts from a good initial guess provided by another algorithm and then tries to optimize the bias in the system.

### Example

```matlab
% Load your measurement data
bobData = [...]; % Replace with your actual data

% Call the m3 function
[MTestResult, Assem, T, epsx, gamma, Lag, Lag0] = tomography_m3(bobData, Assem, epsx);
```

### License

This project is licensed under the MIT License. See the LICENSE file for more details.

### Acknowledgments

This project is linked to the paper "Quantum Assemblage Tomography."
Special thanks to the developers of YALMIP and MOSEK.

### Contact

For questions or feedback, please contact:

    Yuanlong Wang
    Email: wangyuanlong[at]amss.ac.cn
