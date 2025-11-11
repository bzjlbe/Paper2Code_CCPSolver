# Paper2Code_CCPSolver
MATLAB code for "A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics"

Research Background
The nonsmooth modeling method based on the cone complementarity problem (CCP) is one of the most effective approaches for solving spatial frictional contact problems. The conventional method incorporates a relaxation term and a prediction term related to the tangential relative velocity within the normal complementarity condition, making it more manageable to solve and extendable to systems with high tangential relative velocity.

Problem Addressed
When the analyzed system experiences significant changes in tangential relative velocity or has a high friction coefficient, the discrepancy between the prediction term and the relaxation term can lead to unreasonable numerical oscillations in the model, particularly at the initial simulation values. This instability adversely affects the reliability of the simulation.

Proposed Method
This paper proposes an improved numerical method—the prediction-correction method. Building upon the existing method, we introduce an iterative correction mechanism for the predicted value of the tangential relative velocity. This mechanism eliminates artificial errors and the resulting numerical oscillation issues, thereby providing an accurate solution to the spatial frictional contact problem, particularly in cases involving high friction coefficients.

Repository Structure
The code is organized into three main folders:

figs/: Scripts to generate all figures presented in the paper

Analytical/: Analytical solutions for the three numerical examples

funcs/: Core functions for solving the dynamic equations

Core Functions
PM_nonsm_.m: Main function for the three numerical examples

config.m: Generalized-α solver parameters configuration

contact_detection.m: Contact detection function

dynamics_modeling.m: Dynamic modeling script

dynamics_solver.m: Dynamic solution script

cone_complementarity_solver.m: Core cone complementarity iterative solver with prediction-correction mechanism

APGD.m: APGD solution algorithm

Setup
Before running the code, add all folder paths to MATLAB.
Usage To start the program, either create your own main.m or run the files in the figs folder.

Validation
Several numerical examples are provided to validate the effectiveness of the proposed prediction-correction method in handling systems with high friction coefficients and significant tangential relative velocity changes.
