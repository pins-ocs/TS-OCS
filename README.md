# TS-OCS
 Test Suite for Optimal Control Software
 
## Intro
The idea of this repository generates from the work done when writing the paper *"Notes on Numerical Methods for Solving Optimal Control Problems"* submitted to *IEEJ Journal of Industry Applications Special Issue,  Journal IA 2016/03(E) Special Issue on Motion Control and its Related Technologies* the 1st of June 2015.
At present the available code compares three different solution methods for optimal control problem to one example. The selected example is a non linear optimal control problem, yet simple, can be numerically difficult depending on the choosen values of the problem's parameters. Under some conditions, the problem also has the analytical solution that can be used to evaluate the solution methods accuracy.  
The problem has been formulated and solved according the following solution methods:

 * Differential Dynamic Programming using the code freely available at http://www.idsc.ethz.ch/Downloads/DownloadFiles/dpm
 * Indirect method with finite difference approximation of the obtained BVP. The non linear system is solved using a number of numerical solvers written in Matlab code and PINS http://pins.github.io.
 * Direct Method: direct transcription (same finite difference approximation)  into a Non Linear Programming and solution using various optimiser: IPOPT, Matlab fmincon().

In the future this repository will be extended in order to become a test suite for optimal control software with code publicly available and a large set of benchmark tests.

## Instructions

The code can be downloaded and freely used as is. 

## How to cite


## Contacts
Francesco Biral (francesco.biral@unitn.it)
Enrico Bertolazzi (enrico.bertolazzi@unitn.it)
Paolo Bosetti (paolo.bosetti@unitn.it)
