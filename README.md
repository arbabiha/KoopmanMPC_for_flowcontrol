# KoopmanMPC_for_flowcontrol
This project demonstates the application of Koopman-MPC framework for flow control with the example of Burgers equation,
following the paper
"A data-driven Koopman model predictive control framework for nonlinear flows" 
by H. Arbabi, M. Korda and I. Mezic

files in the root folder:

BurgersExample: Runs the Burgers example as explained in the paper, it includes data collection, Extended Dynamic Mode Decomposition for identification of the Koopman linear system and a run of closed-loop controlled system from some initial condition
Feel free to play with theparemeters of the code, in particular, try different observables, embedding dimension, reference signal, initial condition and etc.
The whole program, with the initial paremeter settings, runs on my personal laptop in under 2 minutes.


Before you run the code:

go to subfolder ".\thehood\qpOASES-3.1.0\interfaces\matlab" and run make.m .
This is required to activate the qpOASIS interface for solving the optimization problem


send comments and questions to
arbabiha@gmail.com

H Arbabi

April 2018
