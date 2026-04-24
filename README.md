# File Descriptions

- `figure_1a.f90`  
  Generates the bifurcation diagram of the predator population for the uncoupled system, used to characterize the local dynamics.

- `figure_2a.dat`  
  Adjacency matrix for the globally connected network of 20 patches (G1).

- `figure_2b.dat`  
  Adjacency matrix for network G2.

- `figure_2c.dat`  
  Adjacency matrix for network G3.

- `figure_2d_2e_2f.dat`  
  Contains the parameter values of the coupling strength \( \epsilon_c \) and the predation parameter \( \epsilon_p \).

- `figure_2d_2e_2f_5a_initial_condition.f90`  
  Generates initial conditions sampled from regions of the phase space where predator populations persist over long transient times.

- `figure_2d.f90`  
  Computes the mean transient time in the \( \epsilon_c - \epsilon_p \) parameter plane for network G1.

- `figure_3a.f90`  
  Generates time series of predator populations for network G1.

- `figure_4b.f90`  
  Performs the primary simulations used to compute synchronization errors of the proposed networked system.

- `figure_4c.f90`  
  Computes the Master Stability Function (MSF) and the Maximum Lyapunov Exponent (MLE) to quantify the stability of the synchronized state.

- `figure_5a.f90`  
  Performs global stability analysis of the synchronized state by mapping basins of attraction.

- `figure_5c.py`  
  Visualizes the chaotic saddle using a Poincaré section on the \( C \)-plane for the Resource–Predator (R–P) dynamics.

  **Note** : All the remaining figures are computed with the programs which are already provided. We have used **MATLAB** for plotting the data files generated from the above codes.
