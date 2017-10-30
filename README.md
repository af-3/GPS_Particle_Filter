# GPS_Particle_Filter
SISR particle filter used for tracking a runner in Forest Park, Portland, Oregon

File overview:

basemapImage.py | Creates a visual of the path traveled.

dataPrep.py | The gps data is converted from the gpsx file format to a compatible matlab format.

particleFilter.m | Used for Monte Carlo Simulation for a Seqential Importance Sampling with Resampling (SISR) Particle Filter using synthetic data. Outputs the MSE error between the estimated and true states.

particleFilterRealData.m | Script used for a Seqential Importance Sampling with Resampling (SISR) Particle Filter using data from personal runs. Generates estimates for each state element and plots observed vs estimated path, particle trajectories, and particle weights.

Summary & Report | Introduction of problem, state space model, and results. Report provides more detail.
