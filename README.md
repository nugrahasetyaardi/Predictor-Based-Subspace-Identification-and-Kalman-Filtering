# Predictor-Based-Subspace-Identification-and-Kalman-Filtering
This project is part of Estimation and Learning in Aerospace course. 
The task is to identify the lateral dynamics of a quadrotor given PRBS position set point as input and position and acceleration in NED frame as output. 
The method employed is PBSID (Predictor-Based Subspace Identification) with VARX and VARMAX model set. 
Then the identified model is used to design a DT-DT kalman filter to estimate the lateral velocity.
Please note that this code requires PBSID toolbox installed!
