% House Cleaning
clc
close all

% Import Data
data = importdata("config1.mat");
data2 = importdata("config2-1.mat");
data3 = importdata("config3.mat");

% Constants
R = 287; %[J/kg*K]
P_atm = 83700; %[Pa]
T_atm = 21.8 + 273; %[K]

% Convert Tables To Arrays
both_data = table2array(data);
total_data = table2array(data2);
static_data = table2array(data3);

% Extract Data Needed
delta_P = both_data(:,3);
delta_P_max = max(delta_P);

% Equation for Airspeed
V_air_max = sqrt(2 * delta_P_max * (R*T_atm)/P_atm) %[m/s]

% Error Values
P_atm_error = 50; %[Pa]
T_atm_error = 0.05; %[K]
delta_P_error = 49.76; %[Pa]

% Error Equations
delta_P_partial = (R * T_atm)/(P_atm * sqrt(2 * delta_P_max * R*T_atm/P_atm));
T_atm_partial = (delta_P_max * R)/(P_atm * sqrt(2*delta_P_max*R*T_atm/P_atm));
P_atm_partial = (-delta_P_max * R * T_atm)/((P_atm)^2 * sqrt(2*delta_P_max*R*T_atm/P_atm));

% Error Calculation
V_air_error = sqrt((P_atm_error*P_atm_partial)^2 + (T_atm_error*T_atm_partial)^2 + (delta_P_partial*delta_P_error)^2)