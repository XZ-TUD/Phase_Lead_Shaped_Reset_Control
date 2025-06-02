%% case studies in the paper
clc
clear
close all
s = tf('s');

%% case study 1
omega_c = 2*pi*80;
gamma = -0.3;
omega_alpha = 0;
Cs = (s/950+1)/(s/3000+1)*1/(s/10^4+1);
[num_Cs,den_Cs] = tfdata(Cs);
phi_lead_case1 = func_phi_lead_cal(num_Cs,den_Cs, gamma, omega_alpha, omega_c)

%% case study 2
omega_c = 2*pi*50;
gamma = -0.3;
omega_alpha = 1.602019602515524e+02;
Cs = (s/950+1)/(s/2000+1)*1/(s/10^5+1);
[num_Cs,den_Cs] = tfdata(Cs);
phi_lead_case2 = func_phi_lead_cal(num_Cs,den_Cs, gamma, omega_alpha, omega_c)
