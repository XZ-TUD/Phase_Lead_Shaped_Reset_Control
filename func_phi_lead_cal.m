%% This function refers to Remark 4 in the paper
function [phi_lead] = func_phi_lead_cal(num_Cs,den_Cs, gamma, omega_alpha, omega_c)
s = tf('s');
Cs = tf(num_Cs,den_Cs);
Cs_fres = freqresp(Cs,omega_c);
angle_Cs = angle(Cs_fres);
if omega_alpha ==0
    Lambda = omega_c^2 + omega_alpha^2;
    phi_lambda = atan((cos(2*angle_Cs)+1)/((1+gamma)*pi/(2*(1-gamma))-sin(2*angle_Cs)));
    angle_C1 = phi_lambda-pi/2;
    angle_C1_0 =  atan(4*(1-gamma)/((1+gamma)*pi))-pi/2;
    phi_lead = rad2deg(angle_C1-angle_C1_0);
    % rad2deg(angle_C1_0) uncomment for test
elseif omega_alpha > 0
    phi_s = angle(Cs_fres);
    Theta = exp(-pi*omega_alpha/omega_c);
    Omega = (1-gamma)*(1+Theta)/(1+gamma*Theta);
    Lambda = omega_c^2 + omega_alpha^2;
    k_zeta = omega_c*Omega/Lambda/pi;
    k_gamma = omega_c*cos(2*phi_s)+omega_alpha*sin(2*phi_s)+omega_c;
    phi_alpha = atan(1/((k_gamma*k_zeta)^(-1)-tan(phi_s)));
    angle_C1 = phi_alpha-atan(omega_c/omega_alpha);
    angle_C1_0 =  atan(2*omega_c*k_zeta)-atan(omega_c/omega_alpha);
    phi_lead = rad2deg(angle_C1-angle_C1_0);
else
    print("\omega_\alpha should >0");
end
end