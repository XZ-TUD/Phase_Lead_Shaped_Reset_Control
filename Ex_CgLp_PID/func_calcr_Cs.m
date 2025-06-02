% Developer information:

%%
function [Hr,C_rho_n] = func_calcr_Cs(Cr, Cs, gamma, n, W)
s = tf('s');
[A_R, B_R, C_R, D_R] = ssdata(Cr);
A_rho = eye(size(A_R));
A_rho(1,1) = gamma;
Hr = zeros(1,numel(W));

 for i = 1:numel(W)
     w = W(1,i);     

     Ang_cs = angle(freqresp(Cs*s/s,w));
     Delta = eye(size(A_R)) + exp(pi/w*A_R);
     Delta_r = eye(size(A_R)) + A_rho*exp(pi/w*A_R);
     Omega = Delta - Delta*inv(Delta_r)*A_rho*Delta;
     Lamda = w^2*eye(size(A_R)) + A_R^2;
     Theta_phi = -2*1j*w*exp(1j*Ang_cs)/pi*Omega*(w*eye(size(A_R))*cos(Ang_cs)-A_R*sin(Ang_cs))*inv(Lamda)*B_R;
     if (n == 1)
     Hr(i) = C_R * inv(A_R - 1j*w*eye(size(A_R))) * Theta_phi + C_R * inv(1j*w*eye(size(A_R))-A_R) * B_R + D_R;
     else
     Hr(i) = C_R * inv(A_R-1j*n*w*eye(size(A_R))) * Theta_phi;  
     end
  end
end