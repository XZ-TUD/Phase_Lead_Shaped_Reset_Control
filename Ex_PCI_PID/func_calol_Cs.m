% Developer information:
% Xinxin Zhang
% x.zhang-15@tudelft.nl
% 2024.05.05
% TU Delft, NL

% Instruction:
% This function is to calculate the HOSIDF of an open-loop reset system, defined as Hn 
% There are 10 inputs in this function.
% Cr: the transfer function of the reset controller.
% Ang_Cs [rad]: angle of the shaping filter Cs.
% gamma: reset value.
% n: n-th harmonic, n=2k+1, k \in N.
% W: Frequency [rad/s].
% C1, C2, C3, C4, and P are transfer functions of elements in the reset
% system.

%% 
function [Hn] = func_calol_Cs(C1, Cs, Cr, C2, C3, C4, P, gamma, n, W)
s = tf('s');
Gr = func_calcr_Cs(Cr, Cs, gamma, n, W);
 for i = 1:numel(W)
     w = W(1,i);
     G1(1,i) = freqresp(C1, w);
     G2(n,i) = freqresp(C2*s/s, n*w);
     G3(n,i) = freqresp(C3*s/s, n*w);
     Gp(n,i) = freqresp(P*s/s, n*w);
     G4(n,i) = freqresp(C4*s/s, n*w);
     if (n == 1)
        Hn = G1(1,i)*(Gr + G2(n,i))*G3(n,i)*Gp(n,i)*G4(n,i);
     else
        Hn = G1(1,i)*Gr*G3(n,i)*Gp(n,i)*G4(n,i)*exp(1j*(n-1)*angle(G1(i)));
        
     end
 end
end