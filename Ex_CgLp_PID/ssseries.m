function [A, B, C, D] = ssseries(A1, B1, C1, D1, A2, B2, C2, D2)
% function [A, B, C, D] = ssseries(A1, B1, C1, D1, A2, B2, C2, D2)
% Outputs: SS matrixes of two LTI systems in series (G1 followed by G2).
% Inputs: SS matrixes of G1 and G2.
% DV 2019
A = [A1, zeros(size(A1,1),size(A2,2)); B2*C1, A2];
B = [B1; B2*D1];
C = [D2*C1, C2];
D = D1*D2;