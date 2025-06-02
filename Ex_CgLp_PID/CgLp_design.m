% CgLp_design code
function [scale,offset] = CgLp_design(gamma,order,phase_req,beta,zeta)

% fc = 100 Hz
wc = 2*pi*100;

vals = fminsearch(@(vals) DF_error(order, wc, gamma, beta, zeta, vals, phase_req),[1.0 1.0]);

offset = vals(1);
scale = vals(2);

s = tf('s');
wr = wc/scale;

if(order == 1)
    reset = 1/(s/wr*offset + 1);
    non_reset = (s/wr + 1)/(s/wc/1000 + 1);
    [A1,B1,C1,D1] = ssdata(reset);
    [A2,B2,C2,D2] = ssdata(non_reset);
    [A, B, C, D] = ssseries(A1, B1, C1, D1, A2, B2, C2, D2);
    Arho = eye(size(A,1));
    Arho(1,1) = gamma;
else
    reset = 1/((s/wr*offset)^2 + 2*s*beta/wr*offset + 1);
    non_reset = ((s/wr)^2 + 2*s*zeta/wr + 1)/(s/wc/1000 + 1)/(s/wc/1000 + 1);
    [A1,B1,C1,D1] = ssdata(reset);
    [A2,B2,C2,D2] = ssdata(non_reset);
    [A, B, C, D] = ssseries(A1, B1, C1, D1, A2, B2, C2, D2);
    Arho = eye(size(A,1));
    Arho(1,1) = gamma;
    Arho(2,2) = gamma;
end

W = logspace(0,3,1e3);
CgLpval = hosidfcalc(ss(A,B,C,D), Arho, 1, W);
% copy niranjan "open-loop"这篇文章里的hosidfcalc
% figure;
% subplot(211);semilogx(W/2/pi,20*log10(abs(CgLpval)),'LineWidth',2);hold on;grid on;
% subplot(212);semilogx(W/2/pi,180/pi*(angle(CgLpval)),'LineWidth',2);hold on;grid on;
end