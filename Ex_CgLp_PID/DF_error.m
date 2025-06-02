function error = DF_error(order, wc, gamma, beta, zeta, vals, phase_req)
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
    W = logspace(log10(wc) - 2,log10(wc) + 3,51);
%     CgLpval = hosidfcalc(ss(A,B,C,D), Arho, 1, wc);
%     CgLpvalh = hosidfcalc(ss(A,B,C,D), Arho, 1, wc*1000);
    CgLpval = hosidfcalc(ss(A,B,C,D), Arho, 1, W);
    CgLpvalwc = hosidfcalc(ss(A,B,C,D), Arho, 1, wc);
%     error = abs(20*log10(abs(CgLpval))) + abs((180/pi*angle(CgLpval)) - phase_req) + 0.5*abs(20*log10(abs(CgLpvalh)));
    error = 0.25*sum(abs(20*log10(abs(CgLpval(1:10))))) + 0.25*sum(abs(20*log10(abs(CgLpval(32:51))))) + sum(abs(20*log10(abs(CgLpval(11:31))))) + 31*abs((180/pi*angle(CgLpvalwc)) - phase_req);
end