function [G] = hosidfcalc(sys, Ar, n, freqs)
    % G = hosidfcalc(SYS, AR, N, FREQS, CLOL)
    % Calculated the higher order (n) describing function for a reset system.
    %
    % SYS is the reset element described in state space
    % Ar is the amount of reset you want to achieve (typical 0)
    % n is the describing function order
    % FREQS contains the frequencies the describing function is calculated for
    
    % Kars Heinen - TU Delft - 2018
    
    % to do; replace inv() by 'matlab \' for faster results
    
    % odd orders will be skipped 
    if (mod(n,2) == 0) 
        G = 0;
        return;
    end  
    
    A = sys.a; B = sys.b; C = sys.c; D = sys.d;
    
    G = zeros(1,numel(freqs));

    for i=1:numel(freqs)
        w = freqs(i);

        Lambda = w*w*eye(size(A)) + A^2;
        LambdaInv = inv(Lambda);

        Delta = eye(size(A)) + expm(A*pi/w);
        DeltaR = eye(size(A)) + Ar*expm(A*pi/w);

        GammaR = inv(DeltaR)*Ar*Delta*LambdaInv; 
        
        ThetaD = (-2*w*w/pi)*Delta*(GammaR-LambdaInv);
        
        if (n==1)
            G(i) = C*inv(j*w*eye(size(A)) - A)*(eye(size(A)) + j*ThetaD)*B;
        else
            % J1 and J2 dissappear
            G(i) = C*inv(j*w*n*eye(size(A)) - A)*j*ThetaD*B;
        end
    end
    
    if (n == 1)
        G = G + D;
    end
 end
