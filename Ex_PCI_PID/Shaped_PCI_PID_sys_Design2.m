clc;
clear;
close all;
s = tf('s');
Nh = 2;
freqs = logspace(0,3,1e3);
Ws = 2*pi*freqs;
gamma = -0.3;
wc = 2*pi*70;
wi = wc*0.1;
wf = wc*10;
sys = 6.615e5/(83.57*s^2+279.4*s+5.837e5);
phase_sys = rad2deg(angle(freqresp(sys, wc)));
lpf = 1/(s/wf+1);
phase_lpf = rad2deg(angle(freqresp(lpf, wc)));



PI = wi/s+1;
PI_freqs(1,:) = freqresp(PI, Ws);
phase_PI = rad2deg(angle(freqresp(PI, wc)));

PCI = 1.3095*wi/s+1; 
[A, B, C, D] = ssdata(PCI);
k_pci = (abs(freqresp(PI, wc)))/(abs(hosidfcalc(ss(A, B, C, D), gamma, 1, wc)));

PCI = PCI*k_pci; 
[A, B, C, D] = ssdata(PCI);
PCI_freqs = hosidfcalc(ss(A,B,C,D), gamma, 1, Ws);

% Design Der
ReqPM = 20;
phase_pci = rad2deg(angle(hosidfcalc(ss(A,B,C,D), gamma, 1, wc)));
ReqPh = (-180 + ReqPM - phase_pci  - phase_sys - phase_lpf - phase_PI)*2*pi/360; 
dscale = tan((ReqPh + pi/2)/2); 
wd = wc/dscale;
wt = wc*dscale;
Der = (s/wd+1)/(s/wt+1);
Der_freqs(1,:) = freqresp(Der, Ws);

% PCID
PID_freqs(1,:) = PI_freqs(1,:) .* Der_freqs(1,:);
PCID_freqs(1,:) = PCI_freqs(1,:) .* Der_freqs(1,:);


%% series PCID + sys + lpf
lpf_freqs(1,:) = freqresp(lpf, Ws);
sys_freqs(1,:) = freqresp(s/s, Ws);
PIID_sys_freqs(1,:) = PID_freqs .* sys_freqs.* lpf_freqs.*PI_freqs;
PCI_PID_sys_freqs(1,:) = PCID_freqs .* sys_freqs.* lpf_freqs.*PI_freqs;
sys_PID_freqs(1,:) = freqresp(PI*sys*Der, wc);

% kp
G = abs(sys_PID_freqs.*hosidfcalc(ss(A,B,C,D),gamma,1,wc));
kp = 1/G;
PIID_sys_freqs(1,:) = PIID_sys_freqs(1,:) * kp;
PCI_PID_sys_freqs(1,:) = PCI_PID_sys_freqs(1,:) * kp;

%% Design Shaped PCID

wl = 1e3;
wh = 3e4;
Cs = (s/wl+1)/(s/wh+1);

% Cs = ;
for hn = 1:Nh
    n = 2*hn-1;
    PCI_PID_sys_freqs(hn,:) = hosidfcalc(ss(A,B,C,D), gamma, n, Ws) .* sys_freqs.* lpf_freqs.*PI_freqs* kp .* Der_freqs;
    Shaped_PCI_PID_sys_freqs(hn,:) = func_calcr_Cs((1.3*wi/s+1)*0.76, Cs, gamma, n, Ws) .* sys_freqs.* lpf_freqs.*PI_freqs* kp .* Der_freqs;
%     func_calol_Cs(s/s, Cs, PCI, 0*s/s, kp*Der*PI*lpf, s/s, sys, gamma, n, Ws);
end
%% Plot PCID + sys
h = figure; 
subplot(2,1,1);
semilogx(freqs, mag2db(abs(PIID_sys_freqs)),'LineWidth',2,'Color','#0072bd'); hold on;
semilogx(freqs, mag2db(abs(PCI_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, mag2db(abs(Shaped_PCI_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;


semilogx(freqs, mag2db(abs(PCI_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, mag2db(abs(Shaped_PCI_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#d95319'); hold on;
grid on;
ylabel('Magnitude [dB]');
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
lgd = legend('PI^2D', 'PCI-PID: 1^s^t Harmonic', 'Shaped PCI-PID: 1^s^t Harmonic', 'PCI-PID: 3^r^d Harmonic', 'Shaped PCI-PID: 3^r^d Harmonic');
set(lgd,'fontsize', 10);

subplot(2,1,2);
semilogx(freqs, rad2deg(angle(PIID_sys_freqs)),'LineWidth',2,'Color','#0072bd'); hold on;
semilogx(freqs, rad2deg(angle(PCI_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, rad2deg(angle(Shaped_PCI_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;
xlabel('Frequency [Hz]');
ylabel('Phase [Degree]');grid on;
h.Position = [100 100 800 500];
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
% lgd = legend('PIID', 'PCI-PID');
% set(lgd,'fontsize', 14);
