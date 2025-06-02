clc;
clear;
% close all;
% 7.8 degree phase lead
s = tf('s');
freqs = logspace(0,3,1e3);
Ws = 2*pi*freqs;
gamma = -0.3;
Nh = 2;
wc = 50*2*pi;
wi = wc*0.1;
wf = wc*10;
sys = 6.615e5/(83.57*s^2+279.4*s+5.837e5);
phase_sys = rad2deg(angle(freqresp(sys, wc)));
lpf = 1/(s/wf+1);
phase_lpf = rad2deg(angle(freqresp(lpf, wc)));


%% Design PID
PI = wi/s+1;
PI_freqs(1,:) = freqresp(PI, Ws);
phase_PI = rad2deg(angle(freqresp(PI, wc)));

% Design Der
ReqPM = 50;
ReqPh = (-180 + ReqPM - phase_sys - phase_lpf - phase_PI)*2*pi/360; 
dscale = tan((ReqPh + pi/2)/2); 
wd = wc/dscale;
wt = wc*dscale;
Der = (s/wd+1)/(s/wt+1);
Der_freqs(1,:) = freqresp(Der, Ws);

% pid
PID_freqs(1,:) = PI_freqs(1,:) .* Der_freqs(1,:);

% kp
sys_freqs(1,:) = freqresp(s/s, Ws);
lpf_freqs(1,:) = freqresp(lpf, Ws);
G = abs(freqresp(sys*Der*PI, wc));
kp = 1/G;
PI_freqs(1,:) = PI_freqs(1,:) * kp;
PID_freqs(1,:) = PID_freqs(1,:) * kp .* lpf_freqs;
PID_sys_freqs(1,:) = PID_freqs(1,:) .* sys_freqs;

%% Design CgLp

PhCglpAdd = 20;
[scale,offset] = CgLp_design(gamma,1,PhCglpAdd,1.0,1.0);
wr = wc/scale;
wr_alpha = wr/offset;
FORE = 1/(s/wr_alpha+1);
[A, B, C, D] = ssdata(FORE);
FORE_freqs = hosidfcalc(ss(A,B,C,D), gamma, 1, Ws);

Lead = (s/wr+1)/(s/(100*wc)+1);
Lead_freqs(1,:) = freqresp(Lead, Ws);

% Der2
ReqPh = (-180 + ReqPM - phase_sys - phase_lpf - phase_PI - PhCglpAdd)*2*pi/360; 
dscale2 = tan((ReqPh + pi/2)/2); 
wd2 = wc/dscale2;
wt2 = wc*dscale2;
Der2 = (s/wd2+1)/(s/wt2+1);
G = abs(freqresp(sys*PI*Der2*Lead*lpf, wc)*hosidfcalc(ss(A,B,C,D),gamma,1,wc));
kp2 = (1/G);
PID2 = (wi/s+1) * (s/wd2+1)/(s/wt2+1) * kp2 * lpf;
PID2_freqs(1,:) = freqresp(PID2, Ws);


% CgLp
CgLp_freqs = FORE_freqs .* Lead_freqs; 

% CgLp+PID
CgLp_PID_freqs = CgLp_freqs .* PID2_freqs; 

% CgLp+PID+SYS
CgLp_PID_sys_freqs = CgLp_freqs .* PID2_freqs(1,:) .* sys_freqs; 

C_alpha_P = Lead  * PID2 * s/s;

%% Design Shaped PCID
% Improved steady-state responses
% wl = 500;
% wh = 5000;
% lpf2 = 1/(s/1500+1)^2;
% Cs = (s/wl+1)/(s/wh+1)*lpf2;

wl = 950;
wh = 1e5;
lpf2 = 1/(s/2000+1);
Cs = (s/wl+1)/(s/wh+1)*lpf2;
% Cs = ((s/wl+1)/(s/wh+1))^2*(1/(s/2000+1))^2;

% wl = 950;
% wh = 1e4;
% lpf2 = 1/(s/2000+1);
% Cs = (s/wl+1)/(s/wh+1)*lpf2;
% wl1 = 1e3;
% wh1 = 1e5;
% Cs = (s/wl1+1)/(s/wh1+1);
Cs_freqs(1,:) = freqresp(Cs, Ws);

lpf2_freqs(1,:) = freqresp(lpf2, Ws);


for hn = 1:Nh
    n = 2*hn-1;
    Wns = n*Ws;
    C_alpha_P_freqs(hn,:) = freqresp(C_alpha_P, Wns);
    CgLp_PID_sys_freqs(hn,:) = hosidfcalc(ss(A,B,C,D), gamma, n, Ws) .* C_alpha_P_freqs(hn,:);
    % Improved steady-state responses
%     Shaped_CgLp_PID_sys_freqs(hn,:) = func_calcr_Cs(1/(s/wr_alpha*1.1+1)*1.4, Cs, gamma+0.28, n, Ws) .* C_alpha_P_freqs(hn,:);
    Shaped_CgLp_PID_sys_freqs(hn,:) = func_calcr_Cs(1/(s/wr_alpha*1.1+1)*1.8, Cs, gamma+0.28, n, Ws) .* C_alpha_P_freqs(hn,:);
    Cs_Shaped_CgLp_PID_sys_freqs(hn,:) =  Cs_freqs.* Shaped_CgLp_PID_sys_freqs(hn,:);
end

%% Plot CgLp
h = figure; 
subplot(3,1,1);
semilogx(freqs, mag2db(abs(PID_sys_freqs)),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs2(1,:))),'LineWidth',2,'Color','b'); hold on;
semilogx(freqs, mag2db(abs(Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;
semilogx(freqs, mag2db(abs(Cs_Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#999999'); hold on;
grid on;
ylabel('Magnitude [dB]');
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
lgd = legend('PID', 'CgLp-PID: 1^s^t Harmonic', 'Shaped CgLp-PID: 1^s^t Harmonic');
set(lgd,'fontsize', 10);

subplot(3,1,2);
semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#0072bd'); hold on;
semilogx(freqs, mag2db(abs(Shaped_CgLp_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#d95319'); hold on;
grid on;
ylabel('Magnitude [dB]');
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
% lgd = legend('PID', 'CgLp-PID', 'Shaped CgLp-PID');
lgd = legend('CgLp-PID: 3^r^d Harmonic', 'Shaped CgLp-PID: 3^r^d Harmonic');
set(lgd,'fontsize', 10);


subplot(3,1,3);
semilogx(freqs, rad2deg(angle(PID_sys_freqs)),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, rad2deg(angle(CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, rad2deg(angle(CgLp_PID_sys_freqs2(1,:))),'LineWidth',2,'Color','b'); hold on;
semilogx(freqs, rad2deg(angle(Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;
xlabel('Frequency [rad/s]');
ylabel('Phase [Degree]');grid on;
h.Position = [100 100 800 700];
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
% lgd = legend('PID', 'CgLp-PID', 'Shaped CgLp-PID');
% set(lgd,'fontsize', 14);


%% Plot CgLp
h = figure; 
subplot(2,1,1);
semilogx(freqs, mag2db(abs(PID_sys_freqs)),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs2(1,:))),'LineWidth',2,'Color','b'); hold on;
semilogx(freqs, mag2db(abs(Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;

semilogx(freqs, mag2db(abs(CgLp_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#0072bd'); hold on;
semilogx(freqs, mag2db(abs(Shaped_CgLp_PID_sys_freqs(2,:))),'--','LineWidth',2,'Color','#d95319'); hold on;
grid on;
ylabel('Magnitude [dB]');
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
% lgd = legend('PID', 'CgLp-PID', 'Shaped CgLp-PID');
lgd = legend('PID','CgLp-PID: 1^s^t Harmonic', 'Shaped CgLp-PID: 1^s^t Harmonic','CgLp-PID: 3^r^d Harmonic', 'Shaped CgLp-PID: 3^r^d Harmonic');
set(lgd,'fontsize', 10);

subplot(2,1,2);
semilogx(freqs, rad2deg(angle(PID_sys_freqs)),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, rad2deg(angle(CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, rad2deg(angle(CgLp_PID_sys_freqs2(1,:))),'LineWidth',2,'Color','b'); hold on;
semilogx(freqs, rad2deg(angle(Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;
xlabel('Frequency [Hz]');
ylabel('Phase [Degree]');grid on;
h.Position = [100 100 800 500];
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);
% lgd = legend('PID', 'CgLp-PID', 'Shaped CgLp-PID');
% set(lgd,'fontsize', 14);

h = figure; 
semilogx(freqs, (abs(CgLp_PID_sys_freqs(2,:)))./(abs(CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#edb120'); hold on;
semilogx(freqs, (abs(Shaped_CgLp_PID_sys_freqs(2,:)))./ (abs(Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#d95319'); hold on;
% semilogx(freqs, mag2db(abs(Cs_Shaped_CgLp_PID_sys_freqs(1,:))),'LineWidth',2,'Color','#999999'); hold on;
grid on;
ylabel('Magnitude [dB]');
set(gca,'fontsize', 16); 
set(gca,'FontName','Times New Roman','fontSize', 16);