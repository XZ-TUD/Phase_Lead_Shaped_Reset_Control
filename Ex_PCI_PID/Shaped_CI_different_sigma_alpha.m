clc;
clear;
% close all;
s = tf('s');
Nh = 2;
freqs = logspace(0,4,1e3);
Ws = 2*pi*freqs;
gamma = 0;
wc = 2*pi*80;
wi = wc*0.1;
wf = wc*10;
sys = 6.615e5/(83.57*s^2+279.4*s+5.837e5);
phase_sys = rad2deg(angle(freqresp(sys, wc)));
lpf = 1/(s/wf+1);
phase_lpf = rad2deg(angle(freqresp(lpf, wc)));



PI = wi/s+1;
PI_freqs(1,:) = freqresp(PI, Ws);
phase_PI = rad2deg(angle(freqresp(PI, wc)));

PCI =1/s; 
[A, B, C, D] = ssdata(PCI);
k_pci = (abs(freqresp(PI, wc)))/(abs(hosidfcalc(ss(A, B, C, D), gamma, 1, wc)));

PCI = PCI*0.8787*1.05; 
[A, B, C, D] = ssdata(PCI);
PCI_freqs = hosidfcalc(ss(A,B,C,D), gamma, 1, Ws);

% Design Der
ReqPM = 10;
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
% wl = 950;
% wh = 1e4;
% Cs = (s/wl+1)/(s/wh+1)/(s/3000+1)^2;
%%

% sigma = 0
% wl1 = 1;
% wh1 = 1;
% wl2 = 1;

% sigma = 0.01, 8.1096
% wl1 = 526;
% wh1 = 700;
% wl2 = 1e4;

% sigma = 0.05, 18.1949
% wl1 = 480;
% wh1 = 915;
% wl2 = 1e4;


% sigma = 0.1; 25.8419
% wl1 = 415;
% wh1 = 1050;
% wl2 = 1e4;


% sigma = 0.2; 36.8699
wl1 = 300;
wh1 = 1200;
wl2 = 1e4;


Cs = (s/wl1+1)/(s/wh1+1);%*1/(s/wl2+1);
Cs_freqs(1,:) = freqresp(Cs, Ws);
% figure;
% semilogx(freqs, rad2deg(angle(Cs_freqs(1,:))),'-.','LineWidth',2,'Color','#d95319'); hold on;
% grid on;

% Cs = ;
for hn = 1:Nh
    n = 2*hn-1;
    PCI_PID_sys_freqs(hn,:) = hosidfcalc(ss(A,B,C,D), gamma, n, Ws);
    Shaped_PCI_PID_sys_freqs(hn,:) = func_calcr_Cs(1/s, Cs, gamma, n, Ws);
%     Cs_Shaped_PCI_PID_sys_freqs(hn,:) = Cs_freqs.*func_calcr_Cs((1.3*wi/s+1)*0.76, Cs, gamma, n, Ws) .* sys_freqs.* lpf_freqs.*PI_freqs* kp .* Der_freqs;
    
    %     func_calol_Cs(s/s, Cs, PCI, 0*s/s, kp*Der*PI*lpf, s/s, sys, gamma, n, Ws);
end
%% Plot PCID + sys
aa = 1;
bb = 945;
y1 = rad2deg(angle(PCI_PID_sys_freqs(1,:)));
y2 = rad2deg(angle(Shaped_PCI_PID_sys_freqs(1,:)));
% y3 = 0*ones(1,numel(freqs));
% semilogx(freqs, y1,'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, y2,'-.','LineWidth',2,'Color','#d95319'); hold on;
% fill([freqs(1,aa:bb) fliplr(freqs(1,aa:bb))], [y1(1,aa:bb) fliplr(y2(1,aa:bb))], 'g', 'FaceAlpha', 0.05); hold on; % 'g' for green, FaceAlpha controls transparency

%%
% h = figure; 
subplot(3,1,1);
semilogx(freqs, (y2),'LineWidth',2); hold on;
grid on;
ylabel('$\angle\mathcal{C}_1(\omega)$ [deg] ','Interpreter','latex');
set(gca,'fontsize', 12); 
set(gca,'FontName','Times New Roman','fontSize', 12);

subplot(3,1,2);
semilogx(freqs, mag2db(abs(Shaped_PCI_PID_sys_freqs(1,:))),'LineWidth',2); hold on;
grid on;
ylabel('$|\mathcal{C}_1(\omega)|$  [dB]','Interpreter','latex');
set(gca,'fontsize', 12); 
set(gca,'FontName','Times New Roman','fontSize', 12);
lgd = legend('CI','$\sigma_{\alpha}=0.01$','$\sigma_{\alpha}=0.05$','$\sigma_{\alpha}=0.1$','$\sigma_{\alpha}=0.2$','Interpreter','latex');
set(lgd,'fontsize', 10);

% subplot(3,1,3);
% zeta_n_pci_pid(1,:) = abs(PCI_PID_sys_freqs(2,:))./abs(PCI_PID_sys_freqs(1,:));
% zeta_n_Shaped_PCI_PID(1,:) = abs(Shaped_PCI_PID_sys_freqs(2,:))./abs(Shaped_PCI_PID_sys_freqs(1,:));
% semilogx(freqs, (zeta_n_Shaped_PCI_PID(1,:)-zeta_n_pci_pid(1,:)),'LineWidth',2); hold on;
% % semilogx(freqs, (log10(zeta_n_Shaped_PCI_PID(1,:))./log10(zeta_n_pci_pid(1,:))),'LineWidth',2); hold on;
% grid on;
% xlabel('Frequency [Hz]');
% ylabel('$\zeta_3(\omega) - \zeta_3^{0}(\omega)$ [abs]','Interpreter','latex');
% % ylabel('$log(\zeta_3(\omega))/log(\zeta_3^{0}(\omega))$ [dB]','Interpreter','latex');
% set(gca,'fontsize', 12); 
% set(gca,'FontName','Times New Roman','fontSize', 12);
% % lgd = legend('$\sigma_{\alpha}=0.01$','$\sigma_{\alpha}=0.05$','$\sigma_{\alpha}=0.1$','$\sigma_{\alpha}=0.2$','Interpreter','latex');
% % set(lgd,'fontsize', 10);
% h.Position = [100 100 800 500];

subplot(3,1,3);
semilogx(freqs, mag2db(Shaped_PCI_PID_sys_freqs(2,:)),'LineWidth',2); hold on;
grid on;
xlabel('Frequency [Hz]');
ylabel('$|\mathcal{C}_3(\omega)|$ [dB] ','Interpreter','latex');
set(gca,'fontsize', 12); 
set(gca,'FontName','Times New Roman','fontSize', 12);
lgd = legend('CI','$\sigma_{\alpha}=0.01$','$\sigma_{\alpha}=0.05$','$\sigma_{\alpha}=0.1$','$\sigma_{\alpha}=0.2$','Interpreter','latex');
set(lgd,'fontsize', 10);

h.Position = [100 100 800 700];
%%
% % h = figure; 
% subplot(3,1,1);
% % semilogx(freqs, y1,'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, (y2-y1),'LineWidth',2); hold on;
% grid on;
% ylabel('$\angle \mathcal{C}_1(\omega) - \angle \mathcal{C}_1^{0}(\omega) $ [deg]','Interpreter','latex');
% set(gca,'fontsize', 12); 
% set(gca,'FontName','Times New Roman','fontSize', 12);
% 
% subplot(3,1,2);
% % semilogx(freqs, ,'LineWidth',2,'Color','#0072bd'); hold on;
% semilogx(freqs, (abs(Shaped_PCI_PID_sys_freqs(1,:)) - abs(PCI_PID_sys_freqs(1,:)))./abs(PCI_PID_sys_freqs(1,:)),'LineWidth',2); hold on;
% grid on;
% ylabel('$|\mathcal{C}_1(\omega)| - |\mathcal{C}_1^{0}(\omega)|$ [abs]','Interpreter','latex');
% set(gca,'fontsize', 12); 
% set(gca,'FontName','Times New Roman','fontSize', 12);
% 
% subplot(3,1,3);
% zeta_n_pci_pid(1,:) = abs(PCI_PID_sys_freqs(2,:))./abs(PCI_PID_sys_freqs(1,:));
% zeta_n_Shaped_PCI_PID(1,:) = abs(Shaped_PCI_PID_sys_freqs(2,:))./abs(Shaped_PCI_PID_sys_freqs(1,:));
% semilogx(freqs, (zeta_n_Shaped_PCI_PID(1,:)-zeta_n_pci_pid(1,:))./zeta_n_pci_pid(1,:),'LineWidth',2); hold on;
% grid on;
% xlabel('Frequency [Hz]');
% ylabel('$\zeta_3(\omega) - \zeta_3^{0}(\omega)$ [abs]','Interpreter','latex');
% set(gca,'fontsize', 12); 
% set(gca,'FontName','Times New Roman','fontSize', 12);
% lgd = legend('$\sigma_{\alpha}=0.01$','$\sigma_{\alpha}=0.05$','$\sigma_{\alpha}=0.1$','$\sigma_{\alpha}=0.2$','Interpreter','latex');
% set(lgd,'fontsize', 10);
% h.Position = [100 100 800 500];
% 
% figure;
% cos_angle_Cs(1,:) = cos(angle(Cs_freqs(1,:)));
% semilogx(freqs, cos_angle_Cs(1,:));