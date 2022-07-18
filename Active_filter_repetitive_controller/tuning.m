clc
close all

T=1/50;
Tc=100e-6;
R_plant=2*13e-3;
L_plant=2*1.18e-3;
k_serie=5;
k_parallelo=10;
k_interno=0.99;
R_rete=10e-3;
L_rete=4e-3;


s=tf('s');
filt_rep=exp(-s*(T/2-2*Tc))/(1+k_interno*exp(-s*T/2))*k_serie-k_parallelo;

% Notch=1-2*0.1*s/(s^2+2*0.1*s+(50*2*pi)^2);

Plant=1/(R_plant+L_plant*s);
imp_rete=(R_rete+s*L_rete)*0;
Inverter=exp(-s*2*Tc);

Gol=filt_rep*Plant*Inverter;
figure
plotbode(Gol,1,1/Tc,1);

Gcl=Gol/(1+Gol-imp_rete*Plant);

figure
plotbode(c2d(Gcl,'zoh'),1,1/Tc,0.1);
