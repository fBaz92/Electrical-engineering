%% Dati filtro attivo

clear all
close all
clc


%% Filtro attivo
% Resistenza carico inverter [Ohm]
R_carico=50;
% Capacità condensatore bus DC [F]
C=2.2e-3;
% Tensione iniziale condensatore bus DC [V]
V_cond_iniz=0;

lunghezza=1:5000:0.1;
%% parametri fisici inverter

%diodo
Resistance_D=0.001;
Forw_volt_d=0.8; %prima 0.8
%igbt
Int_res=0.001;
Snubber_res=1e5;



%% Carico Raddrizzatore 
% Resistenza [Ohm]
R_raddr=60;
% Capacità [F]
C_raddr=600e-6;

%% Rete
% Frequenza rete [Hz]
f=50;
% Pulsazione rete [rad/s]
w=2*pi*f;
% Periodo di rete nominale [s]
T=20e-3;

% Ampiezza fondamentale [V]
ampiezza_rete=150/sqrt(3);
% Pendenza variac [-]
pendenza=1000;

% Resistenza rete [Ohm]
R_rete=100e-3;
% Induttanza rete [H]
L_rete=4e-3;

%% Disaccoppiamento
% Resistenza di disaccoppiamento [Ohm]
R_disaccoppiamento=2*13e-3;
% Induttanza di disaccoppiamento [H]
L_disaccoppiamento=2*1.18e-3;
tau_disaccoppiamento=L_disaccoppiamento/R_disaccoppiamento;

%% Riferimenti
% Tensione di riferimento sul bus DC [V]
Vdc_ref=200;
% Corrente di riferimento asse q (riferimento sincrono con la rete) [A]
iq_ref=0;
% Corrente di riferimento per armoniche [A]
i_armoniche_ref=0;

% Costante di tempo del filtro su Vdc_ref [s]
tau_filtro_Vdc_ref=0.15;

%% Campionamento
% Periodo di campionamento [s]
Tc=1/10000;
% Campioni di ritardo [-]
M=round(T/Tc);
%% peaking e notch per selezionare prima armonica
s=tf('s');
eta=0.5;
Peaking_s=2*eta*w*s/(s^2+2*eta*w*s+w^2);
Peaking_z=c2d(Peaking_s,Tc,'prewarp',w);
Notch_s=1-Peaking_s;
Notch_z=1-Peaking_z;

%% Regolazione asse d (idref)
% Guadagno Vdc_ref [S]
Gain_Vdc=0.1;

%% Regolatore PI antiwindup del PLL
% Guadagno proporzionale [rad/(Vs)]
kp_PLL=0.9;
% Guadagno integrale [rad/(Vs)]
ki_PLL=1000;
% Limiti di pulsazione [rad/s]
min_PLL=-1000;
max_PLL=1000;

% Costante di tempo del filtro (PLL) [s]
tau_filtro_PLL=0.008;

%filtro sincronizzazione inversa
[Peaking_dir1_d,Peaking_dir1_s,Peaking_inv_d,Peaking_inv_s]= semipeaking(0.05,50,Tc);

%% Regolatore PI antiwindup (idref)
% Guadagno proporzionale [S]
%kp_PI=0.2;
kp_PI=1.57;
% Guadagno integrale [S]

ki_PI=0.2;
% Limiti di corrente [A]
min_PI=-8;
max_PI=8;

%iq
kp_PI_iq=2.2511;
ki_PI_iq=120;
peak_rete=Peaking(f,1/Tc,0.1);

%% Regolatore risonante
% Costante proporzionale regolatore di corrente [Ohm]

kp_ris=47.5*0.089*2;
% Costante integrale regolatore di corrente [Ohm]

ki_ris=47.5*2;



%%

%Filtro di Sallen Key 
f_taglio=9e3;
w_filtro=2*pi*f_taglio;



%% Modulazione AM

f_am=500; %frequenza portante
k_par_AM=10;

%coefficienti di lagrange
h=Lagrange_coef(f,1/Tc);
%demodulazione
[Notch_z_AM1_5,Notch_z_AM2_5] = Notch(5*f,1/Tc,2*f_am,0.1);
[Notch_z_AM1_7,Notch_z_AM2_7] = Notch(7*f,1/Tc,2*f_am,0.1);
[Notch_z_AM1_11,Notch_z_AM2_11] = Notch(11*f,1/Tc,2*f_am,0.1);
[Notch_z_AM1,Notch_z_AM2] = Notch(0,1/Tc,f_am,0.1);

s=tf('s');
tau_modulaz=1/2/pi/1000;
filtro_pbAM=1/(1+tau_modulaz*s)^6;
filtro_pbAM_d=c2d(filtro_pbAM,Tc,'prewarp',1/tau_modulaz);

Peaking_z_AM5=Peaking(5*f,1/Tc,0.01);
Peaking_z_AM7=Peaking(7*f,1/Tc,0.01);
Peaking_z_AM11=Peaking(11*f,1/Tc,0.01);
%% ritardo inverter

h_rit=Lagrange_coef_rit(1.5);

%% modifica ZHU
Q=0.05*exp(-2*Tc)+0.9*exp(-s*Tc)+0.05;
G_plant=1/(R_disaccoppiamento+s*L_disaccoppiamento);
G_inv=1*exp(-3/2*s*Tc);
G=G_inv*G_plant/(1+G_inv*G_plant);

% calcolo massimo Kr
d=4; 
for i=1:1:40
    [abs_fpb(i),arg_fpb(i)]=bode(Q,2*pi*i*50);
    [abs_g(i),arg_g(i)]=bode(G,2*pi*i*50);
    arg_g(i)=arg_g(i)*pi/180;
    teta_kserie(i)=arg_g(i)+d*Tc*(2*pi*50*i);
    ks_max(i)=2*cos(teta_kserie(i))/abs_g(i);
end

figure

i=50:50:2000;
plot(i,ks_max);

xlim([250,2000])
title('Andamento guadagno serie')
ylabel('Ks ')
xlabel('frequenza [Hz]')
grid on

% 

%% mattavelli

d_mat=4;

for i=1:1:20
    teta_kserie(i)=arg_g(i)+d_mat*Tc*(2*pi*50*i);
    arg_fpb(i)=arg_fpb(i)*pi/180;
    teta_2(i)=teta_kserie(i)-arg_fpb(i);
    teta_3(i)=2*teta_kserie(i)-arg_fpb(i);
    k_asterisco(i)=1+1/abs_g(i)+cos(arg_fpb(i))*(1+1/abs_fpb(i))+cos(teta_3(i))*(1/abs_fpb(i)-1);    
  
   
    ks_max_mat(i)=2*cos(teta_kserie(i))/abs_g(i)/k_asterisco(i)+2*cos(teta_2(i))/abs_g(i)/k_asterisco(i)/abs_fpb(i);
end

figure
i=50:50:1000;
plot(i,ks_max_mat);

xlim([250,1000])
title('Andamento guadagno serie')
ylabel('Ks ')
xlabel('frequenza [Hz]')
grid on


%% taratura regolatori

% pi corrente
G_inv_a=1/(1+1.5*Tc*s);
G_plant2=G_inv_a*G_plant;
ki_id_corr=47.5;
kp_id_corr=0.089*47.5;

G_PI_id=kp_id_corr+ki_id_corr/s;

w_corrente=1.74e3;

G_anello_corr=1/(1+s/w_corrente);
Gv=1.5*150/sqrt(3)/C/Vdc_ref/s;

G_plant3=G_anello_corr*G_inv_a*Gv;


%% filtro pb fase zero

f_vett=1:1:5000;
w_vett=2*pi*f_vett;
c1=[0.05,0.1,0.15];
c0=[0.9,0.8,0.7];

for alfa=0:1:2

for count=1:1:5000
    
    modulo_filtro(alfa+1,count)=2*c1(alfa+1)*cos(2*pi*f_vett(count)*Tc)+c0(alfa+1);
    modulo_filtro(alfa+1,count)=20*log10(modulo_filtro(alfa+1,count));

end
end

figure
plot(f_vett,modulo_filtro(1,:),f_vett,modulo_filtro(2,:),f_vett,modulo_filtro(3,:));
hold on
title('Modulo del Filtro')
ylabel('[dB]')
xlabel('frequenza [Hz]')
grid on