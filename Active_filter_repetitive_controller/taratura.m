clc
clear all

s=tf('s');
%%anello di corrente
R_disaccoppiamento=2*13e-3;
L_disaccoppiamento=2*1.18e-3;
Tc=1/10000;
w0=2*pi*50;
G_dis=1/(R_disaccoppiamento+s*L_disaccoppiamento);
G_inv=exp(-1.5*s*Tc);
G_inv_a=1/(1+s*3/2*Tc);
Gv=G_inv*G_dis;

ki_pi=45;
tau=1/11.109;
kp_pi=tau*ki_pi;

G_pi_C=kp_pi+ki_pi/s;

G_dir_corr=G_pi_C*Gv;

G_anello_corr=G_dir_corr/(1+G_dir_corr);

%anello tensione
Vdc_ref=200;
ampiezza_rete=150/sqrt(3);
c=2.2e-3;
G_dc=3/2/s/c*ampiezza_rete/Vdc_ref;

G_compensare=G_dc*G_anello_corr;

Kp_PIDC=1.5723;


G_p_DC=Kp_PIDC;
G_completa_catena_dir=G_anello_corr*G_p_DC*G_dc;
G_completa=G_completa_catena_dir/(1+G_completa_catena_dir);
filtro_vdc=1/(1+0.15*s);



