function[Peaking_dir1_d,Peaking_dir1_s,Peaking_inv1_d,Peaking_inv1_s]= semipeaking(eta_ref,f1,Tc)

%% semipeaking diretta
w1=2*pi*f1;
fs=1/Tc;
Peaking_dir1_d_s=tf([eta_ref eta_ref^2],[1 2*eta_ref eta_ref^2+(w1)^2]);
Peaking_dir1_q_s=tf([w1*eta_ref],[1 2*eta_ref eta_ref^2+(w1)^2]);

Peaking_dir1_s=[Peaking_dir1_d_s,-Peaking_dir1_q_s;Peaking_dir1_q_s,Peaking_dir1_d_s];
Peaking_dir1_d=c2d(Peaking_dir1_s,1/fs,'tustin');

%% semipeaking inversa


Peaking_inv1_d_s=tf([eta_ref eta_ref^2],[1 2*eta_ref eta_ref^2+(w1)^2]);
Peaking_inv1_q_s=tf([-w1*eta_ref],[1 2*eta_ref eta_ref^2+(w1)^2]);

Peaking_inv1_s=[Peaking_inv1_d_s,-Peaking_inv1_q_s;Peaking_inv1_q_s,Peaking_inv1_d_s];
Peaking_inv1_d=c2d(Peaking_inv1_s,1/fs,'tustin');
