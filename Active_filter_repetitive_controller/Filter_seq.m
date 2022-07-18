function [ k,Peaking_s,Peaking_z,Peaking_sc,Peaking_zc,Notch_s,Notch_z,Peaking_inv1_s,Peaking_inv1_d ] = Filter_seq(f1,fs,eta_ref,ft)


h_max=round(fs/f1/2)-1;


eta=eta_ref*ones(1,h_max);


w1=2*pi*f1;
Ts=1/fs;


for h=1:2:h_max
    w(h)=h*w1;
end

s=tf('s');
w_f=2*pi*ft;
F_f=1/(s^2/w_f^2+2/w_f*s+1); %filtro sallen key
[~,phase_f]=bode(F_f,w); %ritardo filtro
[~,phase_i]=bode(exp(-s*2*Ts),w); %ritardo inverter

Tr_Sallen=zeros(h_max,1);
Tr_inv=zeros(h_max,1);
for h=1:2:h_max
    Tr_Sallen(h)=-phase_f(h)/w(h)*pi/180;
    Tr_inv(h)=-phase_i(h)/w(h)*pi/180;
end
Td=Tr_inv+Tr_Sallen;

% Peaking_s=zeros(h_max,1);
% Peaking_z=zeros(h_max,1);
% Peaking_sc=zeros(h_max,1);
% Peaking_zc=zeros(h_max,1);
% Notch_s=zeros(h_max,1);
% Notch_z=zeros(h_max,1);

for h=1:2:h_max
    %filtri di peaking non compensati
    Peaking_s(h)=2*eta(h)*w(h)*s/(s^2+2*eta(h)*w(h)*s+w(h)^2);
    Peaking_z(h)=c2d(Peaking_s(h),Ts,'prewarp',w(h));
    %filtri di peaking compensati
    Peaking_sc(h)=2*eta(h)*(w(h)*s*cos(w(h)*Td(h))-w(h)^2*sin(w(h)*Td(h)))/(s^2+2*eta(h)*w(h)*s+w(h)^2);
    Peaking_zc(h)=c2d(Peaking_sc(h),Ts,'prewarp',w(h));
    %filtri di notch 
    Notch_s(h)=1-Peaking_s(h);
    Notch_z(h)=1-Peaking_z(h);
end


end
 
