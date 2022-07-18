function [Notch_z_AM1,Notch_z_AM2] = Notch(f1,fc,f_am,eta)


s=tf('s');
w_c=2*pi*f_am;
w_m=2*pi*f1;
w1=w_c+w_m;
w2=w_c-w_m;
Peaking_s_AM1=2*eta*w1*s/(s^2+2*eta*w1*s+w1^2);
Peaking_z_AM1=c2d(Peaking_s_AM1,1/fc,'prewarp',w1);
Notch_z_AM1=1-Peaking_z_AM1;
Peaking_s_AM2=2*eta*w2*s/(s^2+2*eta*w2*s+w2^2);
Peaking_z_AM2=c2d(Peaking_s_AM2,1/fc,'prewarp',w2);
Notch_z_AM2=1-Peaking_z_AM2;