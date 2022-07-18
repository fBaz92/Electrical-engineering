function [Peaking_z_wm] = Peaking(f1,fc,eta)


s=tf('s');
w_m=2*pi*f1;
Peaking_s_wm=2*eta*w_m*s/(s^2+2*eta*w_m*s+w_m^2);
Peaking_z_wm=c2d(Peaking_s_wm,1/fc,'prewarp',w_m);
end