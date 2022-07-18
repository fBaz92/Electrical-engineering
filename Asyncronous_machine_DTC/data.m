%% Induction machine parameters
Pn  = 164e3;    % W,  nominal power
Vn  = 550;      % V,  rms phase-to-phase, rated voltage
fn  = 60;       % Hz, rated frequency

Rs  = 0.0139;   % pu, stator resistance
Lls = 0.0672;   % pu, stator leakage inductance
Rr  = 0.0112;   % pu, rotor resistance, referred to the stator side
Llr = 0.0672;   % pu, rotor leakage inductance, referred to the stator side
Lm  = 2.717;    % pu, magnetizing inductance
Lr = Llr+Lm;    % pu, rotor inductance
Ls = Lls+Lm;    % pu, stator inductance

H = 0.2734;     % s, moment of inertia
F = 0.0106;     % pu,friction coefficient
p = 2;          % pole pairs

Vbase = Vn/sqrt(3)*sqrt(2); % V,     base voltage, peak, line-to-neutral
Ibase = Pn/(1.5*Vbase);     % A,     base current, peak
Zbase = Vbase/Ibase;        % ohm,   base resistance
wbase = 2*pi*fn;            % rad/s, base elec. radial frequency
Tbase = Pn/(wbase/p);       % N*m,   base torque
psin  = (Vn/sqrt(3)*sqrt(2)/wbase);  %Nominal flux

Rss = Rs*Zbase;  % ohm, stator resistance
Xls = Lls*Zbase; % ohm, stator leakage reactance
Rrr = Rr*Zbase;  % ohm, rotor resistance, referred to the stator side
Xlr = Llr*Zbase; % ohm, rotor leakage reactance, referred to the stator side
Xm = Lm*Zbase;   % ohm, magnetizing reactance

%% Control parameters
Ts = 5e-6;        % s,  fundamental sample time
fsw = 2e3;        % Hz, switching frequency 
Tsc = 1/(fsw*10); % s,  control sample time

% Velocity controller parameters
Kp_wr = 65.47;
Ki_wr = 3134.24;