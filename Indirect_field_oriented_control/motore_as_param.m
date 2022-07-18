% FILE MOTORE_AS_PARAM.M (C) L. Zarri
% INIZIALIZZA I PARAMETRI NECESSARI PER LA
% SIMULAZIONE DEL MOTORE ASINCRONO
% Ultimo aggiornamento: Maggio 2008

% PARAMETRI RELATIVI AL MOTORE ASINCRONO

% Induttanze (H)
Ls = 123.95e-3;         % Induttanza statorica
Lr = 123.95e-3;         % Induttanza rotorica
M = 118e-3;             % Induttanza mutua
sigma =1-M*M/(Ls*Lr);   % Coefficiente di accoppiamento magnetico

% Resistenze (Ohm)
Rs=2.78;                % Resistenza statorica
Rr=1.137;               % Resistenza rotorica

% Coppie di poli
p=2;

% Coppia nominale (Nm)
Cn=27;

% Tensione nominale (concatenata rms)
Vn = 400;

% Corrente nominale (A rms)
In=9.40;

% Frequenza nominale di alimentazione (Hz)
fn = 50;

% Velocità meccanica nominale in rpm e in rad/s
Nn = 1440;
Wn = Nn*2*pi/60;

% Scorrimento nominale
sn= 1 - Wn*p/(2*pi*fn);

% Flusso rotorico nominale (Wb)
Frn = sqrt((2./3)*Rr*Cn/(p*sn*2*pi*fn));

% Corrente di asse d nominale (A di picco)
Idn = Frn/M;

% Corrente di asse q nominale (A di picco)
Iqn = (2/3)*Cn*(Lr/M)/(p*Frn);

% Inerzia del motore (kg m^2)
Jm = 0.01330;

% PARAMETRI RELATIVI ALL'INVERTER

% Tensione del bus DC (V)
Edc=700;

% Periodo di commutazione del chopper (s)
Tc = 200e-6;

% PARAMETRI RELATIVI AL CARICO

% Momento d'inerzia del carico (kg m^2)
Jc = 3 * Jm; % Supponiamo che l'inerzia del carico sia il triplo del motore.

% Coefficiente di attrito (Nm s/rad)

b=1e-3;

% Costante di coppia resistente
Kc = Cn/Wn-b; % Supponiamo che la scelta del motore sia perfetta, cioè in condizione di regime,
              % la velocità è pari a quella nominale del motore.
            
% PARAMETRI RELATIVI AL CONTROLLO

Rsn=Rs;
Rrn=Rr;
Lsn=Ls;
Lrn=Lr;
Mn=M;
sigman=1-Mn*Mn/(Lsn*Lrn);

% Corrente massima (A di picco)

Imax = 20;

% STUDIO DEL CONTROLLO

J=Jc+Jm;            % Inerzia totale;

s=tf('s');

% Controllo di corrente

tau_m=J/b;              % Costante meccanica
tau_er=Lr/Rr;           % Costante di tempo rotorica
tau_es=Ls/Rs;           % Costante di tempo statorica

Mot=1/(sigma*Ls*s+Rs);    % FdT del motore brushless senza azione della fem
Del=1/(1+0.5*Tc*s);     % Fdt del convertitore

Gvi=Del*Mot; % Funzione di trasferimento Vref->I

% Parametri del regolatore di corrente (MF = 75°, Wc=2620 rad/s) 
tau_i=sigman*Lsn/Rsn;
Ki_i=7664;
Kp_i=Ki_i*tau_i;

% Controllo di velocità

Wc=2620;
Giw=1.5*p*(M/Lr)*Frn*Del*(1/(J*s+b))*1/(1+s/Wc);

tau_w=0.31;
Ki_w=35.65;

tg=1; %periodo gradino
