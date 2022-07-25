%   Calcolo le caratteristiche di carico di un veicolo in funz 
%   della velocità

%Parametri veicolo

%% Masse

%mc0 curb weight without batteries [kg]
mc0=1342;

%mb: battery mass
mb=350;

%mc: curb weight
mc=mc0+mb;

%mh human mass
mh=80;

%np passengers number
np=1;

%mp: play load
mp=0;

%mg: gross weight
mg=np*mh+mp+mc;

%% VEICOLO
%rt: rapporto di trasmissione (gear ratio) []
rt=4.3;

%rw: Raggio ruota (wheel radius) [m]
%   205/60 R16 92 H   x1/x2 Rx3
x1=205;
x2=60;
x3=16;
rw=(x3*25.4/1000)/2+x1*x2/100000;

%% Aerodinamics %%
%è stato trovato direttamente cd*af=0.652 quindi abbiamo messo af a 1.
%af: frontal area [m^2]
af=1;

%cd: aerodynamic drag coefficient [ ] 
cd=0.652;

%rho: Air density [kg/m^2]
rho=1.2041;

%% Rolling resistance

% SONO STATE SCELTE LE MICHELIN 205/60 DI CLASSE C

%crr_r: coefficient of rolling resistance rear tyre[kg/kg] 
crr_r=0.0085;

%crr_f: coefficient of rolling resistance front tyres [kg/kg]
crr_f=0.0085;

%% Gravity

%g: gravity acceleration
g=9.81;

%% COG
%a: distance between Center Of Gravity and front axle [m]
a=1;
%b: distance between Center Of Gravity and rear axle [m]
b=1;
%l: wheelbase
l=a+b;

%% ROAD
%i: road grad [pu]
i=0;
%alfa: road angle
alfa=atan(i);

%% SPEED VECTOR
%dv: passo con cui faremo i calcoli (vehicle speed step) [m/s]
dv=0.1;

%vehicle vector speed [m/s]
v=[0:dv:40];

%vk vehicle speed vector [km/h]
vk=v*3600/1000;

%% FORCE CALCULATION

% AERODINAMICS
%fa: aerodinamics drag force [N]
fa=1/2*rho*cd*af*v.^2;              %attenzione .^ 
%pa aerodinamic drag power
pa=fa.*v;

% ROLLING RESISTANCE
%frr_f: rolling resistance front tyres
frr_f=mg*g*b/l*cos(alfa)*crr_f;
%frr_r: rolling resistance rear tyres
frr_r=mg*g*a/l*cos(alfa)*crr_r;
%frr: total rolling resistance
frr=(frr_r+frr_f)*ones(size(v)); %ones(v) serve per trasformare la frr in un vettore
%prr: rolling resistance power[W]
prr=frr.*v;

% ROAD GRADE
%fgrade: forceto overcome road slope [N]
fgrade=mg*g*sin(alfa)*ones(size(v));
%pgrade: grade power [W]
pgrade=fgrade.*v;

%%TOTAL
%fl: total load force
fl=fa+frr+fgrade;
%pl: total load power [w]
pl=fl.*v;

%gs: Specifical Energy Consumption [Wh/k]
gs=pl./vk;


%% Grafici %%
figure;
subplot(2,1,1);
plot(vk,fl,vk,fa,vk,frr,vk,fgrade);
grid;
xlabel('Vehicle Speed [km/h]');
ylabel('Force [N]');
legend('Total','Aero','Rolling','Grade');

subplot(2,1,2);
plot(vk,pl/1000,vk,pa/1000,vk,prr/1000,vk,pgrade/1000);
grid
xlabel('Vehicle Speed [km/h]');
ylabel('Power [kW]');
legend('Total','Aero','Rolling','Grade');

figure;
plot(vk,gs);
grid;
xlabel('Vehicle Speed [km/h]');
ylabel('Specific consumption [Wh/km]');
legend('Specific Consumption');



