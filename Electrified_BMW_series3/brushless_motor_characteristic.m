%caratteristica_macchina_brushless_i_variabile_senza_dati_motore.m
% *************************** cerca l'equazione della traiettoria in deflussaggio**********************************************
% per una macchina brushless a MP (derivato dal sincrono con eccitazione)

%%richiamo il file dove viene definito il motore

internal_permanent_magnet_motor
%deltai=1
deltai=[0.2 : 0.2 :2];
In = Inn.*deltai
ir=find(deltai==1);


aM=fluxM/Ld./In

% r: rapporto tra le induttanze di asse d ed asse q (eccentricità dell'ellisse)
r=Ld/Lq
r1=1-1/r;
% eseguo i calcoli per ciascun valore di velocità
n= [10:50:12000];
w=p*2*pi*n/60;
%b: parametro rapporto (semiasse maggiore ellisse/raggio cerchio)
% i: indice del vettore correnti
% k: indice del vettore velocità
for  i=1:length(In)


%
%
%%% tratto a bassa velocità:  da w=0 ad w=wB
%
% x1: angolo thetaB a cui si ha la coppia massima erogabile con ie=ieM
x1(i)=(-aM(i)+sqrt(aM(i)^2+8*(1-1/r)^2))/(4*(1-1/r));

thetaB(i)=acos(x1(i));
thetaB_ang(i)=thetaB(i)*180/pi;
%velocità wB fino a cui si ha la coppia massima
% calcolo la velocità base a cui avviene l'intersezione tra MTC e cerchio
wB(i)=Vo/Ld/In(i)/sqrt(1/r^2*(1-x1(i)^2)+(x1(i)+aM(i))^2);
nB(i)=wB(i)*60/2/pi/p;
%deltaB: angolo delta (tra tensione e fem) nel punto a coppia massima
deltaB(i)=atan(1/r*(sin(thetaB(i)))/(aM(i)+cos(thetaB(i))));
deltaB_ang(i)=deltaB(i)*180/pi;
%fiB: angolo fi tra tensione e corrente nel punto B 
fiB(i)=pi/2-thetaB(i)+deltaB(i);
%coppia erogata nel punto base B
cB(i)=1.5*p*Ld*In(i)^2*(aM(i).*sin(thetaB(i))+(1-1/r).*sin(thetaB(i)).*cos(thetaB(i)));
pB(i)=wB(i)/p*cB(i);
% semiasse (lungo la direzione di iq) dell'ellisse limite di tensione   
b(i,:)=Vo/Lq./In(i)./w;
%% trovo l'indice del vettore w delle velocità per cui w<=wB
iwB=find(w<=wB(i));
iiwB=length(iwB);
c_bassa(i,1:iiwB)=cB(i);
w_bassa(i,1:iiwB)=w(1:iiwB);
a_bassa(i,1:iiwB)=aM(i);
theta_bassa(i,1:iiwB)=thetaB(i);
fi_bassa(i,1:iiwB)=fiB(i);
fluxM_bassa(i,1:iiwB)=fluxM;
pippo(i)=iiwB;
v_bassa(i,1:iiwB)=w_bassa(i,1:iiwB).*sqrt((fluxM_bassa(i,1:iiwB)+Ld.*In(i).*cos(theta_bassa(i,1:iiwB))).^2+(In(i).*Lq.*sin(theta_bassa(i,1:iiwB))).^2);
p_bassa=c_bassa.*w_bassa/p;
%% trovo l'indice del vettore w delle velocità per cui w>wB
ind_media=find(w>wB(i));
iiimb=length(ind_media);
w_media(i,1:iiimb)=w(iiwB+1:iiwB+iiimb);
a_media(i,1:iiimb)=aM(i);
b_media(i,1:iiimb)=Vo/Lq./In(i)./w_media(i,1:iiimb);
xM(i,1:iiimb)=(-r^2*aM(i)+sqrt(-r^2+r^2.*b_media(i,1:iiimb).^2+1-b_media(i,1:iiimb).^2+r^2*aM(i)^2))/(r^2-1);
theta_media(i,1:iiimb)=acos(xM(i,1:iiimb));
c_media(i,1:iiimb)=1.5*p*Ld.*In(i)^2*(aM(i).*sin(theta_media(i,1:iiimb))+(1-1/r).*sin(theta_media(i,1:iiimb)).*cos(theta_media(i,1:iiimb)));
fluxM_media(i,1:iiimb)=fluxM;
%deltaBC: angolo delta (tra tensione e fem) nel tratto di velocità media
deltaBC(i,1:iiimb)=atan(1/r.*(sin(theta_media(i,1:iiimb)))./(a_media(i,1:iiimb)+cos(theta_media(i,1:iiimb))));
fi_media(i,1:iiimb)=(pi/2-theta_media(i,1:iiimb)+deltaBC(i,1:iiimb));
%%%

v_media(i,1:iiimb)=w_media(i,1:iiimb).*sqrt((fluxM_media(i,1:iiimb)+Ld.*In(i).*cos(theta_media(i,1:iiimb))).^2+(In(i).*Lq.*sin(theta_media(i,1:iiimb))).^2);

w_bassa_media(i,:)=[w_bassa(i,1:iiwB) w_media(i,1:iiimb)];
a_bassa_media(i,:)=[a_bassa(i,1:iiwB) a_media(i,1:iiimb)];
theta_bassa_media(i,:)=real([theta_bassa(i,1:iiwB) theta_media(i,1:iiimb)]);
c_bassa_media(i,:)=real([c_bassa(i,1:iiwB) c_media(i,1:iiimb)]);
fluxM_bassa_media(i,:)=[fluxM_bassa(i,1:iiwB) fluxM_media(i,1:iiimb)];
fi_bassa_media(i,:)=real([fi_bassa(i,1:iiwB) fi_media(i,1:iiimb)]);
v_bassa_media(i,:)=real([v_bassa(i,1:iiwB) v_media(i,1:iiimb)]);
end
n_bassa_media=w_bassa_media*60/2/pi/p;
cos_fi_bassa_media=cos(fi_bassa_media);
p_bassa_media=c_bassa_media.*w_bassa_media/p;
%deltaBC_ang=deltaBC*180/pi;

% grafici in funzione della velocità
figure
subplot(3,1,1)
plot(n_bassa_media',c_bassa_media','Linewidth',2)
ylabel('torque [Nm]')
axis([0 9000 0 200])
grid
subplot(3,1,2)
plot(n_bassa_media',v_bassa_media','Linewidth',2)
ylabel('voltage [V]')
axis([0 9000 0 400])
grid
subplot(3,1,3)
plot(n_bassa_media',theta_bassa_media'*180/pi,'Linewidth',2)
ylabel('theta [degrees]')
axis([0 9000 80 180])
grid
legend('0.2','0.4','0.6','0.8','1','1.2','1.4','1.6','1.8','2','2.2','2.4')


figure 
subplot(2,1,1)
plot(n_bassa_media',c_bassa_media','Linewidth',2)
%legend('2.4','2.2','2','1.8','1.6','1.4','1.2','1','0.8','0.6','0.4','0.2')
legend('0.2','0.4','0.6','0.8','1','1.2','1.4','1.6','1.8','2','2.2','2.4')
ylabel('torque [Nm]')
axis([0 12000 0 1500])
grid
subplot(2,1,2)
plot(n_bassa_media',p_bassa_media'/1000,'Linewidth',2)
ylabel('power [kW]')
axis([0 12000 0 500])
grid
xlabel ('speed [rpm]')
%%% tratto a velocità intermedia da wB a wC, 
%   calcolato dall'intersezione ellisse-cerchio nel tratto intermedio con ie=ieM



%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
%%%%%% stima delle perdite %%%%%%%
%%%%%%
%% valutazione delle perdite nel ferro e nel rame  
%coefficienti definiti nel file di definizione del motore
% pfe   perdite nel ferro:
%%%
%%% in condizioni nomiali ipotizzo Pfe=pfe_x*Prated
%l'elemento n.5 è quello con la corrente nominale%%%
pfe_rated=pfe_x.*max(p_bassa(ir,:))
%%%%%%%%%%%%%% sono arrivato qui%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%

%%%  in condizioni nominali ipotizzo perdite per isteresi: Pfe_ist =pfe_i_pfe*Pfe
pfe_ist_rated=pfe_i_pfe*pfe_rated;
%%%  in condizioni nominali ipotizzo perdite per correnti parassite: Pfe_cp =60% Pfe
pfe_cp_rated=pfe_rated-pfe_ist_rated;
%pfe_i=ki*w*B^2.2  
ki=pfe_ist_rated./wB(ir)/B_rated^2.2
kcp=pfe_cp_rated./(wB(ir)*B_rated).^2;
%%%
%pcu perdite nel rame: Pcu=pcu_x*Prated
pcu_rated=pcu_x.*max(p_bassa(ir,:))
kcu=pcu_rated./In(ir).^2;
%%%%%%%%%%%%%%%%%%%%

% calcolo perdite
% ipotizzo B costante, NO saturazione
B=B_rated;
for i=1:length(In)
pfe(i,:)=ki.*w*B^2.2+kcp.*w.^2*B^2;
pcu(i,:)=kcu.*In(i).^2.*w./w;
end
ptot=pfe+pcu;
pin=p_bassa_media+ptot;
eta=p_bassa_media./pin;
figure 
plot(n_bassa_media',eta','Linewidth',2)
ylabel('efficiency ')
axis([0 3000 0 1])
grid
legend('0.2','0.4','0.6','0.8','1','1.2','1.4','1.6','1.8','2','2.2','2.4')




%figure

c_bassa_media_pu=c_bassa_media./cB(ir);
w_bassa_media_pu=w_bassa_media./wB(ir);


%%%disegno la mappa di efficienza sul diagramma coppia-velocità in pu

% Create figure
figure1 = figure;
colormap('jet');

% Create axes e definisco il limite sull'asse z della scala colorata
axes1 = axes('Parent',figure1,'CLim',[0.6 0.93]);
grid(axes1,'on');
hold(axes1,'on');

% Create mesh
%mesh(w_bassa_media_pu',c_bassa_media_pu',eta','Parent',axes1,'Facealpha',0.9,'Edgealpha',0.5);
mesh(w_bassa_media_pu',c_bassa_media_pu',eta','Parent',axes1);
% Create colorbar
colorbar('peer',axes1);
title ('EFFICIENCY MAP')

%%trova l'inviluppo delle coppie per trovare la coppia max. 
% ad alta velocità non è detto che sia la stessa curva a bassa
for i=1:length(w)
c_bassa_media_max(i)=max(c_bassa_media(:,i));
end
hold on
%disegno il profilo di coppia con la corrente nominale e max.


plot(w_bassa_media_pu(ir,:)',c_bassa_media(ir,:)./cB(ir)','--','Linewidth',2)
%,n_bassa_media(length(deltai),:)',c_bassa_media(length(deltai),:)','--','Linewidth',2)
plot(w_bassa_media_pu(1,:)',c_bassa_media_max(1,:)./cB(ir)','-g','Linewidth',3)
xlabel ('speed [pu]')
ylabel ('torque [pu]')
 
hold on

%mesh(n_bassa_media',c_bassa_media',p_bassa_media','Facealpha',0.5,'Edgealpha',0.3)
xlabel ('speed [rpm]')
ylabel ('torque [Nm]')

%%%%cambio di variabili
t_out=c_bassa_media;
n_out=n_bassa_media;
p_in=pin;
losses=ptot;


%%%PREPARAZIONE DATO PER LOOK-UP TABLE SIMULINK%%%%

figure
mesh(n_out,t_out,p_in)
n_data=[10:100:10100];
t_data=[1:1:max(max(t_out))];
% per utilizzare la look-up table 2D di Simulink occorre preparare i dati in
% modo che:
% row index inpout value sia monotono crescente

% column index input value sia monotona crescente
% table data non contengano NaN
% Utilizzo il comando griddata, con l'opzione 'nearest' in modo che i punti
% al di sopra della potenza massima siano uguagliati alla potenza massima.
% tanto in quella zona non ci si va.

xlabel('Speed [rpm]')
ylabel('Torque [Nm]')

%p_in_data=griddata(n_out,t_out,p_in,n_data,t_data','nearest');
losses_data=griddata(n_out,t_out,losses,n_data,t_data','nearest');
figure
mesh(n_data,t_data,losses_data)

losses_data=losses_data';
xlabel('Speed [rpm]')
ylabel('Torque [Nm]')

%preparo la caratteristica limite da utilizzare nelle simulazioni in
%Simulink
n=max(n_out);
c=max(t_out);
hold on
figure 
plot(n,c,n_out(5,1:240),t_out(5,1:240),'Linewidth',3)
grid




% save C:\Users\matte\Desktop\DATI_ROSSI\mappa2016 n_data t_data losses_data n c


