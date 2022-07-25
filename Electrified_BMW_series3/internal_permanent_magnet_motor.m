%IPM da circa 12kW
%specifiche di macchina sincrona brushless a magneti permanenti
%superficiali
clear
close all
%Parametri di macchina riferiti al s.d.r. bifasico rotante, k=2/3
% Ldisp: induttanza di dispersione [H]
Ldisp = 0.0e-3; 
% Lmd : induttanza magnetizzante sull'asse d [H] .......sarebbe: 3/2(Lmo+Lm).......
Lmd = 0.156e-3;
% Lmd : induttanza magnetizzante sull'asse q [H] .......sarebbe: 3/2(Lmo-Lm).......
Lmq = Lmd*2e-3;
% rapporto: rapporto di trasformazione rotore/statore N2/N1
% rapporto = 7;
% Ld: induttanza fittizia asse d [H]
Ld = Lmd + Ldisp
% Lq: induttanza fittizia asse q [H]
Lq = Lmq + Ldisp
% D: differenza tra le induttanze dq [H]
D = Ld-Lq;
% M: mutua induttanza tra fase di rotore e statore [H]
% M = 2/3*rapporto*Lmd;
%Inn: corrente nominale [A]
Inn=677.587;
% p: numero di coppie di poli
p = 4;
%Vdc: tensione di bus dc inverter [V]
%Vdc = 540;
%Vo: modulo max del vettore tensione [V]
%Vo=Vdc/sqrt(3)
%Vdc=420;
%Vo=Vdc/1.73
%%%(scelta di 320Vdc quindi 
Vo=185.44;
% ieM: corrente di eccitazione massima ammissibile
%ieM=8.6
%ieM=4.3;
%ieM=2.1;
% valore massimo del parametro a 
%a=1 -> fluxM=Ld*Id
%a<1 -> fluxM<Ld*Id centro delle ellissi di tensione interno al cerchio limite di corrente
%a>1 -> fluxM>Ld*Id centro delle ellissi di tensione esterno al cerchio limite di corrente

% r: rapporto tra le induttanze di asse d ed asse q (eccentricità dell'ellisse)
r=Ld/Lq
r1=1-1/r;
%fluxM: flusso prodotto dai magneti concatenato con lo statore
%flux M può essere calcolato dai dati di targa cioè: ke, p
%fluxM=V0/wn=VLL_RMS*sqrt(2/3)/p/w_mecc=ke*sqrt(2/3)/p
%%
%fluxM può essere calcolato utilizzando una corrente di eccitazione
%fittizia (o reale nel caso di WRSM)
%fluxM=M*ieM;
fluxM=Ld*Inn*0.9;
%a: rapporto tra flusso dei magneti e reazione sullasse d con la corrente nominale
a=fluxM/(Ld*Inn)
% id0: centro dell'ellisse
id0=-a

%% valutazione dei coefficienti per le perdite nel ferro e nel rame  


% pfe   perdite nel ferro/potenza nominale:
pfe_x=0.03;
% pfe_i_pfe :perdite per isteresi rispetto ale perdite totali nel ferro
pfe_i_pfe=0.2;
%pcu_x perdite nel rame rispetto alla potenza nominale 
pcu_x=0.06;
B_rated=1.35;