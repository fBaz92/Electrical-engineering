% Calcolo dello spettro complesso di un segnale

clear
close all
clc

% Importa i dati (insieme di tre segnali)
% y1=importdata('1805_prova3_3corr.txt');
[FileName,PathName] = uigetfile('*.csv','Select csv file saved from the scope');
y1=csvread(strcat(PathName,FileName),17,1);
temp=csvread(strcat(PathName,FileName),0,0,'B8..B8');
SampleRate=temp(2,2);
DisplayBlockSize=length(y1);
y2(:,1:3)=y1(:,1:3);    % 2:4 per correnti, 3:5 per tensioni

y2=y2*10;   % 10 per correnti, 200 per tensioni
y2(:,3)=-y2(:,1)-y2(:,2);     % Usare se manca uno dei tre segnali

y2(:,1)=y2(:,1)-sum(y2(:,1))/length(y2(:,1));
y2(:,2)=y2(:,2)-sum(y2(:,2))/length(y2(:,2));
y2(:,3)=y2(:,3)-sum(y2(:,3))/length(y2(:,3));
% y2(:,3)=y2(:,2)-y2(:,1);

% Vettori di spazio del sistema trifase
yv(:,1)=2/3*(y2(:,1)*cos(0*2*pi/3)+y2(:,2)*cos(1*2*pi/3)+y2(:,3)*cos(2*2*pi/3));
yv(:,2)=-2/3*(y2(:,1)*sin(0*2*pi/3)+y2(:,2)*sin(1*2*pi/3)+y2(:,3)*sin(2*2*pi/3));   % + o - a seconda dei versi

y=yv(:,1)+1i*yv(:,2);

% Frequenza di campionamento [Hz]
% SampleRate=10e4;
% Dimensione vettore dati importati
% DisplayBlockSize=length(y);

% Figura segnale
x=0:1/SampleRate:DisplayBlockSize*1/SampleRate-1/SampleRate;
% y3(:,1:4)=y1(:,2:5);
figure
plot(x,y1(:,1:3))
axis tight
grid on
title('Segnale')

% Figura segnale importato
figure
plot(x,y2)
axis tight
grid on
title('Segnale importato')

% Figura segnale complesso
figure
plot(x,real(y(1:DisplayBlockSize)),x,imag(y(1:DisplayBlockSize)))
axis tight
grid on
title('Segnale complesso')

% Calcolo FFT
var_size=length(y);
pp=hann(length(yv));    % Applicazione della finestra di Hanning
yyv=zeros(length(yv),2);
for count1=1:2
    for count=1:length(yv)
        yyv(count,count1)=yv(count,count1)*pp(count);
    end
end
yy=yyv(:,1)+1i*yyv(:,2);
ft=fft(yy);
Pyypp=abs(ft)/(length(x))*2;
% Pyypp(1)=NaN;
fp=SampleRate/length(x)*(0:length(x)/2-1);

% Figure
figure
semilogy(fp,Pyypp(1:length(x)/2,1)/max(Pyypp),'r');
hold on
semilogy(-fp(2:length(fp)),Pyypp(length(x):-1:length(x)/2+2,1)/max(Pyypp),'r');
axis([-1e3 1e3 1e-3 1.1])
grid on
xlabel('Frequency [Hz]')
ylabel('log Amplitude [dB]')
title('FFT log')
hold off

figure
plot(fp,Pyypp(1:length(x)/2,1)/max(Pyypp),'r','lineWidth',2);
hold on
plot(-fp(2:length(fp)),Pyypp(length(x):-1:length(x)/2+2,1)/max(Pyypp),'r','lineWidth',2);
axis([-1e3 1e3 1e-3 1.01])
grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [-]')
title('FFT lineare')
hold off
%% 

% Calcolo THD
n=2000;
periodi=50;     % periodi=(100200-200)/n;
THDp=zeros(periodi,3);
j=0;
for i=1:periodi
    selezione=y2( 1+j*n : n+1+j*n , :);
    ftTHD=fft(selezione);
    PyTHD=abs(ftTHD)/(length(n))*2;
    THDp(i,:)=sqrt(sum(PyTHD(3:1:41,:).^2))./PyTHD(2,:)*100;
%     THDp(i)=sqrt(sum(PyTHD(3:2:41).^2+sum(PyTHD(1999-41:2:1999).^2)))/PyTHD(2)*100;
    j=j+1;
end
THD=sum(sum(THDp(:,1:3)/3))/periodi;
fprintf('THD = %f\n\n',THD)
figure
PyTHD(1,:)=0;
bar(fp(1:25:25000),PyTHD(1:length(PyTHD)/2,1)./max(PyTHD(:,1)));
ax=gca;
ax.YScale='log';
hold on
f=bar(-fp((26:25:25000-25)),PyTHD(length(PyTHD)-1:-1:length(PyTHD)/2+2,1)./max(PyTHD(:,1)));
axis([-1e3 1e3 1e-3 1.2])
ax=gca;
ax.YScale='log';
grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [-]')
title('FFT lineare')
hold off