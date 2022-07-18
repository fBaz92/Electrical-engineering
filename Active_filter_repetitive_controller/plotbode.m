function[]=plotbode(G,f)
f1=f(1);
f2=f(length(f));
S=size(G);
Gz=zeros(length(f),S(2));
for count=1:length(f)
    Gz(count,:)=freqresp(G,2*pi*f(count));
end
figure
subplot(2,1,1)
semilogx(f,20*log10(abs(Gz))) 
xlim([f1,f2])
title('Bode Diagram')
ylabel('Magnitude  [dB]')
grid on
subplot(2,1,2)
semilogx(f,angle(Gz)/pi*180) 
xlim([f1,f2])
xlabel('Frequency  [Hz]')
ylabel('Phase  [deg]')
grid on

end