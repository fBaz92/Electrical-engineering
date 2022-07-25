function [ puntip, ElementiR, PropElR ] = DiscretizzaRegione( regione, ndivR)
%discretizza una regione
%   Detailed explanation goes here

LR = regione(2) - regione(1);
dx = LR / ndivR;
puntip = zeros(ndivR + 1, 1);
for ip = 0 : ndivR 
    puntip(ip+1, 1) = regione(1) + LR * ip/ndivR;
end


iini = int16(1:ndivR);
iilast = iini + int16(1) ;


ElementiR = [iini', iilast'];
PropElR = [dx*ones(ndivR , 1), regione(3)*ones(ndivR , 1), regione(4)*ones(ndivR , 1)];




end

