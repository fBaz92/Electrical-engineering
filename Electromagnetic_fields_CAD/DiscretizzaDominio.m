function [ npD, ndivD, puntipD, CodCC, ValCC, ElementiD, PropElD ] = DiscretizzaDominio( nr, Regioni, divReg, CodCCR, ValCCR)
%funzione per discretizzare il dominio di calcolo 
%grandezze in uscita tra paretesi quadre 
%grandezze in ingresso tra parentesi tonde
%Nr=numero di regioni
%   Detailed explanation goes here


ndivD = sum(divReg);
npD = int16(ndivD + 1);
puntipD = zeros(npD, 1);
ElementiD = zeros(ndivD, 2,'int16');
PropElD = zeros(ndivD, 3);

CodCC = zeros(npD, 1, 'int16');
ValCC= zeros(npD, 1);


[ puntip, ElementiR, PropElR ] = DiscretizzaRegione( Regioni(1, :), divReg(1));
ipcorr = divReg(1) + 1;
ielcorr = divReg(1);

puntipD(1:ipcorr) =  puntip;
ElementiD(1:ielcorr, :)  = ElementiR;
PropElD(1:ielcorr, :)  = PropElR;
CodCC(1) = CodCCR(1,1);
CodCC(ipcorr) = CodCCR(1,2);

ValCC(1) = ValCCR(1,1);
ValCC(ipcorr) = ValCCR(1,2);

for iR=2:nr
   [puntip, ElementiR, PropElR] = DiscretizzaRegione( Regioni(iR, :), divReg(iR));
   
   ipnext = ipcorr + divReg(iR) ;
   ielnext = ielcorr + divReg(iR) ;
   puntipD(ipcorr+1: ipnext) =  puntip(2:divReg(iR)+1);
   ElementiD(ielcorr + 1:ielnext, :)  = ElementiR + ielcorr;
   PropElD(ielcorr + 1:ielnext, :)  = PropElR;
   if CodCCR(iR,1) > CodCC(ipcorr)
      CodCC(ipcorr) = CodCCR(iR,1);
      ValCC(ipcorr) = ValCCR(iR,1);
   end
   CodCC(ipnext) = CodCCR(iR,2);
   ValCC(ipnext) = ValCCR(iR,2);
   ipcorr = ipnext;
   ielcorr = ielnext;
      
end
   





end

