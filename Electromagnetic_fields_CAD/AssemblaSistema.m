function [ K_glob, t_glob ] = AssemblaSistema( npD, ndivD, ElementiD, PropElD, CodCC, ValCC, puntipD, problema)
% ASSEMBLAGGIO DELLA MATRICE DEI COEFFICIENTI E DEL TERMINE NOTO DEL
% SISTEMA RISOLVENTE

K_glob=zeros (npD);
t_glob=zeros (npD,1);

for iel = 1 : ndivD
    nodi=[puntipD(ElementiD(iel,1));puntipD(ElementiD(iel,2))];
    [K_El,t_El] = Kt_El(problema, nodi,PropElD(iel,:));
    
    %righe per l'assemblaggio
    
for i = 1:2
  irg= ElementiD(iel,i);
  if CodCC(irg) == 2;
      K_glob(irg,irg) = 1;
      t_glob(irg)=ValCC(irg);
  else
      for j = 1:2
          icl= ElementiD(iel,j);
          K_glob(irg,icl) = K_glob(irg,icl) + K_El (i,j);
      end
      t_glob(irg)=t_glob(irg) + t_El(i);
      if CodCC(irg) == 1
          t_glob(irg) = t_glob(irg) + PropElD(iel,2) + ValCC(irg);
      end
  end
end
end
end
