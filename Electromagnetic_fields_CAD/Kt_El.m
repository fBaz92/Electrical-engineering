function [K_El,t_El] = Kt_El (problema, nodi, Prop)
% CALCOLO DELLA MATRICE DI ELEMENTO K_El E DEL TERMINE NOTO DI ELEMENTOO
% t_El

%problema=1 : PROBLEMA ELETTROSTATICO 1D
%problema=2 : PROBLEMA MAGNETICO 1D
%problema=3 : PROBLEMA ELETTROSTATICO ASSIALSIMMETRICO 1D IN r 
%problema=4 : PROBLEMA MAGNETICO STAZIONARIO ASSIALSIMEMTRICO 1D IN r

switch problema
    case {1 2}
        K_El= Prop (2)/Prop(1) * [1, -1;-1, 1];
        t_El= 0.5 * Prop (1) * Prop (3) * ones(2,1);
    case 3
        rm=sum(nodi) * 0.5
        K_El= rm * Prop (2)/Prop(1) * [1, -1;-1, 1];
        t_El= Prop (1) * Prop (3)/6 * [2, 1;1, 2] * nodi;
    case 4
        disp('non ancora implementato')
end






end


