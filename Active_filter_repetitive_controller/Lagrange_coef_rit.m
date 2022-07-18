function h_rit=Lagrange_coef_rit(D)

    h_rit(1)=-(D-1)*(D-2)*(D-3)/6;
    h_rit(2)=D*(D-2)*(D-3)/2;
    h_rit(3)=-D*(D-1)*(D-3)/2;
    h_rit(4)=D*(D-1)*(D-2)/6;
end
