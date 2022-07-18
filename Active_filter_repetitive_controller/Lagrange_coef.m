function h=Lagrange_coef(f1,fc)

h=NaN(4,round(fc/f1*0.5)); %mi calcolo h per ogni freq fino a fc/2
D=zeros(round(fc/f1*0.5),1);
for count=1:round(fc/f1*0.5)
    D(count)=fc/(count*f1)-fix(fc/(count*f1))+1;
    h(1,count)=-(D(count)-1)*(D(count)-2)*(D(count)-3)/6;
    h(2,count)=D(count)*(D(count)-2)*(D(count)-3)/2;
    h(3,count)=-D(count)*(D(count)-1)*(D(count)-3)/2;
    h(4,count)=D(count)*(D(count)-1)*(D(count)-2)/6;
end
