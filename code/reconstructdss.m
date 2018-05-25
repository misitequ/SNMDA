clc
clear
for i=1:383
   load knownre
   knownre=knownre';
   testsample = knownre(:,i);
   knownre(:,i)=[];
   trainsample = knownre;
Sxp_I(:,i)=computaccuracy(trainsample,testsample);
end
Sxp_I1=abs(Sxp_I);
a=zeros(1,383);
b=Sxp_I;
c1=[a;b];
for i=1:383
    for j=1:i-1
        c1(j,i)=c1(j+1,i);
    end
    c1(i,i)=0;
end
Sxp_I=c1;
ind = Sxp_I>1e-5;
cel = cell(1,383);
load knownre
A = knownre';
for i = 1:383
    a = ind(:,i);
    B = A(:,a);
    xp=l1_ls(B,A(:,i),0.01,0.01);
    cel{1,i} = xp;
end
C = zeros(383,383);
for i = 1:383
    ind1 = ind(:,i);
    cell = cel{:,i};
    C(ind1,i) = cell;
end
dss=abs(C)
redss=(dss+dss')/2;