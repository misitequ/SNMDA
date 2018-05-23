clear
clc
load mfs;
load dss
load remfs;
load redss;
load HMDD;
load knownre;
A=HMDD;
nd=max(A(:,1)); 
nm=max(A(:,2));
[pp,qq]=size(A);
 for i=1:495
     for j=1:495
         if mfs(i,j)==0;
             mfs1(i,j)=remfs(i,j);
        else
             mfs1(i,j)=(mfs(i,j)+remfs(i,j)/2);
        end
     end
 end
 for i=1:383
    for j=1:383
         if dss(i,j)==0;
             dss1(i,j)=redss(i,j);
         else
             dss1(i,j)=(dss(i,j)+redss(i,j)/2);
         end
     end
 end
M1=sum(mfs1);
  for i=1:495
      for j=1:495
          mfs1(i,j)=mfs1(i,j)/(((M1(i)*M1(j))^0.5));
      end
  end
 D1=sum(dss1);   
 for i=1:383
     for j=1:383
         dss1(i,j)=dss1(i,j)/(((D1(i)*D1(j))^0.5));
     end
 end
gamma=0.5;
PM=knownre';
PD=knownre;
P0=knownre';
k= 0;
delta = 1;
while  (delta > 1e-6)
    PM1 = (1-gamma)*mfs1*PM+gamma*P0;
    delta =abs(sum(sum((abs(PM1)-abs(PM)))));
    PM = PM1;
end
delta = 1;
while  (delta > 1e-6)
    PD1 = (1-gamma)*dss1*PD+gamma*P0';
    delta =abs(sum(sum((abs(PD1)-abs(PD)))));
    PD =PD1;
end
F=(PM1+PD1')';
for i=1:nd
    for j=1:nm
        if knownre(i,j)==1
           F(i,j)=-10000;
        end
    end
end
finalprediction=[];
for i=1:nm
    finalprediction=[finalprediction;F(:,i)];
end
finalpredictiondisease=[];
for i=1:nm
    finalpredictiondisease=[finalpredictiondisease,1:nd];
end
finalprediction(:,2)=finalpredictiondisease';
finalpredictionmicrobe=[];
for i=1:nm
    finalpredictionmicrobe=[finalpredictionmicrobe;i.*ones(nd,1)];
end
finalprediction(:,3)=finalpredictionmicrobe;
discard=finalprediction(find(finalprediction(:,1)==-10000),:);
finalprediction=setdiff(finalprediction,discard,'rows');
result=sortrows(finalprediction,-1);
[a1,~,a]=xlsread('diseasename');
[b1,~,b]=xlsread('miRNAname');
m=length(result);
c=zeros(m,1);
f=zeros(m,1);
d={};
e={};
for i=1:m
if(~isempty(find(a1==result(i,2))))
c(i,1)=find(a1==result(i,2));
d{i,1}=a{c(i,1),2};
else
d{i,1}=result(i,2);
end
end
for i=1:m
if(~isempty(find(b1==result(i,3))))
f(i,1)=find(b1==result(i,3));
e{i,1}=b{f(i,1),2};
else
e{i,1}=result(i,3);
end
end
save result result
