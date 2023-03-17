clear all;
clc;

V_all=0.1:0.1:8;
L=610;
DataAll=zeros(length(V_all),L);
DataEAll=zeros(length(V_all),L);

figure()
for k=1:length(V_all)

V=V_all(k);    
theta=2*pi*rand();
alpha=(sqrt(5)-1)/2;


a=0.5;

X=1-a*cos(2*pi*(theta+alpha*linspace(1,L,L)));

J=1;
lambda=ones(1,L);
H=J*diag(ones(1,L-1),-1)+J*diag(ones(1,L-1),1)+V*diag(lambda./X);

% ±ß½çÌõ¼þ
bd=0;
if bd==0
H(1,L)=0;
end
if bd==1
t1=J;
H(1,L)=t1;
H(L,1)=t1;
end

[Ev,E]=eig(H,'vector'); 
Data=zeros(1,L);
for i=1:L
   Data(i)=Ipr(Ev(:,i)); 
end

DataAll(k,:)=-log(Data')/log(L);
DataEAll(k,:)=E';

end

EVIprPlot(V_all,DataEAll,DataAll)


function a2 = Ipr(psi)
%IPR get the Ipr for a vector \sum_i|psi_i|^4
%  
a=psi.*conj(psi);
a2=sum(a.^2);
end

