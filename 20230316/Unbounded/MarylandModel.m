% Maryland Model

clear all;
clc;


theta=2*pi*rand();
alpha=(sqrt(5)-1)/2;

L=100;
X=diag(tan(2*pi*(theta+alpha*linspace(1,L,L))));

J=1;
lambda=1;
H=J*diag(ones(1,L-1),-1)+J*diag(ones(1,L-1),1)+lambda*(X);

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

plot(linspace(1,L,L)/L,E,'.')
hold on;

figure()
mesh(Ev.*conj(Ev))

figure()
plot(X)
