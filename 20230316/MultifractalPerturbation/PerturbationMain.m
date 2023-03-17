
PerturbationAnalytic4()

function []=typicalEigenstates()

% We plot typical eigenstates and the enery spectrum

L=300; 
t0=1;
lambda=0.05;
J0=1.5;
H=H1(L,t0,lambda,J0);

[Ev,E]=eig(H,'vector');

m=L;
psi=Ev(:,m);

% plot(log10(psi.*conj(psi)))
% axis([1,L,-log(2*L),0])
% hold on;

subplot(1,2,1)
plot(psi.*conj(psi))
hold on;
subplot(1,2,2)
plot(E,'.')
hold on;

end

function H=H1(L,t0,lambda,J0)
% We consider the PBC case

% J0+cos term
    omega=(sqrt(5)-1)/2;
    h1=kron(J0*ones(1,L-1)+2*t0*cos(2*pi*omega*linspace(1,L-1,L-1)),[1,0]);
    h2=kron(lambda*ones(1,L-1),[0,1]);
    h=h1+h2;
    h=[h,J0+2*t0*cos(2*pi*omega*L)];
    V=kron(J0+2*t0*cos(2*pi*omega*linspace(1,L,L)),[1,1]);
    H=diag(h,-1)+diag(h,1)+diag(V);
%     H(1,2*L)=lambda;
%     H(2*L,1)=lambda;
    
end

% another perturbation
% In this part�� we check the result with another perturbation theory

function []=PerturbationAnalytic3()
% ��֤��΢���۹�ʽ

%parameters
% m,n Ϊѭ��ָ�꣬��ʵ������
% a,b Ϊѭ��ָ�꣬��ʵ������

% ׼������
L=100;
lambda=0.2;
J0=2.5;
t0=1;
omega=(sqrt(5)-1)/2;


% EigenStates and EigenValues

% localized eigenstates
% E !=0
E1=J0*ones(1,L)+2*t0*cos(2*pi*omega*linspace(1,L,L));
Ev1=kron(eye(L),[1;1])/sqrt(2);
% E=0 EigenStates
Ev2=kron(eye(L),[1;-1])/sqrt(2);

% ͨ�����ڷ��̼�������һ��΢�źͲ��������΢��
% ��д�� ΢�ž���
h2=kron(lambda*ones(1,L-1),[0,1]);
h2=[h2,0];
V=diag(h2,-1)+diag(h2,1);
% V(2*L,1)=lambda;
% V(1,2*L)=lambda;

% ���ڷ��̾���
V1=zeros(L,L);
for m=1:L
    for n=1:L
       psi1=Ev2(:,n);
       psi2=Ev2(:,m)';
       V1(m,n)=psi2*V*psi1;       
    end
end
%PBC case �����΢�ź󣬴��ڼ�
%OBC case ����������free hopping

% �����ڷ��̵õ���������׽���
[Ev,E]=eig(V1,'vector');
% �õ��㼶������������һ����������
Ev2a=zeros(2*L,L);
for m=1:L
    for n=1:L
        Ev2a(:,m)=Ev2a(:,m)+Ev(n,m)*Ev2(:,n);
    end
end
 
% ������ѡȡһ���ض��Ĳ���������������һ��΢��

M=10; % ѡȡ��M������������΢�ż���
psi0=Ev2a(:,M);
psi=psi0;
% ���΢����
A=[];
for m=1:L
   Vkl=Ev1(:,m)'*V*psi;
   a=Vkl/(0-E1(m));
   A=[A,a];
   psi=psi+a*Ev1(:,m); 
end

psi=psi/sqrt(sum(psi.*conj(psi)));
plot(psi.*conj(psi))
hold on;
% ��һ����ED�Ľ���Ա�
% 
% H=H1(L,t0,lambda,J0);
% [Evm,E]=eig(H,'vector');
% psi=Evm(:,M);
% plot(psi.*conj(psi))
% -log(Ipr(psi))/log(2*L)

plot(E,psi0.*conj(psi0))
% figure()
% plot(A)


end

function []=PerturbationAnalytic4()
% ��֤��΢���۹�ʽ

%parameters
% m,n Ϊѭ��ָ�꣬��ʵ������
% a,b Ϊѭ��ָ�꣬��ʵ������

% ׼������
L=1000;
lambda=0.1;
J0=2.5;
t0=1;
omega=(sqrt(5)-1)/2;


% EigenStates and EigenValues

% localized eigenstates
% E !=0
E1=J0*ones(1,L)+2*t0*cos(2*pi*omega*linspace(1,L,L));
Ev1=kron(eye(L),[1;1])/sqrt(2);
% E=0 EigenStates
Ev2=kron(eye(L),[1;-1])/sqrt(2);

% ͨ�����ڷ��̼�������һ��΢�źͲ��������΢��
% ��д�� ΢�ž���
h2=kron(lambda*ones(1,L-1),[0,1]);
h2=[h2,0];
V=diag(h2,-1)+diag(h2,1);
% V(2*L,1)=lambda;
% V(1,2*L)=lambda;

% ���ڷ��̾���
V1=zeros(L,L);
for m=1:L
    for n=1:L
       psi1=Ev2(:,n);
       psi2=Ev2(:,m)';
       V1(m,n)=psi2*V*psi1;       
    end
end
%PBC case �����΢�ź󣬴��ڼ�
%OBC case ����������free hopping

% �����ڷ��̵õ���������׽���
[Ev,E]=eig(V1,'vector');
% �õ��㼶������������һ����������
Ev2a=zeros(2*L,L);
for m=1:L
    for n=1:L
        Ev2a(:,m)=Ev2a(:,m)+Ev(n,m)*Ev2(:,n);
    end
end
 
% ������ѡȡһ���ض��Ĳ���������������һ��΢��
Data=zeros(1,L);
for k=1:L
    M=k; % ѡȡ��M������������΢�ż���
    psi0=Ev2a(:,M);
    psi=psi0;
    % ���΢����

    for m=1:L
       Vkl=Ev1(:,m)'*V*psi;
       a=Vkl/(0-E1(m));
       psi=psi+a*Ev1(:,m); 
    end

    psi=psi/sqrt(sum(psi.*conj(psi)));
    Data(1,k)=-log(Ipr(psi))/log(2*L)
end

plot(E,Data,'*')
% axis([min(E),max(E),0,1])

hold on;
H=H1(L,t0,lambda,J0);
[Evm,Ea]=eig(H,'vector');
Data=zeros(1,2*L);
for M=1:2*L
    psi=Evm(:,M);
    Data(1,M)=-log(Ipr(psi))/log(2*L);
end
plot(Ea,Data,'.')

ylim([0,1])

figure()
plot(real(E),imag(E),'.')
hold on;
plot(real(Ea),imag(Ea),'o')

end

function a2 = Ipr(psi)
%IPR get the Ipr for a vector \sum_i|psi_i|^4
%  
a=psi.*conj(psi);
a2=sum(a.^2);
end