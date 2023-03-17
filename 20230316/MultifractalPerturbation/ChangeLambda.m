%  energy spectrum when tuning the lambda


% ChangeLambda1()

% TypicalMfStates()
% perturbation()
zeropoints()

function []=TypicalMfStates()

    L=300;
    t0=1;
    
    lambda=0.05;
    H=H1(L,t0,lambda);
    [Ev,E]=eig(H,'vector');    
    
    psi=Ev(:,L);
    
%     yyaxis left
    plot(psi.*conj(psi))
    hold on;
    
    %plot the hopping term
%     yyaxis right
%     x=[1,kron(linspace(2,L-1,L-2),[1,1]),L];
%     length(x)
%     y0=diag(H,1);
%     y=kron(y0,[1,1]);
%       
%     plot(2*x,y(1:2*L-2),'.')
%     ylim([-0.1,0.1])
end

function []=ChangeLambda1()

    L=300;
    t0=1;
    
    lambda_all=0:0.5:2;

    figure()    
    for m=1:length(lambda_all)
        lambda=lambda_all(m);
        H=H1(L,t0,lambda);
        [Ev,E]=eig(H,'vector');    
        plot(E,'.')
        hold on;
    end
    
end

function H=H1(L,t0,lambda)

    omega=(sqrt(5)-1)/2;
    h1=kron(2*t0*cos(2*pi*omega*linspace(1,L-1,L-1)),[1,0]);
    h2=kron(lambda*ones(1,L-1),[0,1]);
    h=h1+h2;
    h=[h,2*t0*cos(2*pi*omega*L)];
    V=kron(2*t0*cos(2*pi*omega*linspace(1,L,L)),[1,1]);
    H=diag(h,-1)+diag(h,1)+diag(V);
   
    
end

function H=H2(L,t0,lambda,a)

    omega=(sqrt(5)-1)/2;
    h1=kron(a+2*t0*cos(2*pi*omega*linspace(1,L-1,L-1)),[1,0]);
    h2=kron(lambda*ones(1,L-1),[0,1]);
    h=h1+h2;
    h=[h,a+2*t0*cos(2*pi*omega*L)];
    V=kron(2*t0*cos(2*pi*omega*linspace(1,L,L)),[1,1]);
    H=diag(h,-1)+diag(h,1)+diag(V);
 
end

function H=H3(L,t0,lambda,J0)

% J0+cos term
    omega=(sqrt(5)-1)/2;
    h1=kron(J0*ones(1,L-1)+2*t0*cos(2*pi*omega*linspace(1,L-1,L-1)),[1,0]);
    h2=kron(lambda*ones(1,L-1),[0,1]);
    h=h1+h2;
    h=[h,J0+2*t0*cos(2*pi*omega*L)];
    V=kron(J0+2*t0*cos(2*pi*omega*linspace(1,L,L)),[1,1]);
    H=diag(h,-1)+diag(h,1)+diag(V);
   
    
end

function []=perturbation()
% we consider the effect of the perturbation
    L=20;
    t0=1;
    lambda_all=0:0.01:2;
    
    
    data1=zeros(2*L,length(lambda_all));
    data2=zeros(2*L,length(lambda_all));
        
    for m=1:length(lambda_all)
        lambda=lambda_all(m);
        a=0.1*lambda;
        H=H2(L,t0,lambda,a);
        [Ev,E]=eig(H,'vector');    
        
        data1(:,m)=E';
        for k=1:2*L
            psi=Ev(:,k);
            data2(k,m)=-log(Ipr(psi))/(log(2*L));
        end
    end
    
    L1=2*L;
    L2=length(lambda_all);
    Data=zeros(L1,L2);

    for i=1:L2
       Data(:,i)=ones(L1,1)*0.1*(i-1);  
    end
    figure()
    for i=1:L2
        scatter(Data(:,i),data1(:,i),'.','cdata',(data2(:,i)))
        hold on;
    end      
    
end

function []=zeropoints()
% we consider the effect of the perturbation
    L=200;
    t0=1;
    lambda_all=0:0.01:2;
    J=2.1;
    
    data1=zeros(2*L,length(lambda_all));
    data2=zeros(2*L,length(lambda_all));
        
    for m=1:length(lambda_all)
        lambda=lambda_all(m);
        a=J;
        H=H3(L,t0,lambda,a);
        [Ev,E]=eig(H,'vector');    
        
        data1(:,m)=E';
        for k=1:2*L
            psi=Ev(:,k);
            data2(k,m)=-log(Ipr(psi))/(log(2*L));
        end
    end
    
    L1=2*L;
    L2=length(lambda_all);
    Data=zeros(L1,L2);

    for i=1:L2
       Data(:,i)=ones(L1,1)*0.1*(i-1);  
    end
    figure()
    for i=1:L2
        scatter(Data(:,i),data1(:,i),'.','cdata',(data2(:,i)))
        hold on;
    end      
    
    colorbar()
    H
    
end

function a2 = Ipr(psi)
%IPR get the Ipr for a vector \sum_i|psi_i|^4
%  
a=psi.*conj(psi);
a2=sum(a.^2);
end