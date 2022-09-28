function EqAdj2()
    % arXiv:2208.13221
    % Singh J.
    % Shooting method for solving two-point boundary value problems 
    % in ODEs numerically
    %
    [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas();
    [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter); 
    [N,cols]=size(X);
    disp('Initial values:')
    disp(X(1,:))
    disp('Final values:')
    disp(X(N,:))
    clf
    plot(t,X(:,1),'r',t,X(:,4),'g--')
    legend('f(t)','theta(t)');
end

function [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas()
    n=5;
    m=3;
    t0=0;
    tf=5;
    y01=0;
    y02=1;
    y04=1;
    y0=[y01;y02;y04];
    yf=[0;0];
    % A tolerance
    eps=1.0e-6;
    % Maximum number of iterations
    maxIter=50;
    % Initial approximation
    c=[-2;0];
    % Number of nodes on [t0,tf]
    N=101; %21
end

function y=f(t,x)
    % Differential system of the BVP
    % x'(t)=f(t,x(t))
    k=0.71;
    y=zeros(5,1);
    y(1)=x(2);
    y(2)=x(3);
    y(3)=x(2)^2-x(1)*x(3);
    y(4)=x(5);
    y(5)=-k*x(1)*x(5);
end

function [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter)
    y01=y0(1);y02=y0(2);y04=y0(3);
    %x0=zeros(5,1);
    %x0(1)=y01;
    %x0(2)=y02;
    %x0(4)=y04;
    iter=0;
    sw=true;
    z0=zeros(10,1);
    z0(1)=y01;
    z0(2)=y02;
    z0(4)=y04;
    while sw 
        iter=iter+1;
        %x0(3)=c(1);
        %x0(5)=c(2);
        %[t,Z]=ode45(@f,[t0,tf],x0);
        %[N,cols]=size(Z);
        %F=yf-[Z(N,2);Z(N,4)];   
        z0(3)=c(1);
        z0(5)=c(2);        
        z0(8)=1;
        z0(10)=0;
        [t,Z1]=ode45(@g,[t0,tf],z0);
        [N,cols]=size(Z1);
        F=yf-[Z1(N,2);Z1(N,4)];   
        a11=Z1(N,7);
        a21=Z1(N,9);
        z0(8)=0;
        z0(10)=1;
        [t,Z2]=ode45(@g,[t0,tf],z0);
        [N,cols]=size(Z2);
        a12=Z2(N,7);
        a22=Z2(N,9);
        D=[a11,a12;a21,a22];  
        dc=inv(D)*F;
        nrm=norm(dc,'inf');
        c=c+dc;
        fprintf('iter: %d c(1): %f c(2): %f\n',iter,c(1),c(2))
        if nrm<eps || iter>=maxIter 
            sw=false;
        end
        if nrm<eps 
            ind=0;
        else
            ind=1;
        end
    end
    x0=[y01;y02;c(1);y04;c(2)];
    [t,X]=ode45(@f,[t0,tf],x0);
end

function y=g(t,x)
    k=0.71;
    y=zeros(10,1);
    y(1)=x(2);
    y(2)=x(3);
    y(3)=x(2)^2-x(1)*x(3);
    y(4)=x(5);
    y(5)=-k*x(1)*x(5);
    y(6)=x(7);
    y(7)=x(8);
    y(8)=-x(3)*x(6)+2*x(2)*x(7)-x(1)*x(8);
    y(9)=x(10);
    y(10)=-k*x(5)*x(6)-k*x(1)*x(10);
end

