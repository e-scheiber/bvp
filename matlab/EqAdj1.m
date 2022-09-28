function EqAdj1()
    % arXiv:2208.13221
    % Singh J.
    % Shooting method for solving two-point boundary value problems 
    % in ODEs numerically
    %
    [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas();
    [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter);
    [N,cols]=size(X);
    Y=sol(t,N);
    nrm=norm(X-Y','inf');
    fprintf('||X-sol||_inf: %g\n',nrm);
    disp('Initial values:')
    disp(X(1,:))
    disp('Final values:');
    disp(X(N,:))
    clf
    plot(t,X(:,1),'r',t,X(:,2),'g--')
    legend('x(t)','dx(t)/df')
end

function [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas();
    n=2;
    m=1;
    t0=0;
    tf=1;
    y0=0;
    yf=2;
    % A tolerance
    eps=1.0e-6;
    % Maximum number of iterations
    maxIter=100;
    % Initial approximation
    c=2;
    % Number of nodes on [t0,tf]
    N=21;
end

function y=f(t,x)
    % Differential system of the BVP
    % x'(t)=f(t,x(t))
    y=zeros(2,1);
    y(1)=x(2);
    y(2)=2*x(1)*x(2);
end

function y=g(t,x)
    y=zeros(4,1);
    y(1)=x(2);
    y(2)=2*x(1)*x(2);
    y(3)=x(4);
    y(4)=2*x(2)*x(3)+2*x(1)*x(4);
end

function [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter)
    x0=zeros(4,1);
    x0(1)=y0;
    x0(4)=1;
    iter=0;
    sw=true;
    while sw 
        iter=iter+1;
        x0(2)=c;
        [t,Z]=ode45(@g,[t0,tf],x0);
        [N,cols]=size(Z);
        dc=(yf-Z(N,1))/Z(N,3);
        nrm=norm(dc,'inf');
        c=c+dc;
        fprintf('iter: %d c: %f\n',iter,c)
        if nrm<eps || iter>=maxIter 
            sw=false;
        end
        if nrm<eps 
            ind=0;
        else
            ind=1;
        end
    end
    x0=[y0;c];
    [t,X]=ode45(@f,[t0,tf],x0);
end

function Y=sol(t,N)
    a=1.0768740;
    %a=3.6435972
    s=@(t)[a*tan(a*t),a^2./cos(a*t)^2];
    Y=zeros(2,N);
    for i=1:N
      Y(:,i)=s(t(i));
    end
end