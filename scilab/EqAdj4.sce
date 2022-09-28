function EqAdj4()
    // arXiv:2208.13221
    // Singh J.
    // Shooting method for solving two-point boundary value problems 
    // in ODEs numerically
    //
    [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas()
    [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter) 
    disp('Initial values:')
    disp(X(:,1))
    disp('Final values:')
    disp(X(:,N))
    clf 
    plot(t',X(1,:)','r',t',X(2,:)','g--')
    legend(['$x_1(t)$','$x_2(t)$'])
endfunction

function [n,m,t0,tf,y0,yf,c,N,eps,maxIter]=datas()
    n=2
    m=1
    t0=1
    tf=2
    y0=2
    yf=2.5
    // A tolerance
    eps=1.0e-6
    // Maximum number of iterations
    maxIter=50
    // Initial approximation
    c=1
    // Number of nodes on [t0,tf]
    N=21
endfunction

function y=f(t,x)
    // Differential system of the BVP
    // x'(t)=f(t,x(t))
    y=zeros(2,1)
    y(1)=x(2)
    y(2)=2*x(1)^3-6*x(1)-2*t^3
endfunction

function y=g(t,x)
    y=zeros(4,1)
    y(1)=x(2)
    y(2)=2*x(1)^3-6*x(1)-2*t^3
    y(3)=x(4)
    y(4)=6*(x(1)^2-1)*x(3)
endfunction

function [t,X]=nonLShM(n,m,t0,tf,y0,yf,c,N,eps,maxIter)
    x0=zeros(4,1)
    x0(1)=y0
    x0(4)=1
    t=linspace(t0,tf,N)
    iter=0
    sw=%t
    while sw do
        iter=iter+1
        x0(2)=c
        Z=ode('rk',x0,t0,t,g)
        dc=(yf-Z(1,N))/Z(3,N)
        nrm=norm(dc,'inf')
        c=c+dc
        printf('iter: %d c: %f\n',iter,c)
        if nrm<eps || iter>=maxIter then
            sw=%f
        end
        if nrm<eps then
            ind=0
        else
            ind=1
        end
    end
    x0=[y0;c]
    X=ode('rk',x0,t0,t,f)
endfunction


