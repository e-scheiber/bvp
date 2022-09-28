using ODE
using LinearAlgebra
using DelimitedFiles

function EqAdj4()
    # arXiv:2208.13221
    # Singh J.
    # Shooting method for solving two-point boundary value problems 
    # in ODEs numerically
    #
    n=2
    m=1
    t0=1.0
    tf=2.0
    y0=2.0
    yf=2.5
    # A tolerance
    eps=1.0e-6
    # Maximum number of iterations
    maxIter=50
    # Initial approximation
    c=-0.3
    # Number of nodes on [t0,tf]
    N=21

    
    function f(t,x)
        # Differential system of the BVP
        # y'(x)=f(t,x)
        y::Array{Float64,1}=zeros(2)
        y[1]=x[2]
        y[2]=2*x[1]^3-6*x[1]-2*t^3
        return y
    end

    function g(t,x)
        y::Array{Float64,1}=zeros(4)
        y[1]=x[2]
        y[2]=2*x[1]^3-6*x[1]-2*t^3
        y[3]=x[4]
        y[4]=6*(x[1]^2-1)*x[3]
        return y
    end

    function solver(f,g,n,m,t0,tf,y0,yf,c,N,eps,maxIter)
        x0::Array{Float64,1}=zeros(4)
        x0[1]=y0
        x0[4]=1.0
        t=LinRange(t0,tf,N)
        iter=0
        sw=true
        while sw 
            iter=iter+1
            x0[2]=c
            #Z=ode('rk',x0,t0,t,g)
            (tout,Z)=ode45(g,x0,t)
            #println(size(Z))
            N=length(tout)
            #[tout[:],Z[:]]
            dc=(yf-Z[N][1])/Z[N][3]
            nrm=norm(dc,Inf)
            c=c+dc    
            println("iter: $iter c: $c")
            if nrm<eps || iter>=maxIter 
                sw=false
            end
            if nrm<eps 
                ind=0
            else
                ind=1
            end
        end
        x0=[y0,c]
        (tout,X)=ode45(f,x0,t)
        return tout,X
    end

    l=solver(f,g,n,m,t0,tf,y0,yf,c,N,eps,maxIter)
    t=l[1]
    X=l[2]
    
    N=length(t)
    println("Initial values:")
    display(X[1][:])
    println("Final values:")
    display(X[N][:])
    
    # Graphic through gnuplot
    y1::Array{Float64,1}=zeros(N)
    y2::Array{Float64,1}=zeros(N)
    for i=1:N
        y1[i]=X[i][1]
        y2[i]=X[i][2]
    end
      
    # Please set the folder where to save the results
    path="e:/mk/Julia/bvp1"
    writedlm(path*"/eq1.csv",[t y1 y2],", ")
end

