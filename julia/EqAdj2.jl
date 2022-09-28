using ODE
using LinearAlgebra
using DelimitedFiles
 
function EqAdj2()
    # arXiv:2208.13221
    # Singh J.
    # Shooting method for solving two-point boundary value problems 
    # in ODEs numerically
    #
    n=5
    m=3
    t0=0.0
    tf=5.0
    y01=0.0
    y02=1.0
    y04=1.0
    y0::Array{Float64,1}=[y01,y02,y04]
    yf::Array{Float64,1}=[0.0,0.0]
    # A tolerance
    eps=1.0e-6
    # Maximum number of iterations
    maxIter=50
    # Initial approximation
    c::Array{Float64,1}=[0.0,0.0]
    # Number of nodes on [t0,tf]
    N=101 #21
    
    function f(t,x::Array{Float64,1})
        # Differential system of the BVP
        # y'(x)=f(t,x)
        k=0.71
        y::Array{Float64,1}=zeros(5)
        y[1]=x[2]
        y[2]=x[3]
        y[3]=x[2]^2-x[1]*x[3]
        y[4]=x[5]
        y[5]=-k*x[1]*x[5]
        return y
    end

    function g(t,x::Array{Float64,1})
        k=0.71
        y::Array{Float64,1}=zeros(10)
        y[1]=x[2]
        y[2]=x[3]
        y[3]=x[2]^2-x[1]*x[3]
        y[4]=x[5]
        y[5]=-k*x[1]*x[5]
        y[6]=x[7]
        y[7]=x[8]
        y[8]=-x[3]*x[6]+2*x[2]*x[7]-x[1]*x[8]
        y[9]=x[10]
        y[10]=-k*x[5]*x[6]-k*x[1]*x[10]
        return y
    end

    function solver(f,g,n,m,t0,tf,y0,yf,c,N,eps,maxIter)
        y01=y0[1]
        y02=y0[2]
        y04=y0[3]
         t=LinRange(t0,tf,N)
        iter=0
        sw=true
        z0::Array{Float64,1}=zeros(10)
        z0[1]=y01
        z0[2]=y02
        z0[4]=y04
         while sw 
            iter=iter+1
            z0[3]=c[1]
            z0[5]=c[2]        
            z0[8]=1.0
            z0[10]=0.0
            (tout,Z1)=ode45(g,z0,t)
            N1=length(tout)
            F=yf-[Z1[N1][2],Z1[N1][4]]
            a11=Z1[N1][7]
            a21=Z1[N1][9]
            z0[8]=0.0
            z0[10]=1.0
            (tout,Z2)=ode45(g,z0,t)
            N2=length(tout)
            a12=Z2[N][7]
            a22=Z2[N][9] 
            D=[a11 a12;a21 a22]
            dc=inv(D)*F
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
        x0=[y01,y02,c[1],y04,c[2]]
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

