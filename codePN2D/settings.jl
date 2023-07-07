__precompile__

include("quadratures/Quadrature.jl")

using SphericalHarmonicExpansions,SphericalHarmonics,TypedPolynomials,GSL
using MultivariatePolynomials

mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    Ny::Int64;
    # number spatial cells
    NCellsX::Int64;
    NCellsY::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    c::Float64;
    d::Float64;
    # grid cell width
    dx::Float64
    dy::Float64

    # time settings
    # end time
    tEnd::Float64;
    # time increment
    dt::Float64;
    # CFL number 
    cfl::Float64;
    
    # degree PN
    nPN::Int64;

    # spatial grid
    x
    xMid
    y
    yMid

    # problem definitions
    problem::String;

    # physical parameters
    sigmaA::Float64;

    # rank
    r::Int;
    rMaxTotal::Int;
    epsAdapt::Float64;

    B0::Float64;

    function Settings(Nx::Int=501,Ny::Int=501,r::Int=100)

        # test case settings
        problem = "Beam" # "Linesource", "Beam"

        # spatial grid setting
        NCellsX = Nx - 1;
        NCellsY = Ny - 1;

        if problem == "Linesource"
            a = -1.5; # left boundary
            b = 1.5; # right boundary

            c = -1.5; # lower boundary
            d = 1.5; # upper boundary

            B0 = 1.0;

            # physical parameters
            sigmaA = 0.1;      

            epsAdapt = 0.05;
        else
            a = -1.0; # left boundary
            b = 1.0; # right boundary

            c = -1.0; # lower boundary
            d = 1.0; # upper boundary

            B0 = 1;

            # physical parameters
            sigmaA = 0.5;     

            epsAdapt = 5*1e-4;
        end
         
        rMaxTotal = 100;

        # spatial grid
        x = collect(range(a,stop = b,length = NCellsX));
        dx = x[2]-x[1];
        x = [x[1]-dx;x]; # add ghost cells so that boundary cell centers lie on a and b
        x = x.+dx/2;
        xMid = x[1:(end-1)].+0.5*dx
        y = collect(range(c,stop = d,length = NCellsY));
        dy = y[2]-y[1];
        y = [y[1]-dy;y]; # add ghost cells so that boundary cell centers lie on a and b
        y = y.+dy/2;
        yMid = y[1:(end-1)].+0.5*dy

        # time settings
        tEnd = 0.5;
        # tEnd = 1.5;
        cfl = 0.7 # CFL condition
        dt = cfl*dx;
        
        # number PN moments
        nPN = 29 # use odd number

        # B0 = 50;

        # build class
        new(Nx,Ny,NCellsX,NCellsY,a,b,c,d,dx,dy,tEnd,dt,cfl,nPN,x,xMid,y,yMid,problem,sigmaA,r,rMaxTotal,epsAdapt,B0);
    end
end

function IC(obj::Settings,x,y)
    u = zeros(length(x),length(y),(obj.nPN+1)^2);
    if obj.problem == "Linesource"
        x0 = 0.0;
        y0 = 0.0;
        
        s1 = 0.01
        s2 = 0.03^2
        floor = 1e-14
        for j = 1:length(x);
            for i = 1:length(y);
                u[j,i,1] = max(floor,1e6/(4.0*pi*s2) *exp(-((x[j]-x0)*(x[j]-x0)+(y[i]-y0)*(y[i]-y0))/4.0/s2))/4.0/pi;
            end
        end
    elseif obj.problem == "Beam"
        # setup quadrature
        qorder = obj.nPN+1; # must be even for standard quadrature
        qtype = 1; # Type must be 1 for "standard" or 2 for "octa" and 3 for "ico".
        Q = Quadrature(qorder,qtype)
        # Construct Gauss quadrature
        mu,gaussweights = gausslegendre(qorder)

        weights = Q.weights
        Norder = (obj.nPN+1)^2
        nq = length(weights);
        nq = Q.nquadpoints;
        nx = obj.NCellsX;
        ny = obj.NCellsY;
        @polyvar xx yy zz
        phi_beam = pi/2;                               # Angle of beam w.r.t. x-axis.
        mu_beam = 0;  

        psi = zeros(nq);
        Omega1 = 1/sqrt(2); 
        Omega3 = 1/sqrt(2); 
        pointsxyz = Q.pointsxyz;

        for k = 1:nq 
            psi[k] = 1e6*normpdf(pointsxyz[k,1],Omega1,.1).*normpdf(pointsxyz[k,3],Omega3,.1);
        end

        O,M = ComputeTrafoMatrices(Q,Norder,obj.nPN)

        test = [];
        for i = 1:nx
            for j = 1:ny
                pos_beam = [0.0,0.0,0.0];
                space_beam = normpdf(obj.xMid[i],pos_beam[1],.1).*normpdf(obj.yMid[j],pos_beam[2],.1)
                u[i,j,:] = Float64.(M*psi)*space_beam
                #u = Float64.(M*psi)*space_beam
                #test = [test; u[1]];
            end
        end
    end
    
    return u;
end