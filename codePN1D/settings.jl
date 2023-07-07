__precompile__
mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    Ny::Int64;
    # number spatial cells
    NCells::Int64;
    NCellsY::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    d::Float64;
    e::Float64;
    # grid cell width
    dx::Float64;
    dy::Float64;

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

    # physical parameters
    sigmaA::Float64;
    epsilon::Float64;

    # constants
    h::Float64;
    c::Float64;
    kB::Float64;
    aRad::Float64;

    # low rank parameters
    r::Int;
    rMaxTotal::Int;
    epsAdapt::Float64;
    B0::Float64;

    problem::String;

    function Settings(Nx::Int=1001)
        # spatial grid setting
        NCells = Nx - 1;
        a = -10.0; # left boundary
        b = 10.0;#0.05;#0.08; # right boundary
        d=-1.0;
        e=1.0;

        Ny=1001;
        
        NCellsY=Ny - 1;

        
        # time settings
        tEnd = 3.16;
        
        # number PN moments
        nPN = 500;
       # nPN = 501;

        y = collect(range(d,stop = e,length = NCellsY));
        dy = y[2]-y[1];
        y = [y[1]-dy;y]; # add ghost cells so that boundary cell centers lie on a and b
        y = y.+dy/2;
        yMid = y[1:(end-1)].+0.5*dy
    
        x = collect(range(a,stop = b,length = NCells));
        dx = x[2]-x[1];
        x = [x[1]-dx;x]; # add ghost cells so that boundary cell centers lie on a and b
        x = x.+dx/2;
        xMid = x[1:(end-1)].+0.5*dx

        # physical parameters
        epsilon = 1.0;

        # constants
        h = 6.62607015e-34*100^2*1000 # cm g / s
        
        kB = 1.38064852e-16; 

        # low-rank parameters
        r = 20;
        rMaxTotal = 100;

        # physical parameters
        sigmaA = 1.0
        c = 299792458.0 * 100.0;                # speed of light in [cm/s]
        sigmaSB = 5.6704 * 1e-5;                # Stefan Boltzmann constant in [erg/cm^2/s/K^4]
        aRad = 4.0 * sigmaSB / c;               # radiation constant [erg/(cm^3 K^4)]

        cfl = 0.99; # CFL condition

        dt = cfl*dx/epsilon^2;

        epsAdapt = 1e-2;

        # initial condition B
        B0 = 50;

        # test case settings
        problem = "SuOlsonTestcase" # SuOlsonPlaneSource, SuOlsonTestcase

        if problem == "SuOlsonPlaneSource"
            tEnd = 8;
            epsAdapt = 1e-1;
            B0 = 1;
        end

        # build class
        # new(Nx,NCells,a,b,dx,tEnd,dt,cfl,nPN,x,xMid,sigmaA,epsilon,h,c,kB,aRad,r,rMaxTotal,epsAdapt,B0,problem);
        new(Nx,Ny,NCells,NCellsY,a,b,d,e,dx,dy,tEnd,dt,cfl,nPN,x,xMid,y,yMid,sigmaA,epsilon,h,c,kB,aRad,r,rMaxTotal,epsAdapt,B0,problem);
    end

end

function IC(obj::Settings,x)
    y = zeros(length(x),obj.nPN);
    if obj.problem == "SuOlsonPlaneSource"
        x0 = 1.0
        s1 = 0.03
        s2 = s1^2
        floor = 1e-4
        for j = 1:length(x);
            y[j,1] = max(floor,1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2))
        end
    end
    return y;
end