__precompile__

using ProgressMeter
using LinearAlgebra
using LegendrePolynomials
using QuadGK
using PyCall

struct Solver
    # spatial grid of cell interfaces
    x::Array{Float64};

    # Solver settings
    settings::Settings;

    # preallocate memory for performance
    outRhs::Array{Float64,2};
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # flux matrix PN system
    A::Array{Float64,2};
    # Roe matrix
    AbsA::Array{Float64,2};

    # stencil matrices
    Dx::Tridiagonal{Float64, Vector{Float64}};
    Dxx::Tridiagonal{Float64, Vector{Float64}};

    # physical parameters
    sigmaA::Float64;

    # constructor
    function Solver(settings)
        x = settings.x;
        nx = settings.NCells;
        dx = settings.dx;

        outRhs = zeros(nx,settings.nPN);

        # setup flux matrix
        gamma = ones(settings.nPN);

        # setup gamma vector
        gamma = zeros(settings.nPN);
        for i = 1:settings.nPN
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        
        # setup flux matrix
        A = zeros(settings.nPN,settings.nPN);

        for i = 1:(settings.nPN-1)
            n = i-1;
            A[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
        end

        for i = 2:settings.nPN
            n = i-1;
            A[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
        end

        # setup Roe matrix
        S = eigvals(A)
        V = eigvecs(A)
        AbsA = V*abs.(diagm(S))*inv(V)

        # set up spatial stencil matrices
        Dx = Tridiagonal(-ones(nx-1)./dx/2.0,zeros(nx),ones(nx-1)./dx/2.0); # central difference matrix
        Dxx = Tridiagonal(ones(nx-1)./dx/2.0,-ones(nx)./dx,ones(nx-1)./dx/2.0); # stabilization matrix
        Dx[1,1] = 0.0; Dx[1,2] = 0.0; Dx[end,end] = 0.0; Dx[end,end-1] = 0.0;
        Dxx[1,1] = 0.0; Dxx[1,2] = 0.0; Dxx[end,end] = 0.0; Dxx[end,end-1] = 0.0;

        new(x,settings,outRhs,gamma,A,AbsA,Dx,Dxx,settings.sigmaA);
    end
end

py"""
import numpy
def qr(A):
    return numpy.linalg.qr(A)
"""

function Solve(obj::Solver)

    # store needed variables from settings class
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    dy = obj.settings.dy;
    tEnd = obj.settings.tEnd;
    epsilon = obj.settings.epsilon;
    epsInv = 1/obj.settings.epsilon;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.settings.nPN;
    nx = obj.settings.NCells;
    ny = obj.settings.NCellsY;

    # Set up initial condition
    u = IC(obj.settings,obj.settings.xMid);
    B = obj.settings.B0 * ones(nx);

    # allocate intermediate memory
    uNew = zeros(nx,N);
    uHalf = zeros(nx,N);
    BNew = zeros(nx);

    # Precompute flux and opacity matrices
    sig = obj.settings.sigmaA;
    A = obj.A;
    AbsA = obj.AbsA;

    e1 = [1.0;zeros(N-1)];

    # source term 
    Q = zeros(nx);
    if obj.settings.problem == "SuOlsonTestcase"
        for j = 1:nx
            if obj.settings.xMid[j] >= -0.5 && obj.settings.xMid[j] <= 0.5
                Q[j] = 1.0/obj.settings.aRad/sig;
            end
        end
    end

    # setup coupling matrix
    C = zeros(2,2);
    C[1,1] = epsilon^2 + sig*dt; C[2,2] = epsilon^2 + sig*dt;
    C[1,2] = sig*dt; C[2,1] = sig*dt;
    C ./= (epsilon^2 + 2*sig*dt);

    massInTime = zeros(2,nt);
    c = 1;  # speed of particles

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        uHalf .= u .- dt * epsInv .* obj.Dx*u*A' .+ dt  * epsInv .* obj.Dxx*u*AbsA' .+ dt*Q*e1'
        uNew[:,2:end] = uHalf[:,2:end]./(1 .+ dt* epsInv^2 .*sig);

        for j = 1:nx
            uB = C*[uHalf[j,1];B[j]];
            uNew[j,1] = uB[1];
            BNew[j] = uB[2];
        end
        u .= uNew;
        B .= BNew;

        massInTime[1,n] = t;
        massInTime[2,n] = sum(u[:,1]/c + B) * dx;        ;

        t += dt;
        next!(prog) # update progress bar
    end

    vec=zeros(ny);
    mat=zeros(N,ny);
    for i=1:ny
        vec=collectPl(-1+(i-1)*dy;lmax=N-1);
        mat[:,i] = vec;
    end
    
    f=zeros(nx,ny);
    f=u*mat;
    

    # return end time and solution
    return t, u, B, massInTime,f;

end

function SolveDLRA(obj::Solver)

    # store needed variables from settings class
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    dy = obj.settings.dy;
    tEnd = obj.settings.tEnd;
    epsilon = obj.settings.epsilon;
    epsInv = 1/obj.settings.epsilon;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.settings.nPN;
    nx = obj.settings.NCells;
    ny = obj.settings.NCellsY;

    r = obj.settings.r; # rank

    # Set up initial condition
    u = IC(obj.settings,obj.settings.xMid);
    B = obj.settings.B0 * ones(nx);
   
    # allocate intermediate memory
    BNew = zeros(nx);

    # Precompute flux and opacity matrices
    sig = obj.settings.sigmaA;
    A = obj.A;
    AbsA = obj.AbsA;

    e1 = [1.0;zeros(N-1)];

    # source term 
    Q = zeros(nx)
    if obj.settings.problem == "SuOlsonTestcase"
        for j = 1:nx
            if obj.settings.xMid[j] >= -0.5 && obj.settings.xMid[j] <= 0.5
                Q[j] = 1.0/obj.settings.aRad/sig;
            end
        end
        # initialization to rank r for DLRA using the source term
        X,S,V = svd(dt*Q*e1');
        X = X[:,1:r];
        S = zeros(r,r);
        V = V[:,1:r];
    else
        # initialization to rank r for DLRA using the source term
        X,S,V = svd(u);
        X = X[:,1:r];
        S = diagm(S[1:r]);
        V = V[:,1:r];
    end

    # setup coupling matrix
    C = zeros(2,2);
    C[1,1] = epsilon^2 + sig*dt; C[2,2] = epsilon^2 + sig*dt;
    C[1,2] = sig*dt; C[2,1] = sig*dt;
    C ./= (epsilon^2 + 2*sig*dt);

    rankInTime = zeros(2,nt);
    massInTime = zeros(2,nt);
    c = 1; # speed of particles

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;

        S_tmp = S;
        
        # K step

        K = X * S;
        VAV = V'*A'*V;
        VAbsAV = V'*AbsA'*V;

        K = K .- dt* epsInv .* obj.Dx*K*VAV .+  dt*epsInv .* obj.Dxx*K*VAbsAV .+ dt*Q*(e1'* V);

        # augment basis to 2r

        K = [K X];
        XNew,_ = qr!(K); XNew = Matrix(XNew); XNew = XNew[:,1:(2*r)];
        MUp = XNew'*X;

        # L step
    
        L = V * S';
        XDxX = X'*obj.Dx'*X;
        XDxxX = X' * obj.Dxx' *X;

        L = L .- dt * epsInv .* A*L*XDxX .+ dt*epsInv .* AbsA*L*XDxxX .+ dt*e1*(Q'*X);

        # augment basis to 2r

        L = [L V];
        VNew,_= qr!(L); VNew = Matrix(VNew); VNew = VNew[:,1:(2*r)];
        NUp = VNew'*V;

        # S step

        S = MUp*S*NUp';
        XNewDxXNew = XNew'*obj.Dx* XNew;
        XNewDxxXNew = XNew'*obj.Dxx*XNew;
        VNewAVNew = VNew' *A'*VNew;
        VNewAbsAVNew = VNew'*AbsA'*VNew;

        S = S .- dt  * epsInv .* XNewDxXNew*S* VNewAVNew .+ dt *epsInv .* XNewDxxXNew*S*VNewAbsAVNew .+ dt* (XNew' *Q)*(e1' * VNew);

        # update scalar flux and energy

        u0_tmp = X * (S_tmp * V[1,:]) .- dt * epsInv .* obj.Dx*X*S_tmp*(V'* A[:,1]) .+ dt*epsInv .* obj.Dxx*X*S_tmp*(V' * AbsA[:,1]) .+ dt*Q;

        u0 = zeros(nx)
        for j = 1:nx
            uB = C*[u0_tmp[j];B[j]];
            u0[j] = uB[1];
            BNew[j] = uB[2];
        end
       
        # scattering step
        L = VNew * S';
          for i = 1:size(L,2)
            L[2:end,i] = L[2:end,i]./(1 .+ dt* epsInv^2 .*sig);
          end
        VNew,ST = qr!(L); VNew = Matrix(VNew); VNew = VNew[:,1:(2*r)]; S = Matrix(ST)';

        # append u0 to basis
        X,_ = qr([u0 XNew]); X = Matrix(X); X = X[:,1:(2*r+1)];
        V,_ = qr([e1 VNew]); V = Matrix(V); V = V[:,1:(2*r+1)]; # can remove e1

        # correct coefficient matrix
        S = (X'*XNew)*S*(VNew'*(I - e1*e1')*V) + X'*u0*e1'*V

        ################## truncate ##################
        X,S,V = truncate(obj,X,S,V)

        # update rank
        r = size(S,1);

        B .= BNew;

        massInTime[1,n] = t;
        massInTime[2,n] = sum(X*S*V[1,:] / c + B) * dx;

        t += dt;
        next!(prog) # update progress bar
    end
    u .= X*S*V';

    vec=zeros(ny);
    mat=zeros(N,ny);
    for i=1:ny
        vec=collectPl(-1+(i-1)*dy;lmax=N-1);
        mat[:,i] = vec;
    end
    
    f=zeros(nx,ny);
    f=u*mat;

    # return end time and solution
    return t, u, B, rankInTime, massInTime, f;
end


function truncate(obj,X,S,W)
    r0 = size(S,1);
    rMaxTotal = obj.settings.rMaxTotal;
    e1 = [1.0;zeros(obj.settings.nPN-1)];
    K = X*S;
    Kcons = K[:,1]; Krem = K[:,2:end];
    Wcons = W[:,1]; Wrem = W[:,2:end];
    Xcons = Kcons ./ norm(Kcons); Scons = norm(Kcons);
    Xrem,Srem = py"qr"(Krem);
    U,Sigma,V = svd(Srem);
    # truncate
    rmax = -1;

    tmp = 0.0;
    tol = obj.settings.epsAdapt*norm(Sigma);

    rmax = Int(floor(size(Sigma,1)/2));

    for j=1:2*rmax
        tmp = sqrt(sum(Sigma[j:2*rmax]).^2);
        if(tmp<tol)
            rmax = j;
            break;
        end
    end

    rmax = min(rmax,rMaxTotal);
    r1 = max(rmax,2);

    # if 2*r was actually not enough move to highest possible rank
    if rmax == -1
        rmax = rMaxTotal;
    end

    Srem = Diagonal(Sigma[1:r1]);
    Xrem = Xrem * U[:,1:r1];
    Wrem = Wrem * V[:,1:r1];
    What = [e1 Wrem];
    Xhat = [Xcons Xrem];
    Xnew,R1 = py"qr"(Xhat);
    Wnew,R2 = py"qr"(What);
    Snew = R1*[Scons zeros(1,r1); zeros(r1,1) Srem]*R2';
    return Xnew, Snew, Wnew;
end