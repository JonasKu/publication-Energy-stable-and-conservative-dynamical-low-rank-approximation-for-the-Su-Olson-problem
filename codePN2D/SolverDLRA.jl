__precompile__

using ProgressMeter
using LinearAlgebra
using LegendrePolynomials
using QuadGK
using SparseArrays
using SphericalHarmonicExpansions,SphericalHarmonics,TypedPolynomials,GSL
using MultivariatePolynomials
using Einsum
using PyCall

include("PNSystem.jl")

struct SolverDLRA
    # spatial grid of cell interfaces
    x::Array{Float64};
    y::Array{Float64};

    # Solver settings
    settings::Settings;

    # preallocate memory for performance
    outRhs::Array{Float64,3};
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # Roe matrix
    AbsAx::Array{Float64,2};
    AbsAz::Array{Float64,2};
    # normalized Legendre Polynomials
    P::Array{Float64,2};
    # quadrature points
    mu::Array{Float64,1};
    w::Array{Float64,1};

    # functionalities of the PN system
    pn::PNSystem;

    L1x::SparseMatrixCSC{Float64, Int64};
    L1y::SparseMatrixCSC{Float64, Int64};
    L2x::SparseMatrixCSC{Float64, Int64};
    L2y::SparseMatrixCSC{Float64, Int64};

    ghostIdx::Vector{Int64};
    boundaryIdx::Vector{Int64};

    # constructor
    function SolverDLRA(settings)
        x = settings.x;
        y = settings.y;

        # setup flux matrix
        gamma = zeros(settings.nPN+1);
        for i = 1:settings.nPN+1
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end

        # construct PN system matrices
        pn = PNSystem(settings);
        SetupSystemMatrices(pn);

        outRhs = zeros(settings.NCellsX,settings.NCellsY,pn.nTotalEntries);

        # setup Roe matrix
        S = eigvals(pn.Ax);
        V = eigvecs(pn.Ax);
        AbsAx = V*abs.(diagm(S))*inv(V);

        S = eigvals(pn.Az);
        V = eigvecs(pn.Az);
        AbsAz = V*abs.(diagm(S))*inv(V);

        # compute normalized Legendre Polynomials
        Nq=200;
        (mu,w) = gauss(Nq);
        P=zeros(Nq,settings.nPN);
        for k=1:Nq
            PCurrent = collectPl(mu[k],lmax=settings.nPN-1);
            for i = 1:settings.nPN
                P[k,i] = PCurrent[i-1]/sqrt(gamma[i]);
            end
        end

        # setupt stencil matrix
        nx = settings.NCellsX;
        ny = settings.NCellsY;
        N = pn.nTotalEntries;
        L1x = spzeros(nx*ny,nx*ny);
        L1y = spzeros(nx*ny,nx*ny);
        L2x = spzeros(nx*ny,nx*ny);
        L2y = spzeros(nx*ny,nx*ny);

        # setup index arrays and values for allocation of stencil matrices
        II = zeros(3*(nx-2)*(ny-2)); J = zeros(3*(nx-2)*(ny-2)); vals = zeros(3*(nx-2)*(ny-2));
        counter = -2;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 3;
                # x part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i+1,j);
                indexMinus = vectorIndex(nx,i-1,j);

                II[counter+1] = index;
                J[counter+1] = index;
                vals[counter+1] = 2.0/2/settings.dx; 
                if i > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dx;
                end
                if i < nx
                    II[counter+2] = index;
                    J[counter+2] = indexPlus;
                    vals[counter+2] = -1/2/settings.dx; 
                end
            end
        end
        L1x = sparse(II,J,vals,nx*ny,nx*ny);

        II .= zeros(3*(nx-2)*(ny-2)); J .= zeros(3*(nx-2)*(ny-2)); vals .= zeros(3*(nx-2)*(ny-2));
        counter = -2;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 3;
                # y part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i,j+1);
                indexMinus = vectorIndex(nx,i,j-1);

                II[counter+1] = index;
                J[counter+1] = index;
                vals[counter+1] = 2.0/2/settings.dy; 

                if j > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dy;
                end
                if j < ny
                    II[counter+2] = index;
                    J[counter+2] = indexPlus;
                    vals[counter+2] = -1/2/settings.dy; 
                end
            end
        end
        L1y = sparse(II,J,vals,nx*ny,nx*ny);

        II = zeros(2*(nx-2)*(ny-2)); J = zeros(2*(nx-2)*(ny-2)); vals = zeros(2*(nx-2)*(ny-2));
        counter = -1;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 2;
                # x part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i+1,j);
                indexMinus = vectorIndex(nx,i-1,j);

                if i > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dx;
                end
                if i < nx
                    II[counter+1] = index;
                    J[counter+1] = indexPlus;
                    vals[counter+1] = 1/2/settings.dx;
                end
            end
        end
        L2x = sparse(II,J,vals,nx*ny,nx*ny);

        II .= zeros(2*(nx-2)*(ny-2)); J .= zeros(2*(nx-2)*(ny-2)); vals .= zeros(2*(nx-2)*(ny-2));
        counter = -1;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 2;
                # y part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i,j+1);
                indexMinus = vectorIndex(nx,i,j-1);

                if j > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dy;
                end
                if j < ny
                    II[counter+1] = index;
                    J[counter+1] = indexPlus;
                    vals[counter+1] = 1/2/settings.dy;
                end
            end
        end
        L2y = sparse(II,J,vals,nx*ny,nx*ny);

        # determine boundary indices and ghost cells
        BCType = "periodic"
        boundaryIdx = [];
        ghostIdx = [];
        for i = 2:nx-1
            j = 1;
            ghostIdx = [ghostIdx; vectorIndex(nx,i,j)]
            if BCType != "periodic"
                j = 2;
            else
                j = ny-1;
            end
            boundaryIdx = [boundaryIdx; vectorIndex(nx,i,j)]
            j = ny;
            ghostIdx = [ghostIdx; vectorIndex(nx,i,j)]
            if BCType != "periodic"
                j = ny-1;
            else
                j = 2;
            end
            boundaryIdx = [boundaryIdx; vectorIndex(nx,i,j)]
        end

        for j = 2:ny-1
            i = 1;
            ghostIdx = [ghostIdx; vectorIndex(nx,i,j)]
            if BCType != "periodic"
                i = 2;
            else
                i = nx-1;
            end
            boundaryIdx = [boundaryIdx; vectorIndex(nx,i,j)]
            i = nx;
            ghostIdx = [ghostIdx; vectorIndex(nx,i,j)]
            if BCType != "periodic"
                i = nx-1;
            else
                i = 2;
            end
            boundaryIdx = [boundaryIdx; vectorIndex(nx,i,j)]
        end

        new(x,y,settings,outRhs,gamma,AbsAx,AbsAz,P,mu,w,pn,L1x,L1y,L2x,L2y,ghostIdx,boundaryIdx);
    end
end

py"""
import numpy
def qr(A):
    return numpy.linalg.qr(A)
"""

function SetupIC(obj::SolverDLRA)
    u = zeros(obj.settings.NCellsX,obj.settings.NCellsY,obj.pn.nTotalEntries);
    u[:,:,1] = IC(obj.settings,obj.settings.xMid,obj.settings.yMid)
    return u;
end

function Solve(obj::SolverDLRA)
    # Get rank
    s = obj.settings;
    # Set up initial condition and store as matrix
    nx = obj.settings.NCellsX;
    ny = obj.settings.NCellsY;

    # Set up initial condition
    u = zeros(nx * ny, obj.pn.nTotalEntries);
    uMat = IC(obj.settings,obj.settings.xMid,obj.settings.yMid);
    for i = 1:obj.pn.nTotalEntries
        u[:,i] = vec(uMat[:,:,i])
    end
    B = obj.settings.B0 * ones(nx*ny);

    nT = Int(ceil(s.tEnd/s.dt));
    dt = s.dt;

    # source term 
    Q = zeros(nx*ny)
    if obj.settings.problem == "SuOlsonTestcase"
        for j = 1:(nx*ny)
            if obj.settings.xMid[j] >= -0.5 && obj.settings.xMid[j] <= 0.5
                Q[j] = 1.0/obj.settings.aRad/sig
            end
        end
    end

    # setup coupling matrix
    sig = obj.settings.sigmaA
    C = zeros(2,2);
    C[1,1] = 1+ sig*dt; C[2,2] = 1 + sig*dt;
    C[1,2] = sig*dt; C[2,1] = sig*dt;
    C ./= (1 + 2*sig*dt);

    prog = Progress(nT,1)
    t = 0.0;

    mass = zeros(2,nT)

    for n=1:nT

        # impose boundary conditions
        u[obj.ghostIdx,:] .= u[obj.boundaryIdx,:];
        B[obj.ghostIdx] .= B[obj.boundaryIdx];

        # compute mass evolution
        mass[1,n] = t;
        mass[2,n] = sum(u[:,1] .+ B) * obj.settings.dx * obj.settings.dy;
        mass[2,n] -= sum(u[obj.ghostIdx] .+ B[obj.ghostIdx])* obj.settings.dx * obj.settings.dy;
        
        u .=  u .- dt*(obj.L2x*u*obj.pn.Ax .+ obj.L2y*u*obj.pn.Az .+ obj.L1x*u*obj.AbsAx .+ obj.L1y*u*obj.AbsAz);        

        u[:,2:end] .= u[:,2:end]./(1 .+ dt .*sig);

        for j = 1:(nx*ny)
            uB = C*[u[j,1];B[j]];
            u[j,1] = uB[1];
            B[j] = uB[2];
        end
        t += dt;

        next!(prog) # update progress bar
    end

    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*u,B,mass;

end

function SolveDLRA(obj::SolverDLRA)

    # store needed variables from settings class
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.pn.nTotalEntries;
    nx = obj.settings.NCellsX;
    ny = obj.settings.NCellsY;

    r = obj.settings.r; # rank

    # Set up initial condition
    u = zeros(nx * ny, obj.pn.nTotalEntries);
    uMat = IC(obj.settings,obj.settings.xMid,obj.settings.yMid);
    for i = 1:obj.pn.nTotalEntries
        u[:,i] = vec(uMat[:,:,i])
    end
    B = obj.settings.B0 * ones(nx*ny);
   
    # allocate intermediate memory
    BNew = zeros(nx*ny);

    # Precompute flux and opacity matrices
    sig = obj.settings.sigmaA;

    e1 = [1.0;zeros(N-1)];

    # source term 
    Q = zeros(nx*ny)

    # initialization to rank r for DLRA
    X,S,V = svd(u);
    X = X[:,1:r];
    S = diagm(S[1:r]);
    V = V[:,1:r];
    
    # setup coupling matrix
    C = zeros(2,2);
    C[1,1] = 1 + sig*dt; C[2,2] = 1 + sig*dt;
    C[1,2] = sig*dt; C[2,1] = sig*dt;
    C ./= (1 + 2*sig*dt);

    rankInTime = zeros(2,nt);
    massInTime = zeros(2,nt);
    c = 1 # speed of particles

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;

        X[obj.ghostIdx,:] .= X[obj.boundaryIdx,:];
        B[obj.ghostIdx] .= B[obj.boundaryIdx];

        S_tmp = S;
        
        # K step

        K = X * S;
        VAxV = (V'*obj.pn.Ax*V)';
        VAzV = (V'*obj.pn.Az*V)';
        VAbsAxV = (V'*obj.AbsAx*V)';
        VAbsAzV = (V'*obj.AbsAx*V)';

        K .= K .- dt .* obj.L2x*K*VAxV .- dt .* obj.L2y*K*VAzV .-  dt .* obj.L1x*K*VAbsAxV .-  dt .* obj.L1y*K*VAbsAzV .+ dt*Q*(e1'* V);

        # augment basis to 2r

        K = [K X];
        XNew,_ = py"qr"(K);
        MUp = XNew'*X;

        # L step
    
        L = V * S';
        XL2xX = (X'*obj.L2x*X)';
        XL2yX = (X'*obj.L2y*X)';
        XL1xX = (X' * obj.L1x *X)';
        XL1yX = (X' * obj.L1y *X)';

        L .= L .- dt .* obj.pn.Ax*L*XL2xX .- dt .* obj.pn.Az*L*XL2yX .- dt.* obj.AbsAx*L*XL1xX .- dt .* obj.AbsAz*L*XL1yX .+ dt*e1*(Q'*X);

        # augment basis to 2r
        L = [L V];
        VNew,_= py"qr"(L);
        NUp = VNew'*V;

        # S step

        S = MUp*S*NUp';

        XNewL2xXNew = XNew'*obj.L2x* XNew;
        XNewL2yXNew = XNew'*obj.L2y* XNew;
        XNewL1xXNew = XNew'*obj.L1x*XNew;
        XNewL1yXNew = XNew'*obj.L1y*XNew;

        VNewAxVNew = (VNew' *obj.pn.Ax*VNew)';
        VNewAzVNew = (VNew' *obj.pn.Az*VNew)';
        VNewAbsAxVNew = (VNew'*obj.AbsAx*VNew)';
        VNewAbsAzVNew = (VNew'*obj.AbsAz*VNew)';

        S = S .- dt .* XNewL2xXNew*S* VNewAxVNew .- dt .* XNewL2yXNew*S* VNewAzVNew .- dt .* XNewL1xXNew*S*VNewAbsAxVNew .- dt .* XNewL1yXNew*S*VNewAbsAzVNew .+ dt* (XNew' *Q)*(e1' * VNew);
        
        u0_tmp = X * (S_tmp * V[1,:]) .- dt .* obj.L2x*X*S_tmp*(V'* obj.pn.Ax[:,1]) .- dt .* obj.L2y*X*S_tmp*(V'* obj.pn.Az[:,1]) .- dt .* obj.L1x*X*S_tmp*(V' * obj.AbsAx[:,1]) .- dt .* obj.L1y*X*S_tmp*(V' * obj.AbsAz[:,1]);
        
        u0 = zeros(nx*ny)
        for j = 1:(nx*ny)
            uB = C*[u0_tmp[j];B[j]];
            u0[j] = uB[1];
            BNew[j] = uB[2];
        end
       
        # scattering step
        L = VNew * S';
        for i = 1:size(L,2)
            L[2:end,i] = L[2:end,i]./(1 .+ dt .*sig);
        end
        VNew,ST = py"qr"(L); S = ST';

        # append u0 to basis
        X,_ = py"qr"([u0 XNew]); 
        V,_ = py"qr"([e1 VNew]); # can remove e1

 
        S = (X'*XNew)*S*(VNew'*(I - e1*e1')*V) + X'*u0*e1'*V;

        massInTime[1,n] = t;
        massInTime[2,n] = sum(X*S*V[1,:] / c .+ BNew)* obj.settings.dx * obj.settings.dy;
        massInTime[2,n] -= sum(u0[obj.ghostIdx] / c .+ BNew[obj.ghostIdx])* obj.settings.dx * obj.settings.dy;

        ################## truncate ##################
        X,S,V = truncate(obj,X,S,V)

        # update rank
        r = size(S,1);

        B .= BNew;

        # println(minimum(B));

        t += dt;

        next!(prog) # update progress bar
    end
    u .= X*S*V';
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*u, B, rankInTime, massInTime;
end

function vectorIndex(nx,i,j)
    return (i-1)*nx + j;
end

function Vec2Mat(nx,ny,v)
    m = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            m[i,j] = v[(i-1)*nx + j]
        end
    end
    return m;
end

function truncate(obj,X,S,W)
    # println("---------------------------------")
    r0 = size(S,1);
    rMaxTotal = obj.settings.rMaxTotal;
    e1 = [1.0;zeros(obj.pn.nTotalEntries-1)];
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
    # println("norm sigma = ",norm(Sigma))
    # println("tol = ",tol)

    rmax = Int(floor(size(Sigma,1)/2));

    for j=1:2*rmax
        tmp = sqrt(sum(Sigma[j:2*rmax]).^2);
        # println("-> tmp = ",tmp, " at j = ",j)
        if tmp < tol
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