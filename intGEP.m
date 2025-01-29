function [ eigGFuns, eigGCoords ] = intGEP(coeffs,X,f,N,epsilon,m)
%intGEP Evaluate generalized eigenfunctions, associated with spectral 
%parameters X, of an integral operator on [-1,1]. Epsilon is the smoothing
%parameter in the order m rational convolution kernel used to evaluate the 
%approximation eigGFun. f can be any function with nontrivial projection 
%onto the generalized eigenspace associated with X.

%% Set up kernel parameters
[poles,res]=rational_kernel(m,'equi');          %rk poles and residues

%% Multiplication by a(x) and kernel G(x,y) with rank k representation
a=coeffs{1};
G=chebfun2(coeffs{2}); [Gc,Gd,Gr]=cdr(G);
Gd_inv=diag(1./diag(Gd));                   %invert Gd for Woodbury solve

%% Construct "wavepacket" with narrow-band spectral content approximating generalized eigenfunction

% Construct discretizations
[ccn,ccw]=chebpts(N,1);

%Discretized integral operator with low-rank kernel and right-hand side
M_a=a(ccn);           %multiplicative term
M_Gr=ccw.*Gr(ccn)';   %kernel row basis (y) with quad weights
M_Gc=Gc(ccn);         %kernel column basis (x)
f_vals=f(ccn);        %Right hand side

%Solve diagonal plus low-rank with Woodbury formula at each X(i)
for n=1:length(X)
    u=zeros(N,1);
    for j=1:m
        z=X(n)+epsilon*poles(j);
        u1=f_vals./(M_a-z);
        u2=M_Gc./(M_a-z);
        u=u+res(j)*(u1-u2*((Gd_inv+M_Gr*u2)\(M_Gr*u1)));
    end
    psi_pts=imag(u)/pi;
    psi_norm=imag(ccw*(u.*conj(f_vals)))/pi;
end
eigGCoords=sqrt(psi_norm);
eigGFuns=chebfun(psi_pts);
end

