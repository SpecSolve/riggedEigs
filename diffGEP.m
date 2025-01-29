function [ eigGFuns, eigGCoords, eigGFuns_map ] = diffGEP(a,X,f,N,epsilon,m)
%diffGEP Evaluate generalized eigenfunctions, associated with spectral 
%parameters X, of a differential operator on the real line. Epsilon is the 
%smoothing parameter in the order m rational convolution kernel used to 
%evaluate the approximation eigGFun. f can be any function with nontrivial 
%projection onto the generalized eigenspace associated with X.

%% set up kernel and discretization parameters
[poles,res]=rational_kernel(m,'equi');          %rk poles and residues
scale=10;                                       %scale map to unit circle
M=200;                                          %max bandwidth

%% Set up chebfuns for variable coeffs and rhs
map=@(t) scale*1i*(1-exp(1i*t))./(1+exp(1i*t));
singular_weight=chebfun(@(t) 0.5*exp(-1i*t).*(1+exp(1i*t)).^2,[-pi,pi],'trig')/scale;
F=chebfun(@(t) f(map(t)),[-pi,pi],'trig');
G=chebfun(@(t) scale*2*exp(1i*t).*f(map(t))./(1+exp(1i*t)).^2,[-pi,pi],'trig');
var_coeffs=zeros(2*M+1,length(a));
for k=1:length(a)
    var_coeffs(:,k)=trigcoeffs(chebfun(@(t) a{k}(map(t)),[-pi,pi],'trig'),2*M+1);
end

%% Construct "wavepacket" with narrow-band spectral content approximating generalized eigenfunction

% Construct differentiation matrix
D=1i*spdiags((-N:N)',0,2*N+1,2*N+1);
weight_coeffs=trigcoeffs(singular_weight,4*N+1);
SW=sptoeplitz(weight_coeffs(2*N+1:end),weight_coeffs(2*N+1:-1:1));
D=SW*D;

% Construct variable coeff diff op matrix
for k=1:length(a)
    if 2*N>M
        var_coeffsT=[sparse(2*N-M,1);
                    var_coeffs(:,k);
                    sparse(2*N-M,1)];
        C=sptoeplitz(var_coeffsT(2*N+1:end),var_coeffsT(2*N+1:-1:1));
    else
        C=sptoeplitz(var_coeffs(2*N+1:end,k),var_coeffs(2*N+1:-1:1,k));
    end
    if k==1
        L=C; % extra "if" needed, otherwise get memory error
    else
        L=L+C*(D^(k-1));
    end
end

% Shifts
I=speye(2*N+1);

% Discretize right-hand-side
F_coeffs=sqrt(2*pi)*trigcoeffs(F,2*N+1);
G_coeffs=sqrt(2*pi)*trigcoeffs(G,2*N+1);

% Solve banded, shifted linear systems at each X(i)
for i=1:length(X)
    u_coeffs=zeros(size(F_coeffs));
    for ii=1:m
        u_coeffs=u_coeffs+res(ii)*((L-(X(i)+poles(ii)*epsilon)*I)\F_coeffs);
    end
    psi_coeffs=u_coeffs;
    psi_norm=imag(G_coeffs'*u_coeffs)/pi;
end
psi_coeffs=[psi_coeffs(N+1:end); psi_coeffs(1:N)];
psi_vals=fft(psi_coeffs); psi_vals=imag([psi_vals(N+1:end); psi_vals(1:N)])/pi;
eigGFuns_map=chebfun(psi_vals,[-pi,pi],'trig');
eigGCoords=sqrt(psi_norm);
map_inv=@(x) -1i*log((1+1i*x/10)./(1-1i*x/10));
eigGFuns=@(x) eigGFuns_map(real(map_inv(x)));
end