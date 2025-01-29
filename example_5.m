%% Example 5. In this example, we compute wave packet approximations to the 
%   generalized eigenfunctions of the Laplacian on the unbounded strip 
%   [-inf inf] x [-1 1] with homogenous Dirichlet boundary conditions on 
%   the boundary, u(x,-1)=u(x,1)=0.
%
%  These experiments correspond to Example 5 in the accompanying paper
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.


%% Problem parameters: run this section to set problem parameters

X= pi*(pi-1); %pi^2-0.2 % pi*(pi+1)         % spectral parameter
fx=@(x) exp(-x.^2);                         % project f(x,y) = fx(x)fy(y)
fy=@(y) cos(pi*y/2) + sin(pi*y);            % project f(x,y) = fx(x)fy(y)
epsilon=0.1;                                % smoothing parameter
m=4;                                        % smoothing kernel order
N1=900;                                     % 2*N1+1 Fourier modes in x
N2=110;                                     % N2 UltraS modes in y

%% Discretize Laplacian: run this section to build discretizations

% differentiation matrices on x-slices
scale=15;
map=@(t) scale*1i*(1-exp(1i*t))./(1+exp(1i*t));
singular_weight=chebfun(@(t) 0.5*exp(-1i*t).*(1+exp(1i*t)).^2,[-pi,pi],'trig')/scale;
Dx=1i*spdiags((-N1:N1)',0,2*N1+1,2*N1+1);
weight_coeffs=trigcoeffs(singular_weight,4*N1+1);
SW=sptoeplitz(weight_coeffs(2*N1+1:end),weight_coeffs(2*N1+1:-1:1));
Dx=SW*Dx;
D2x=Dx*Dx;

% differentiation and conversion matrices on y-slices
D2y=ultraS.diffmat(N2,2);
S02=ultraS.convertmat(N2,0,1);

% discretize right-hand side on x slices
FxFun=chebfun(@(t) fx(map(t)),[-pi,pi],'trig');
FxCoeffs=sqrt(2*pi)*trigcoeffs(FxFun,2*N1+1);

% discretize right hand side y slices x slices
FyFun=chebfun(@(y) fy(y));
FyCoeffs=S02*chebcoeffs(FyFun,N2);

% assemble right-hand side matrix from low-rank factors
FMat=FxCoeffs*transpose(FyCoeffs);

% boundary bordering for DBCs on y slices
u1=(-1).^(0:(N2-1)); u2=ones(1,N2);
D2y=[u1; u2; D2y(1:N2-2,:)];
S02=[zeros(2,N2); S02(1:N2-2,:)];
FMat=[FMat(:,1:size(FMat,2)-2) zeros(size(FMat,1),2)];

% assemble for linear system via Kronecker identity
I=speye(2*N1+1);
FVec=FMat(:);
H=-kron(S02,D2x)-kron(D2y,I);
H=real(H);
Ikron=kron(S02,I);

% reorder for reduced bandwidth
p=symrcm(H); H=H(p,p);
Ikron=Ikron(p,p); FVec=FVec(p);

%% Compute and plot generalized eigenfunction

% rational convolution kernel for mth order convergence
[poles,res] = rational_kernel(m,'equi');

C=zeros(size(H,2),1);
for j=1:m
    C=C+res(j)*((H-(X+poles(j)*epsilon)*Ikron)\FVec);
end
pp(p)=1:length(p); C=C(pp);   % restore original indexing

% transform coeffs to values and map back to strip
C=reshape(C,2*N1+1,N2);
U=chebtech1.coeffs2vals(transpose(C));
U=trigtech.coeffs2vals(transpose(U))/sqrt(2*pi);
Ufun=chebfun2(transpose(U),[-pi pi -1 1]);
map=@(x) -1i*log((1+1i*x/scale)./(1-1i*x/scale));

% plot generalized eigenfunction on truncated strip
figure(1)
gef1=chebfun2(@(x,y) (imag(Ufun(map(x),y))/pi), [-15 15 -1 1]);
plot(gef1), view(2)
ax = gca; ax.FontSize = 14;

%% Compute measure analytically for 2D Laplace in strip to illustrate 
%   singular points corresponding to multiplicity changes.

nu = 1e-2;
X1 = linspace(0,0.25*pi^2-nu,200);
X2 = linspace(0.25*pi^2+nu,pi^2-nu,200);
X3 = linspace(pi^2+nu,15,200);
rho1 = @(x) (x>pi^2/4).*( exp(-(x-pi^2/4)/4) ).^2 ./ sqrt(x-pi^2/4) / 2;
rho2 = @(x) (x>pi^2).*( exp(-(x-pi^2)/4) ).^2 ./ sqrt(x-pi^2)/ 2;
rho = @(x) rho1(x) + rho2(x);

Y = linspace(0,5,250);

figure(2)
plot(X1,rho(X1),'k','LineWidth',2), hold on
plot(X2,rho(X2),'k','LineWidth',2)
plot(X3,rho(X3),'k','LineWidth',2)
plot(0.25*pi^2*ones(size(Y)),Y,'k:','LineWidth',2)
plot(pi^2*ones(size(Y)),Y,'k:','LineWidth',2)
ax = gca; ax.FontSize = 14;
