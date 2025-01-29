%% Example 1. In this example, we compute
%
%       (a) wave packet approximations to the generalized eigenfunctions of
%       the multiplication operator [Pu](x) = (x^3 - x)u(x) on [-1 1].
%
%       (b) error in the approximate eigencoordinates of a test function
%       phi = (1+x)cos(pi x) as smoothing parameter decreases.
%
%  These experiments correspond to Example 1 in the accompanying paper
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.


%% Example 1: (a) Generalized eigenfunctions of multiplication operator

X = 0.1; %0.01;                             % evaluation pts
f=chebfun(@(x) (2+x).*cos(2*pi*x));         % project f(x)
a=chebfun(@(x) x.^3-x);                  	% multiplication op
k=chebfun2(@(x,y) 0.*x + 0.*y);             % integral part = 0
coeffs={@(x) a(x), @(x,y) k(x,y)};         	% mult. op. coeffs for intGEP
epsilon=0.01;                            	% smoothing parameter
m=1;                                      	% smoothing kernel order
N = 1e6;                                    % discretization: N Cheb. pts.

[eigGFuns,eigGCoords]=intGEP(coeffs,X,f,N,epsilon,m);

% generalized eigenfunction of multiplication operator at spectral point X
x=chebfun('x');
r = roots(a(x)-X);
da = abs(diff(a));
g=f(r(1))*dirac(x-r(1))/da(r(1));
for j = 2:length(r)
    g = g + f(r(j))*dirac(x-r(j))/abs(da(r(j)));
end

figure(1)
plot(eigGFuns,'LineWidth',2)    % plot wave packet approximations
hold on
plot(g,'r','LineWidth',2)       % plot true Dirac (red arrows at support)
ax = gca; ax.FontSize = 14;
hold on

% error in computed coordinate of phi
phi = chebfun(@(x) (1+x).*cos(pi*x)); % 1+exp(-x.^2)); %
coordErr=abs(eigGFuns'*phi-g'*phi);
fprintf('relative error in generalized coordinates is %.2e.\n',coordErr);

%% Example 1: (b) Convergence of generalized eigencoordinates (weak conv.)

X = 0.1; %0.01; %0.001                      % evaluation pts
m=1;                                      	% smoothing kernel order

% generalized eigenfunction of multiplication operator at spectral point X
x=chebfun('x');
r = roots(a(x)-X);
da = abs(diff(a));
g=f(r(1))*dirac(x-r(1))/da(r(1));
for j = 2:length(r)
    g = g + f(r(j))*dirac(x-r(j))/abs(da(r(j)));
end

% compute error in generalized eigenfunction coordinate as epsilon->0
epsilon=logspace(-3,0,10);
relErr=zeros(size(epsilon));
for j=1:length(epsilon)
    [eigGFuns,eigGCoords]=intGEP(coeffs,X,f,N,epsilon(j),m);
    relErr(j)=abs(eigGFuns'*phi-g'*phi)/abs(g'*phi);
end

figure(2)
set(gca,'ColorOrderIndex',1)
loglog(epsilon,relErr,'o')
hold on
set(gca,'ColorOrderIndex',1)
loglog(epsilon,relErr,'--','LineWidth',2)
ax = gca; ax.FontSize = 14;