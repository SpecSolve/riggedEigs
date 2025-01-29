%% Example 3. In this example, we compute
%
%       (a) wave packet approximations to the generalized eigenfunctions of
%       an integral operator 
%    
%           [Au](x) = (x^3 - x)u(x) + int_{-1}^1 exp(-x^2-y^2)u(y)dy 
%
%       on the unit interval [-1 1].
%
%       (b) the spectral measure of the integral operator A.
%
%  These experiments correspond to Example 3 in the accompanying paper
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.
%
%   Note that (b) uses intMeas() from the original SpecSolve repository:
%
%               https://github.com/SpecSolve/SpecSolve


%% Example 3: (a) Generalized eigenfunctions of integral operator

X = 0.1; %0.01;                             % evaluation pts
f=chebfun(@(x) (2+x).*cos(2*pi*x));         % project f(x)
a=chebfun(@(x) x.^3-x);                  	% multiplication op
k=chebfun2(@(x,y) exp(-x.^2-y.^2));         % integral perturbation
coeffs={@(x) a(x), @(x,y) k(x,y)};         	% int. op. coeffs for intGEP
epsilon=0.01;                            	% smoothing parameter
m=4;                                      	% smoothing kernel order
N = 1e6;                                    % discretization: N Cheb. pts.

[eigGFuns,eigGCoords]=intGEP(coeffs,X,f,N,epsilon,m);

% roots of cubic (singular support for unperturbed mult. op. gen. eig.)
x=chebfun('x');
r = roots(a(x)-X);

figure(1)   % plot wave packet approximations to generalized eigenfunction
plot(r*ones(1,1e3),linspace(-1e2,1.5e2,1e3),':k','LineWidth',2),
hold on
set(gca,'ColorOrderIndex',2)
plot(eigGFuns,'LineWidth',2);
ax = gca; ax.FontSize = 14;

%% Example 3: (b) Spectral measure of integral operator

f=@(x) (2+x).*cos(2*pi*x);                      % measure wrt f(x)
X=-0.5:0.0005:0.5;                              % evaluation pts
epsilon = 0.001;                                % smoothing parameter
m = 4;                                          % smoothing kernel order
mu=intMeas(coeffs,f,X,0.001,'Order',4);

figure(2)
plot(X,mu,'LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([X(1) X(end)])