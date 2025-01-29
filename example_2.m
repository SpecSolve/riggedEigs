%% Example 2. In this example, we compute
%
%       (a) wave packet approximations to the generalized eigenfunctions of
%       the differential operator [Au](x) = -u''(x) on [-inf inf].
%
%       (b) maximum pointwise error in the generalized eigenfunctions over
%       the compact intveral [-10 10] as a smoothing parameter decreases.
%
%  These experiments correspond to Example 2 in the accompanying paper
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.


%% Example 2: (a) wave packet approximations of differential operator
X=1;                                            % evaluation pts
f=@(x) exp(-x.^2)/sqrt(pi);                     % project f(r)
a={@(x) 0.*x, @(x) 0.*x, @(x) -1+0*x};          % differential op
epsilon=0.01;                                   % smoothing parameter
m=8;                                            % smoothing kernel order
N = 2e5;                                        % Disc: 2N+1 Fourier modes

[gef,~,~]=diffGEP(a,X,f,N,epsilon,m);

% true eigenfunctions
g=@(x) cos(sqrt(X)*x);
grid2=linspace(-30,30,600);

figure(2)   % plot generalized eigenfunctions with normalized amplitude
plot(grid2,g(grid2)*max(abs(gef(grid2))),'-k','LineWidth',2)
hold on
set(gca,'ColorOrderIndex',5)
plot(grid2,gef(grid2),'o','MarkerSize',3,'LineWidth',2)
ax = gca; ax.FontSize = 14;
hold on
xlim([grid2(1) grid2(end)])
ylim([-max(abs(gef(grid2))) max(abs(gef(grid2)))])

%% Example 2: (b) point-wise convergence of wave packet on [-10 10]
g=@(x) cos(sqrt(X)*x);
grid2=linspace(-10,10,500);
m=1;

epsilon=logspace(-2,0,10);
relErr=zeros(size(epsilon));
for j=1:length(epsilon)
    [gef,~]=diffGEP(a,X,f,N,epsilon(j),m);
    relErr(j)=max(abs(gef(grid2)/max(gef(grid2))-g(grid2)));
end

figure(3)
set(gca,'ColorOrderIndex',1)
loglog(epsilon,relErr,'o'), hold on
set(gca,'ColorOrderIndex',1)
loglog(epsilon,relErr,'--','LineWidth',2)
ax = gca; ax.FontSize = 14;
