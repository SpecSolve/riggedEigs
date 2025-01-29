%% Example 4. In this example, we compute wave packet approximations to the
%       generalized eigenfunctions of Schrodinger operators on [-inf inf]:
%
%       (a) potential v(x) = -5 * cos(x/4) * exp(-x^2/32)
%
%       (b) potential v(x) = -10 / (2 + x^2)
%
%   We compare with a scheme that uses domain truncation + perfectly
%   matched layers. These experiments correspond to Example 4 in the paper
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.

%% Example 4: (a) potential v(x) = -5 * cos(x/4) * exp(-x^2/32)

X=5;                                            % evaluation pts
f=@(x) exp(-x.^2)/sqrt(pi);                     % project f(r)
v = @(x) cos(x/2).*(-5*exp(-x.^2/32));          % potential
a={v, @(x) 0.*x, @(x) -1+0*x};                  % differential op
epsilon=0.01;                                   % smoothing parameter
m=4;                                            % smoothing kernel order
N = 2e5;                                        % Disc: 2N+1 Fourier modes

% compute with diffGEP
[gef,~,~]=diffGEP(a,X,f,N,epsilon,m);
grid2=linspace(-30,30,1e3);

% compare with domain truncation + perfectly matched layer (PML)
L = 31;                                         % truncate domain to [-L L]
dpml = 1.25;                                 	% PML depth
alpha = 1;                                      % PML strength
s = @(x) 0.*(abs(x)<=L) + alpha*(abs(x)>L).*((abs(x)-L)/dpml).^2;   %PML
spml = chebfun(@(x) 1./(1+1i*s(x)),[-L-dpml L+dpml],'splitting','on');
dspml = diff(spml);                             % PML derivative

% grid
N = 250*L;                                      % number of gridpts
x = linspace(-L-dpml,L+dpml,N+2);               % computational grid (w/boundary)
xpts = x(2:end-1);                          	% computational grid (interior)
h = (2*(L+dpml)) / (N+1);                      	% grid spacing 

% finte-difference discretization on grid
D1 = spdiags([-ones(N,1) zeros(N,1) ones(N,1)], -1:1, N,N) / (2*h);
D2 = spdiags([ones(N,1) -2*ones(N,1) ones(N,1)], -1:1, N,N) / h^2;
ML = spdiags(spml(xpts'),0,N,N);
MLD = spdiags(dspml(xpts'),0,N,N);
H1 = - (ML*D2 + MLD*D1)  - X*speye(N);
vv = @(x) 0.*(abs(x)>L) + (abs(x)<=L).*v(x);
MV = spdiags(vv(xpts'),0,N,N);
H = H1 + MV;

% right-hand side on grid
f = @(x) exp(-x.^2)/sqrt(pi);
b = 0.*(abs(xpts')>L) + (abs(xpts')<=L).*f(xpts');

% solve discretized system
u = H \ b;

figure(1)   % plot potential function
plot(xpts',v(xpts'),'LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([-30 30])

figure(2)   % plot generalized eigenfunctions with normalized amplitude
set(gca,'ColorOrderIndex',2)
plot(xpts',(abs(xpts')<=L).*imag(u)/max(abs(imag(u))),'ok','MarkerSize',4,'LineWidth',2)
hold on
plot(grid2,gef(grid2)/max(abs(gef(grid2))),'LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([grid2(1) grid2(end)])

%% Example 4: (b) potential v(x) = -10 / (2 + x^2) 

X=5;                                            % evaluation pts
f=@(x) exp(-x.^2)/sqrt(pi);                     % project f(r)
v = @(x) -10./(2+x.^2);                         % potential (a)
a={v, @(x) 0.*x, @(x) -1+0*x};                  % differential op
epsilon=0.01;                                   % smoothing parameter
m=4;                                            % smoothing kernel order
N = 2e5;                                        % Disc: 2N+1 Fourier modes

% compute with diffGEP
[gef,~,~]=diffGEP(a,X,f,N,epsilon,m);
grid2=linspace(-30,30,1e3);

% compare with domain truncation + perfectly matched layer (PML)
L = 31;                                         % truncate domain to [-L L]
dpml = 1.25;                                 	% PML depth
alpha = 1;                                      % PML strength
s = @(x) 0.*(abs(x)<=L) + alpha*(abs(x)>L).*((abs(x)-L)/dpml).^2;   %PML
spml = chebfun(@(x) 1./(1+1i*s(x)),[-L-dpml L+dpml],'splitting','on');
dspml = diff(spml);                             % PML derivative

% grid
N = 250*L;                                      % number of gridpts
x = linspace(-L-dpml,L+dpml,N+2);               % computational grid (w/boundary)
xpts = x(2:end-1);                          	% computational grid (interior)
h = (2*(L+dpml)) / (N+1);                      	% grid spacing 

% finte-difference discretization on grid
D1 = spdiags([-ones(N,1) zeros(N,1) ones(N,1)], -1:1, N,N) / (2*h);
D2 = spdiags([ones(N,1) -2*ones(N,1) ones(N,1)], -1:1, N,N) / h^2;
ML = spdiags(spml(xpts'),0,N,N);
MLD = spdiags(dspml(xpts'),0,N,N);
H1 = - (ML*D2 + MLD*D1)  - X*speye(N);
MV = spdiags(v(xpts'),0,N,N);
H = H1 + MV;

% right-hand side on grid
f = @(x) exp(-x.^2)/sqrt(pi);
b = 0.*(abs(xpts')>L) + (abs(xpts')<=L).*f(xpts');

% solve discretized system
u = H \ b;

figure(1)   % plot potential function
plot(xpts',v(xpts'),'LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([-30 30])

figure(2)   % plot generalized eigenfunctions with normalized amplitude
set(gca,'ColorOrderIndex',2)
plot(xpts',(abs(xpts')<=L).*imag(u)/max(abs(imag(u))),'ok','MarkerSize',4,'LineWidth',2)
hold on
plot(grid2,gef(grid2)/max(abs(gef(grid2))),'LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([grid2(1) grid2(end)])
