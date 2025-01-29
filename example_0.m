%% Example 0. This example analytically computes the spectral measure of
%
%       1) the multiplication operator [Pu](x) = (x^3 - x)u(x) on [-1 1].
%       2) the differential operator [Du](x) = -u''(x) on [-inf, inf] 
%
%   The branches of the multiplier inverse used to construct the spectral
%   measure are also plotted, as in Section 2.3 of the accompanying paper:
%
%       "Computing Generalized Eigenfunctions in a Rigged Hilbert Space"
%
%   by Matthew J. Colbrook, Andrew Horning, and Tianyiwa Xie.


%% 1) multiplication operator

% first compute and plot inverse branches of mult. op. x^3-x
p = chebfun(@(x) x.^3-x);
pinv = @(y) roots(p-y);

ygrid = linspace(min(p),max(p),5000);
pgrid = pinv(ygrid);
pinv1 = pgrid(2,1:2500);
pinv2 = pgrid(1,1:2500);
pinv3 = pgrid(2,2501:end);
pinv4 = pgrid(1,2501:end);

figure(1)   % plot inverse branches
C=colororder; C(3,:) = C(1,:); colororder(C)
plot(ygrid(1:2500),pinv1,'LineWidth',2), hold on
plot(ygrid,[pinv2 pinv3],'LineWidth',2)
plot(ygrid(2501:end),pinv4,'LineWidth',2)
plot(zeros(1,100),linspace(-1,1,100),':k','LineWidth',2)
plot(min(p)*ones(1,100),linspace(-1,1,100),':k','LineWidth',2)
plot(max(p)*ones(1,100),linspace(-1,1,100),':k','LineWidth',2)
ax = gca; ax.FontSize = 14;

% then compute spectral measure with respect to (2+x)*cos(2*pi*x)
f = chebfun(@(x) (2+x)*cos(2*pi*x));
dp = diff(p);

meas1 = f(pinv1).^2 ./ abs(dp(pinv1));
meas2 = f(pinv2).^2 ./ abs(dp(pinv2));
meas3 = f(pinv3).^2 ./ abs(dp(pinv3));
meas4 = f(pinv4).^2 ./ abs(dp(pinv4));

figure(2)   % plot spectral measure
semilogy(ygrid(1:2500),meas1+meas2,'k','LineWidth',2), hold on
semilogy(ygrid(2501:end),meas3+meas4,'k','LineWidth',2)
plot(zeros(1,100),logspace(-2,10,100),':k','LineWidth',2)
plot(min(p)*ones(1,100),logspace(-2,10,100),':k','LineWidth',2)
plot(max(p)*ones(1,100),logspace(-2,10,100),':k','LineWidth',2)
ax = gca; ax.FontSize = 14;
ylim([1e-2 1e6])

%% 2) differential operator

% first compute and plot inverse branches of Fourier multiplier x^2
p = chebfun(@(x) x.^2, [-2 2]);
pinv = @(y) roots(p-y);

ygrid = linspace(min(p),max(p),1e5);
pgrid = pinv(ygrid);
pinv1 = pgrid(2,:);
pinv2 = pgrid(1,:);

figure(3)   % plot inverse branches
plot(ygrid,pinv1,'LineWidth',2), hold on
plot(ygrid,pinv2,'LineWidth',2)
plot(min(p)*ones(1,100),linspace(-2,2,100),':k','LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([-1 4])

% then compute spectral measure with respect to exp(-pi*x^2)
f = chebfun(@(x) exp(-pi*x.^2), [-inf inf]);
dp = diff(p);

meas1 = f(pinv1).^2 ./ abs(dp(pinv1));
meas2 = f(pinv2).^2 ./ abs(dp(pinv2));

figure(4)   % plot spectral measure
semilogy(ygrid,meas1+meas2,'k','LineWidth',2), hold on
plot(zeros(1,100),logspace(-15,5,100),':k','LineWidth',2)
ax = gca; ax.FontSize = 14;
xlim([-1 4])
ylim([1e-12 1e2])
