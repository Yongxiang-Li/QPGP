% The full likelihood and the first composite likelihood

close all;
clear;
clc;
rng('default');

addpath('routine');

sigma = 1;
theta = 1;
omega = 0.96;
P = 50;
N = 4000; 
SNR = -18;
% generate QPGP signal
delta = 10^(-SNR/20);
noiseSignal = randQPGP(N, P, sigma, delta, theta, omega, @period_sin_gauss_cov); 

% init analysis
period = (10:1:400)';

% QPGP
lob = [0.1 0.1 0.9];
upb = [3 3 0.99];

modelQPGP = fit_QPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
modelQPGP1 = fit_QPGP_L1(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
modelQPGP2 = fit_QPGP_L1Normalized(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);

figure;
subplot(3,1,1)
plot(period, modelQPGP.likelihood)
title('Full Log-likelihood $\ell$' ,'interpreter','latex')
xlabel('\itp')
ylabel('$\ell$' ,'interpreter','latex')
axis tight;
subplot(3,1,2)
plot(period, modelQPGP1.likelihood)
title('First Composite Log-likelihood $\ell_1$','interpreter','latex')
xlabel('\itp')
ylabel('$\ell_1$' ,'interpreter','latex')
axis tight;
subplot(3,1,3)
plot(period, modelQPGP2.likelihood)
title('First Normalized Composite Log-likelihood $\ell_1$','interpreter','latex')
xlabel('\itp')
ylabel('$\ell_1/(kp)$' ,'interpreter','latex')
axis tight;


