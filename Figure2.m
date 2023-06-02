% Typical waveforms of pseudo-periodic signals generated from QPGP

close all;
clear;
clc;
rng('default');

addpath('routine');

sigma = 1;    theta = 5;    P = 100;    N = 1000;  
% synthetic QPGP signal a
omega = 0.5;    delta = 0;
signal_a = randQPGP(N, P, sigma, delta, theta, omega, @period_sin_gauss_cov); 
% synthetic QPGP signal b
omega = 0.9;    delta = 0;
signal_b = randQPGP(N, P, sigma, delta, theta, omega, @period_sin_gauss_cov);
% synthetic QPGP signal b
omega = 0.9;    SNR = -5;
delta = 10^(-SNR/20);
signal_c = randQPGP(N, P, sigma, delta, theta, omega, @period_sin_gauss_cov);

figure;
subplot(3,1,1)
plot(signal_a)
title('(a) Synthetic signal with $\omega$ = 0.5, $\delta$ = 0(SNR= $\infty$dB)' ,'interpreter','latex')
xlabel('\itn')
ylabel('Amplitude')
subplot(3,1,2)
plot(signal_b)
title('(b) Synthetic signal with $\omega$ = 0.9, $\delta$ = 0(SNR = $\infty$dB)' ,'interpreter','latex')
xlabel('\itn')
ylabel('Amplitude')
subplot(3,1,3)
plot(signal_c)
title('(c) Synthetic signal with $\omega$ = 0.9, $\delta$ = 1.78(SNR = -5dB)' ,'interpreter','latex')
xlabel('\itn')
ylabel('Amplitude')

