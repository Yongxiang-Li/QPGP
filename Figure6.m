% Sunspot prediction result with QPGP

close all; clear; clc;

addpath('routine');

load sunspot.dat
signal = sunspot(:,2);
% train
Yhat = zeros(288,1);
Yvar = zeros(288,1);
start = 201;
lob = [0.1 0.1 0.5];
upb = [1 2 0.9];
period = (2:1:40)';
modelQPGP0 = fit_QPGP(period, signal(1:200), @regpoly0, @period_sin_gauss_cov, lob, upb);
for i = start:288
    modelQPGP = modelQPGP0;
    modelQPGP.Y = signal(1:i-1);
    % predict
    modelQPGP = pred_QPGP(modelQPGP, i);
    Yhat(i) = modelQPGP.Yhat;
    Yvar(i) = modelQPGP.YVar;
end

x = sunspot(:,1);
figure;plot(x(150:288), signal(150:288),'r')
axis tight;
hold on
plot(x(start:288),Yhat(start:288),'k') 
y1 = (Yhat-2*Yvar)';
y2 = (Yhat+2*Yvar)';
x = sunspot(:,1)';
pic01 = fill([x(start:288),fliplr(x(start:288))],[y1(start:288),fliplr(y2(start:288))],'k');
set(pic01,'edgealpha', 0, 'facealpha',0.4);
xlabel('Year')
ylabel('Number of Sunspots')
legend({'Sunspot Data','Prediction','Confidence Interval'}, 'location', 'NorthWest')






