function [ Q, C, V, S, m] = nrc( noiseSignal, N )
%NRC Summary of this function goes here
%   Detailed explanation goes here
    
    Y =  cut_signal(N,noiseSignal);
    m = size(Y,2);%  k = floor(m/n);
    Yhat = mean(Y,2)';   l = N*m;
    C = Yhat*Yhat'/length(Yhat);
    V = sum(noiseSignal(1:l).^2)/l - C;
    Q = C - V/(m-1);
    S = m*V/(m-1);
end
