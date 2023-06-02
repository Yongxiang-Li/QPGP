function statistic_MLPE = mlpe(signal,period)
    K = length(signal);
    N = floor(K/period);
    m = ceil(K/period);
    signal(K+1:m*period) = 0;
    Y = reshape(signal,period,m);
    q = [sum(Y(1:mod(K,period),:),2)/(N+1);sum(Y(mod(K,period)+1:period,:),2)/N]';
    Es = N*sum(q.^2)+sum(q(1:K-N*period).^2);
    Phi0 = sum(signal.^2);
    statistic_MLPE = Es-period/K*Phi0;

