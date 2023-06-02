function [ R ] = period_sin_gauss_cov( delta, theta, p, X, Y )
% p: is the period
% 
    if nargin == 5
        D = X - Y';
        D = pi * D / p;
        sinD = theta*sin(D);
        R = exp(-sinD.^2);
    else
        Y = X;
        D = X - Y';
        D = pi * D / p;
        sinD = theta*sin(D);
        R = exp(-sinD.^2);
        R = R + delta^2 * eye(size(R));
    end
end


