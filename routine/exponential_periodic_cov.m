function [ R ] = exponential_periodic_cov( rho, theta, p, X, Y )
% p: is the period
% 
    if nargin ~= 5
        Y = X;
    end
    D = X - Y';
    lambda = 1/(rho*p);
    K1 = exp(-lambda^2 * D.^2 / 2);
    
    D = pi * D / p;
    sinD = theta*sin(D);
    K2 = exp(-sinD.^2);
    R =  K1 .* K2;
end


