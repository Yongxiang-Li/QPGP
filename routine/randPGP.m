function [ Y ] = randPGP(N, p, sigma, delta, theta, corr)
    %RANDQPGP generate QPGP signal
    k = ceil(N/p);
    [Uo, So] = eig(ones(k));    So = diag(So);
    r0 = corr(delta, theta, p, (1:p)', 1);
    r0(1) = r0(1) + eps;    
    Sr = real(fft(r0));    Sr(Sr<0) = 0;
    Z = randn(p, k);
    S = sigma^2 * (Sr * So' + delta^2);
%     Y = fft(conj(((fft(Z*Uo).*sqrt(S))/sqrt(p))*Uo'))/sqrt(p); % YY = Ur*(((Ur'*Z*Uo).*sqrt(S))*Uo');
    Y = ifft(((fft(Z*Uo).*sqrt(S)))*Uo'); % YY = Ur*(((Ur'*Z*Uo).*sqrt(S))*Uo');
    Y = real(Y(1:N))';
    if nargout == 0,    figure; plot(Y);    end
end

