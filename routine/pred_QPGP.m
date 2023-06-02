function [model, Yhat] = pred_QPGP(model, X) 
%  Predection for QPGP
% example: [dmodel, Yhat] = pred_SPGP(dmodel, 1:dmodel.period) 
    model.X = X;    
    n = length(model.Y);    q = size(model.Gamma,2);    p = model.period;
    delta = model.thetahat(1);     theta = model.thetahat(2);    omega = model.thetahat(3);
  	
    if p < n
        k = floor(n/p);    p1 = n - p*k;
        Y = zeros(p,k);    Y(1:k*p) = model.Y(1:k*p);    Ys = model.Y(k*p+1:n);
        G = model.Gamma((1:p)');    Gs = G(1:p1,:);
        [So, Uo] = kmseig(k, omega);
        r0 = model.corr(delta, theta, p, (1:p)', 1);    r0(1) = r0(1) + eps;
        Sr = real(fft(r0));    Sr(Sr<0) = 0; 
        S = Sr * So' + delta^2;    Lambda = 1./S;
        wt = omega.^(abs((1:k)'-ceil(X'/p)));
        r = model.corr(delta, theta, p, (1:p)', X);    FFTYUo = fft(Y*Uo);
        UoW = Uo'*omega.^(k:-1:1)';    UoWt = Uo'*wt;    UoOnes = Uo'*ones(k,1);  
        FFTR = fft(r);    FFTG = fft(G);
        SGY = real(sum(conj(FFTR).*((FFTYUo.*Lambda)*UoWt))/p);
        SGF = sum(conj(FFTR).*((Lambda*(UoWt.*UoOnes)).*FFTG))/p;
        f = model.regr(mod(X,p));
        mu = (f-SGF')*model.betahat + SGY';
        
        Yd = ifft(Sr.*((fft(Y*Uo).*Lambda)*UoW));
    	Yd = real(Ys - Yd(1:p1)); % Yd = Ys - real(UrRd'*(UrYUo.*Lambda)*UoW);
        Fd = ifft((Sr.*(Lambda*(UoW.*UoOnes))).*fft(G));
      	Fd = real(Gs - Fd(1:p1,:)); % Fd = Fs - real(UrRd'*((Lambda*(UoW.*UoOnes)).*UrGamma)); 
        rs = model.corr(delta, theta, p, (1:p1)', X) .* omega.^abs(ceil((1:p1)'/p)+k - ceil(X'/p)); % \boldsymbol{\gamma}_{\bullet}
        SSG = ifft((Sr.*(Lambda*(UoWt.*UoW))).* FFTR);
        rd = real(rs - SSG(1:p1,:)); % \boldsymbol{\gamma}_{\bullet}
        
        SGG = real(sum(conj(FFTR).*((Lambda*(UoWt.*UoWt)).*FFTR))/p);
        sigma1 = 1 + delta^2 - SGG';
        
        
        % PI matrix
        if p1 >0
            r0(1) = r0(1) + delta^2; 
            r0 = r0 - ifft((Lambda*(UoW.^2)).*(Sr.^2));
            L = toep_chol(r0(1:p1)); % S = chol((Pi+eps*eye(size(Pi))), 'lower');
            rd_L = L \ rd;
            Yhat = mu + rd_L' * (L \ (Yd-Fd*model.betahat));
            sigma2 = sum(rd_L .* rd_L)';
            Yvar = model.sigmahat * sqrt(sigma1 - sigma2);
        else
            Yhat = mu;
            Yvar = model.sigmahat * sqrt(sigma1);
        end
    else
        C = model.corr(delta, theta, p, (1:n)', (1:n)') .* omega.^abs(ceil((1:n)'/p)-ceil((1:n)/p));
        C = C + delta^2*eye(size(C));    L = chol(C, 'lower');
        F = repmat(model.Gamma, ceil(n/p), 1);    F = F(1:n,:);
        gamma = L' \ (L \ (model.Y - F*model.betahat));
        Cr = model.corr(delta, theta, p, X, (1:n)') .* omega.^abs(ceil(X/p) - ceil((1:n)/p));
        f = model.regr(mod(X,p));
        Yhat = f * model.betahat + Cr*gamma;
        Cr_L = L \ Cr';
        Yvar = model.sigmahat * sqrt((1 + delta^2 - sum(Cr_L.*Cr_L)'));
    end
    model.Yhat = Yhat;    model.YVar = Yvar;
end % figure; plot(model.X, model.Yhat)
