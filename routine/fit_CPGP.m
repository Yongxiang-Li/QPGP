function [ fit ] = fit_CPGP( P, Y, regr, corr, lob, upb, theta0 )
% Circulant Periodic Gaussian Process Modeling
% Inputs:
%       P: searching grids of the period (column vector)
%       Y: periodic data (column vector)
%       regr: regressor, e.g. @regpoly0
%       corr: correlation function, e.g. @period_sin_gauss_cov
%       lob, upb: lower and upper bounds of the parameter \delta and \theta
%       theta0: initial value for theta without optimizing theta
% Outputs:
%       fit: fitted CPGP model
    
    F = regr((1:length(Y))'/length(Y));
    data = struct('corr',corr, 'regr',regr, 'F',F, 'P', P, 'Y',Y, 'FF',F'*F, 'YY',Y'*Y, 'FY',F'*Y);
    index = (P == round(P));    data.intP = P(index);
  	q = size(F,2);
    [data.Yb, data.Ys] = segment_mean(Y, data.intP);
    data.Gb = ones(max(data.intP), length(data.intP), q); % mean of G
    data.Gs = ones(max(data.intP), length(data.intP), q); % tail of G
    if q > 1
        for j = 2 : q
            [data.Gb(:,:,j), data.Gs(:,:,j)] = segment_mean(F(:,j),data.intP);
        end
    end
    data.Gb = permute(data.Gb, [1,3,2]);
    data.Gs = permute(data.Gs, [1,3,2]);
    
    if nargin == 7
      	theta = theta0;
        if isempty(lob) && isempty(upb) % without optimization
            [~, fit] = likelihood_at_interger(theta0, data); 
        else
            [theta, ~, fit, ~] = boxmin(@likelihood_at_interger, theta0, lob, upb, data);
        end
    else % init theta0 via grid search
        u = linspace(0,1,5);    [x1, x2] = meshgrid(u(2:end-1));
        thetas = lob + (upb-lob).*[x1(:),x2(:)];
        objs = nan(size(thetas,1),1);
        fits = [];
        for i = 1 : size(thetas,1)
            [objs(i), fit] = likelihood_at_interger(thetas(i,:), data);
            fits = [fits; fit]; % figure; plot(fit.likelihood, '-*');
        end 
        [~, index] = min(objs);
        theta0 = thetas(index,:);
        [theta, ~, fit, ~] = boxmin(@likelihood_at_interger, theta0, lob, upb, data);
    end
    
    if any((P ~= round(P)))
        fit = likelihood_at_decimal(theta, fit);
    end
    fit.theta0 = theta0;
    fit.thetahat = theta;
    [~, index] = max(fit.likelihood);
	fit.betahat = fit.beta(index,:);
    fit.sigmahat = fit.sigma(index);
    fit.period = fit.P(index);
end

function [Yb, Ys] = segment_mean(Y, P)
    n = length(Y);    Y(end+max(P)) = 0; 
    Yseg = nan(max(P), length(P)); % sum of Y
    Ys = zeros(max(P), length(P)); % tail of Y
    for i = length(P) : -1 : 1
        if ~isnan(Yseg(1,i)), continue; end
        p = P(i);    p1 = n - p*floor(n/p);
        Yseg(1:p,i) = sum(reshape(Y(1:p*ceil(n/p)), p, ceil(n/p)), 2);
        Ys(1:p1,i)  = Y(n-p1+1:n);
        for j = 1 : i-1
            if ~isnan(Yseg(1,j)), continue; end
            q = P(j);    q1 = n - q*floor(n/q);
            if mod(p,q)==0
                Yseg(1:q,j) = sum(reshape(Yseg(1:p,i), q, p/q), 2);
                Ys(1:q1,j)  = Y(n-q1+1:n);
            end
        end
    end
    Yb = (Yseg - Ys) ./ floor(n./P');
end

function [fit] = likelihood_at_decimal(para, fit) % likelihood at decimal
    delta = para(1);     theta = para(2:end); % Initialize
    n = length(fit.Y);    q = size(fit.F,2);   
    for i = 1 : length(fit.P)
        if ~isnan(fit.likelihood(i)), continue; end
        T = fit.P(i);    p = str2num(regexprep(rats(T, 30),'/.*', ''));
        k = floor(n/p);    p1 = n - p*k;
        if p>n
            r = fit.corr(delta, theta, T, (1:n)', 1);    r(1) = r(1) + delta^2;
            L = toep_chol(r);
            Ft = L \ fit.F;    Yt = L \ fit.Y; 
            [Q G] = qr(Ft,0);    B = G \ (Q'*Yt);
            rho = Yt - Ft*B;     sigma2 = sum(rho.^2)/n;
            fit.likelihood(i) = -(n*log(sigma2) + 2*sum(log(diag(L))) + n + n*log(2*pi))/2;
        else
            Yb = mean(reshape(fit.Y(1:k*p), p, k), 2); % mean of Y
            Ys = fit.Y(1+p*k:end); % tail of Y
            Gb = nan(p, q); % mean of G
            if q == 1
                Gb = mean(reshape(fit.F(1:k*p), p, k), 2);
            else
                for j = 1 : q
                    Gb(:,j) = mean(reshape(fit.F((1:k*p),j), p, k), 2);
                end
            end
            Gs = fit.F(1+p*k:end,:); % tail of G
            GG = fit.F'*fit.F-Gs'*Gs;    YY = fit.Y'*fit.Y-Ys'*Ys;    GY = fit.F'*fit.Y-Gs'*Ys;
            R0.c = fit.corr(delta, theta, p, (1:p)', 1);    R0.eig = abs(fft(R0.c))+eps;
            Rd.c = k*R0.c/delta^2 + eye(p,1);    Rd.eig = real(fft(Rd.c));
            RRd.eig = (k/delta^2) * (R0.eig ./ Rd.eig);
            fftGb = fft(Gb);    fftYb = fft(Yb);
            Gbl = ifft(Rd.eig .\ fftGb);    Ybl = ifft(Rd.eig .\ fftYb);
            SGG = (GG+k*Gb'*(Gbl-Gb))/delta^2;
            SGY = (GY+k*(Gbl-Gb)'*Yb)/delta^2;
            SYY = (YY+k*Yb'*(Ybl-Yb))/delta^2;
            if p1==0
                B = SGG \ SGY; % beta
                sigma2 = (SYY - B'*SGG*B)/n; % sigma
                fit.likelihood(i) = -(n*log(sigma2*delta^2) + sum(log(Rd.eig)) + n + n*log(2*pi))/2;
            else
                RRdGb = ifft(RRd.eig .* fftGb);    RRdYb = ifft(RRd.eig .* fftYb);
                rpi = R0.c - ifft(RRd.eig .* R0.eig) + delta^2*eye(p,1);
                S = toep_chol(rpi(1:p1));
                Gd = Gs - RRdGb(1:p1,:); % Gamma dot
                Yd = Ys - RRdYb(1:p1); % Y dot
                Gdl = S \ Gd;    GGd = Gdl'*Gdl;
                Ydl = S \ Yd;    YYd = Ydl'*Ydl;
                GYd = Gdl'*Ydl;
                B = (SGG + GGd) \ (SGY + GYd); % beta
                sigma2 = ((SYY+YYd) - B'*(SGG+GGd)*B)/n; % sigma
                fit.likelihood(i) = -(n*log(sigma2) + 2*k*p*log(delta) + sum(log(Rd.eig)) + sum(2*log(diag(S))) + n + n*log(2*pi))/2;
            end
        end
       	fit.sigma(i) = sqrt(sigma2);
        fit.beta(i,:) = B;
    end
end % figure; plot(fit.P,fit.likelihood)

function [obj, fit] = likelihood_at_interger(para, data) % likelihood at interger
    delta = para(1);     theta = para(2:end); % Initialize
    n = length(data.Y);    q = size(data.F,2);    
    likelihood = nan(length(data.P),1);
    sigma = nan(size(likelihood));
    beta = nan(length(data.P), q);
    for i = 1 : length(data.intP)
        p = data.intP(i);    k = floor(n/p);    p1 = n - p*k;
        Yb = data.Yb(1:p,i);    Ys = data.Ys(1:p1,i);
        Gb = data.Gb(1:p,:,i);    Gs = data.Gs(1:p1,:,i);
        GG = data.FF-Gs'*Gs;    YY = data.YY-Ys'*Ys;    GY = data.FY-Gs'*Ys;
        R0.c = data.corr(delta, theta, p, (1:p)', 1);    R0.eig = abs(fft(R0.c))+eps;
        Rd.c = k*R0.c/delta^2 + eye(p,1);    Rd.eig = real(fft(Rd.c));
        RRd.eig = (k/delta^2) * (R0.eig ./ Rd.eig);
        fftGb = fft(Gb);    fftYb = fft(Yb);
        Gbl = ifft(Rd.eig .\ fftGb);    Ybl = ifft(Rd.eig .\ fftYb);
        SGG = (GG+k*Gb'*(Gbl-Gb))/delta^2;
      	SGY = (GY+k*(Gbl-Gb)'*Yb)/delta^2;
     	SYY = (YY+k*Yb'*(Ybl-Yb))/delta^2;
        if p1==0
            B = SGG \ SGY; % beta
            sigma2 = (SYY - B'*SGG*B)/n; % sigma
            likelihood(p==data.P) = (n*log(sigma2*delta^2) + sum(log(Rd.eig)) + n + n*log(2*pi))/2;
        else
            RRdGb = ifft(RRd.eig .* fftGb);    RRdYb = ifft(RRd.eig .* fftYb);
            rpi = R0.c - ifft(RRd.eig .* R0.eig) + delta^2*eye(p,1);
            S = toep_chol(rpi(1:p1));
            Gd = Gs - RRdGb(1:p1,:); % Gamma dot
            Yd = Ys - RRdYb(1:p1); % Y dot
            Gdl = S \ Gd;    GGd = Gdl'*Gdl;
            Ydl = S \ Yd;    YYd = Ydl'*Ydl;
            GYd = Gdl'*Ydl;
            B = (SGG + GGd) \ (SGY + GYd); % beta
            sigma2 = ((SYY+YYd) - B'*(SGG+GGd)*B)/n; % sigma
            likelihood(p==data.P) = (n*log(sigma2) + 2*k*p*log(delta) + sum(log(Rd.eig)) + sum(2*log(diag(S))) + n + n*log(2*pi))/2;
        end
        sigma(p==data.P) = sqrt(sigma2);
        beta(p==data.P,:) = B;
    end % figure; plot(data.P(~isnan(likelihood)), -likelihood(~isnan(likelihood)), '-*')
    obj = min(likelihood);
    if  nargout > 1
      fit = struct('P', data.P, 'Y', data.Y, 'F',data.F, 'sigma',sigma, ...
          'beta',beta, 'corr',data.corr, 'regr',data.regr, 'likelihood',-likelihood);
    end
end

function  [t, f, fit, perf] = boxmin(objfunc, t0, lo, up, data)
%BOXMIN  Minimize with positive box constraints

    % Initialize
    [t, f, fit, itdata] = start(objfunc, t0, lo, up, data);
    if  ~isinf(f)
      % Iterate
      p = length(t);
      if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
      for  k = 1 : kmax
        th = t;
        [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data);
        [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data);
      end
    end
    perf = struct('nv',itdata.nv, 'perf',itdata.perf(:,1:itdata.nv));
end

function  [t, f, fit, itdata] = start(objfunc, t0, lo, up, data)
% Get starting point and iteration dataameters

    % Initialize
    t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
    D = 2 .^ ([1:p]'/(p+2));
    ee = find(up == lo);  % Equality constraints
    if  ~isempty(ee)
      D(ee) = ones(length(ee),1);   t(ee) = up(ee); 
    end
    ng = find(t < lo | up < t);  % Free starting values
    if  ~isempty(ng)
      t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
    end
    ne = find(D ~= 1);

    % Check starting point and initialize performance info
    [f  fit] = objfunc(t,data);   nv = 1;
    itdata = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
      'perf',zeros(p+2,200*p), 'nv',1);
    itdata.perf(:,1) = [t; f; 1];
    if  isinf(f)    % Bad dataameter region
      return
    end

    if  length(ng) > 1  % Try to improve starting guess
      d0 = 16;  d1 = 2;   q = length(ng);
      th = t;   fh = f;   jdom = ng(1);  
      for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
          tt = tk .* v; 
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 1];
          if  ff <= fk 
            tk = tt;  fk = ff;
            if  ff <= f
              t = tt;  f = ff;  fit = fitt; jdom = j;
            end
          else
            itdata.perf(end,nv) = -1;   break
          end
        end
      end % improve

      % Update Delta  
      if  jdom > 1
        D([1 jdom]) = D([jdom 1]); 
        itdata.D = D;
      end
    end % free variables

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = explore(objfunc, t, f, fit, itdata, data)
% Explore step

    nv = itdata.nv;   ne = itdata.ne;
    for  k = 1 : length(ne)
      j = ne(k);   tt = t;   DD = itdata.D(j);
      if  t(j) == itdata.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
      elseif  t(j) == itdata.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
      else
        atbd = 0;  tt(j) = min(itdata.up(j), t(j)*DD);
      end
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 2];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
      else
        itdata.perf(end,nv) = -2;
        if  ~atbd  % try decrease
          tt(j) = max(itdata.lo(j), t(j)/DD);
          [ff  fitt] = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 2];
          if  ff < f
            t = tt;  f = ff;  fit = fitt;
          else
            itdata.perf(end,nv) = -2;
          end
        end
      end
    end % k

    itdata.nv = nv;
end

function  [t, f, fit, itdata] = move(objfunc, th, t, f, fit, itdata, data)
% Pattern move

    nv = itdata.nv;   ne = itdata.ne;   p = length(t);
    v = t ./ th;
    if  all(v == 1)
      itdata.D = itdata.D([2:p 1]).^.2;
      return
    end

    % Proper move
    rept = 1;
    while  rept
      tt = min(itdata.up, max(itdata.lo, t .* v));  
      [ff  fitt] = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 3];
      if  ff < f
        t = tt;  f = ff;  fit = fitt;
        v = v .^ 2;
      else
        itdata.perf(end,nv) = -3;
        rept = 0;
      end
      if  any(tt == itdata.lo | tt == itdata.up), rept = 0; end
    end

    itdata.nv = nv;
    itdata.D = itdata.D([2:p 1]).^.25;
end