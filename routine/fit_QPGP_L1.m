function [ fit ] = fit_QPGP_L1( P, Y, regr, corr, lob, upb, theta0 )
% Periodic Gaussian Process
     
    data = struct('corr',corr, 'regr',regr, 'P', P, 'Y',Y);
    
    if nargin == 7
        [f, fit] = objfunc(theta0, data);
    else
       	% init dataameters
        u = linspace(0,1,6);    [x1, x2, x3] = meshgrid(u(2:end-1), u(2:end-1), u(2:end-1));
        thetas = lob + (upb-lob).*[x1(:),x2(:),x3(:)];
        objs = nan(size(thetas,1),1);
        fits = [];
        for i = 1 : size(thetas,1)
            [objs(i), fit] = objfunc(thetas(i,:), data);
            fits = [fits; fit]; % figure; plot(P, fit.likelihood);
        end % likelihoods = [fits(:).likelihood]; figure; plot(max(likelihoods'))
        [~, index] = min(objs);
        theta0 = thetas(index,:);
        [theta, f, fit, perf] = boxmin(theta0, lob, upb, data);
    end
    fit.theta0 = theta0;
    fit.thetahat = theta;
    [~, index] = max(fit.likelihood);
	fit.betahat = fit.beta(index,:);
    fit.sigmahat = fit.sigma(index);
    fit.period = fit.P(index);
    fit.Gamma = data.regr((1:fit.period)'); 
end

function [obj, fit] = objfunc(para, data) % using smt package
    delta = para(1);     theta = para(2);     omega = para(3);
    n = length(data.Y);    q = size(data.regr(0),2); 
    likelihood = nan(length(data.P),1);
    sigma = nan(size(likelihood));
    beta = nan(length(data.P), q);
    for i = 1 : length(data.P)
        p = data.P(i);    k = floor(n/p);
        Y = zeros(p,k);    Y(1:k*p) = data.Y(1:k*p);   
        Gamma = data.regr((1:p)'); 
        [So, Uo] = kmseig(k, omega);
        r0 = data.corr(delta, theta, p, (1:p)', 1);
        r0(1) = r0(1) + eps;    Sr = real(fft(r0));    Sr(Sr<0) = 0; 
        S = Sr * So' + delta^2;    Lambda = 1./S;
        UrGamma = fft(Gamma)/sqrt(p);    UoOnes = Uo'*ones(k,1);    UrYUo = fft(Y*Uo)/sqrt(p);
        ChiYY = real(sum(sum(conj(UrYUo).*(UrYUo.*Lambda))));
        ChiFF = real(UrGamma'*((Lambda*(UoOnes.^2)).*UrGamma));
        ChiFY = real(UrGamma'*(UrYUo.*Lambda)*UoOnes);
        
        Beta = ChiFF \ ChiFY; % beta
        sigma2 = (ChiYY - Beta'*ChiFF*Beta)/(p*k); % sigma
        likelihood(i) = ((p*k)*log(sigma2) + sum(sum(log(real(S)))) + (p*k) + (p*k)*log(2*pi))/2;
       
        sigma(i) = sqrt(sigma2);
        beta(i,:) = Beta;
    end % figure; plot(data.P, -real(likelihood))
    obj = min(likelihood);
    if  nargout > 1
      fit = struct('P', data.P, 'Y', data.Y, 'sigma',sigma, ...
          'beta',beta, 'corr',data.corr, 'regr',data.regr, 'likelihood',-likelihood);
    end
end

function  [t, f, fit, perf] = boxmin(t0, lo, up, data)
%BOXMIN  Minimize with positive box constraints

    % Initialize
    [t, f, fit, itdata] = start(t0, lo, up, data);
    if  ~isinf(f)
      % Iterate
      p = length(t);
      if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
      for  k = 1 : kmax
        th = t;
        [t, f, fit, itdata] = explore(t, f, fit, itdata, data);
        [t, f, fit, itdata] = move(th, t, f, fit, itdata, data);
      end
    end
    perf = struct('nv',itdata.nv, 'perf',itdata.perf(:,1:itdata.nv));
end

function  [t, f, fit, itdata] = start(t0, lo, up, data)
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

function  [t, f, fit, itdata] = explore(t, f, fit, itdata, data)
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

function  [t, f, fit, itdata] = move(th, t, f, fit, itdata, data)
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