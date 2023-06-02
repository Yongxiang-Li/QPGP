function [ fit ] = fit_PGP( P, Y, regr, corr, lob, upb )
% Toeplitz-based Periodic Gaussian Process
    F = regr((1:length(Y))'/length(Y));
    data = struct('corr',corr, 'regr',regr, 'F',F, 'P', P, 'Y',Y);
    
    % init dataameters
    u = linspace(0,1,6);    [x1, x2] = meshgrid(u(2:end-1));
    thetas = lob + (upb-lob).*[x1(:),x2(:)];
    objs = nan(size(thetas,1),1);
    fits = [];
    for i = 1 : size(thetas,1)
        [objs(i), fit] = objfunc(thetas(i,:), data);
        fits = [fits; fit]; % figure; plot(fit.likelihood);
    end % likelihoods = [fits(:).likelihood]; figure; plot(max(likelihoods'))
    [~, index] = min(objs);
    theta0 = thetas(index,:);
    [theta, f, fit, perf] = boxmin(theta0, lob, upb, data);

    fit.theta0 = theta0;
    fit.thetahat = theta;
    [~, index] = max(fit.likelihood);
	fit.betahat = fit.beta(index,:);
    fit.sigmahat = fit.sigma(index);
    fit.period = fit.P(index);
end


function [obj, fit] = objfunc(para, data)
    delta = para(1);     theta = para(2:end); % Initialize
    n = length(data.Y);    q = size(data.F,2);    
    likelihood = nan(length(data.P),1);  costfun = likelihood;
    sigma = nan(size(likelihood));
    beta = nan(length(data.P), q);
    for i = 1 : length(data.P)
        p = data.P(i)';
        R = data.corr(delta, theta, p, (1:n)', (1:n)');    R = R + delta^2*eye(n);
        L = chol(R, 'lower');
        Ft = L \ data.F;    Yt = L \ data.Y; 
        [Q G] = qr(Ft,0);    B = G \ (Q'*Yt);
        rho = Yt - Ft*B;     sigma2 = sum(rho.^2)/n;
        likelihood(i) = (n*log(sigma2) + 2*sum(log(diag(L))) + n + n*log(2*pi))/2; 
        sigma(i) = sqrt(sigma2);
        beta(i,:) = B;
    end % figure; plot(-likelihood)
    obj = min(likelihood);
    if  nargout > 1
      fit = struct('P', data.P, 'Y', data.Y, 'F',data.F, 'sigma',sigma, ...
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