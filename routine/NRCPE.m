function [period, Q] = NRCPE(signal, K, gamma)

    if nargin < 3, gamma = 0.5; end
    range = (1:K)';
    Q = nan(size(range));
    for j = range'
        Q(j) = nrc1(signal, j);
    end
    envlp = cummin(Q);
    [~,index] = max(Q);
    
    period = [];
    for g = gamma(:)'
        c = g*(Q(index)-Q(1))/(envlp(index)-envlp(1));
        [~,p] = max(Q - c*envlp);
        period = [period  p];
    end
end