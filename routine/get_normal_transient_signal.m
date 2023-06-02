function [ signal, t ] = get_normal_transient_signal( N, f0, zeta0, T, fs, slip )   
    P = T*fs;
    t = (-3*P:3*P)';
    orginSig = exp(-(zeta0/sqrt(1-zeta0^2))*(2*pi*f0*(t/fs)).^2).*cos(2*pi*f0*(t/fs));
    signal = zeros(N,1);
    i = 0;
    index = max(i-2*P,1):min(i+2*P,N);
    signal(index) = signal(index) + orginSig(index-i+3*P+1);
    while i < N
        p = round((1+randn*slip)*T*fs);
        i = i + p;
        index = max(i-2*P,1):min(i+2*P,N);
        signal(index) = signal(index) + orginSig(index-i+3*P+1);
    end
    t = (1:N)/fs;
