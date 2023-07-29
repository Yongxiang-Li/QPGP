% Simulation study

close all;
clear;
clc;

addpath('routine', 'routine/fastF0Nls');

% Vibration Signal
f0 = 7;
zeta0 = 0.02;
T0 = 1; % true period
S = 100;    % number of simulations
SNR = -16;
delta = 10^(-SNR/20);

for N = 1000:1000:8000
    filename = ['N_VIBRATION(SNR=',num2str(SNR),';N=',num2str(N),').mat'];
    if exist(filename, 'file')
        continue;
    end
    signals = [];    results = [];
    save(filename,'signals','results')
    for s = 1 : S
        % generate signal
        fs = 80+unidrnd(120);
        [originSignal, t] = get_normal_transient_signal(N, f0, zeta0, T0, fs, 0.01); % get signal
        osc = rand(N,1);
        originSignal = osc.*originSignal;
        originSignal = originSignal/std(originSignal);
        noise = randn(size(originSignal));
        noise = noise*rms(originSignal)/rms(noise)/10^(SNR/20); % 20*log10(rms(originSignal)/rms(noise))
        noiseSignal = originSignal + noise;
        signal = struct('SNR',SNR,'N',N,'noiseSignal',noiseSignal,'orginSignal',originSignal,'period',fs*T0);
        
        % init analysis
        result = struct();
        period = (10:1:300)';
        result.Searchingperiod = period;
        % QPGP
        lob = [0.1 0.1 0.8];    upb = [10 15 0.999999];
        tic
        modelQPGP = fit_QPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
        modelQPGP.time = toc;
        % EPGP
        tic
        lob = [0.1 0.1 0.1];      upb = [10 15 10];
        modelEPGP = fit_EPGP(period, noiseSignal, @regpoly0, @exponential_periodic_cov, lob, upb);
        EPGPtime = toc;
        modelEPGP.time = EPGPtime;
        % CPGP
        lob = [0.1 0.1];
        upb = [10 15];
        tic
        modelCPGP = fit_CPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
        modelCPGP.time = toc;
        % TPGP
        tic
        modelTPGP = fit_TPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
        TPGPtime = toc;
        modelTPGP.time = TPGPtime;
        % MLPE
        modelMLPE = struct();
        Mlpe = nan(size(period));
        tic
        for p = 1 : length(period)
            [Mlpe(p)] = mlpe( noiseSignal, period(p));
        end
        modelMLPE.likelihood = Mlpe;
        [~, modelMLPE.period] = max(Mlpe);
        modelMLPE.time = toc;
        % NRCPE
        modelNRCPE = struct();
        tic
        [NRCperiod, nrcQ] = NRCPE(noiseSignal, length(period), 0.5);
        modelNRCPE.period = NRCperiod;
        modelNRCPE.nrcQ = nrcQ;
        modelNRCPE.time = toc;
        % FNLS
        modelFNLS = struct();
        modelFNLS.maxNoHarmonics = 20;
        modelFNLS.f0Bounds = [1/period(end), 1/period(1)];
        tic
        [modelFNLS] = FNLS( noiseSignal, N, modelFNLS.maxNoHarmonics, modelFNLS.f0Bounds );
        modelFNLS.time = toc;

        % save results
        result.modelQPGP = modelQPGP;
        result.modelEPGP = modelEPGP;
        result.modelTPGP = modelTPGP;
        result.modelCPGP = modelCPGP;
        result.modelFNLS = modelFNLS;
        result.modelNRCPE = modelNRCPE;
        result.modelMLPE = modelMLPE;
        results = [results; result];
        signals = [signals signal];
    end
    save(filename,'signals','results')
end

count_NRC = zeros(20,1);
count_QPGP = zeros(20,1);
count_EPGP = zeros(20,1);
count_PGP = zeros(20,1);
count_MLPE = zeros(20,1);
count_FNLS = zeros(20,1);
SNR = -16;
N = 1000:1000:8000;
periodRecord = zeros(100,1);
p = length(N);
T_NRC = zeros(p,1);
T_QPGP = zeros(p,1);
T_EPGP = zeros(p,1);
T_PGP = zeros(p,1);
T_MLPE = zeros(p,1);
T_FNLS = zeros(p,1);
for j = 1:p
    load(['N_VIBRATION(SNR=',num2str(SNR),';N=',num2str(N(j)),').mat'])
    period_NRC = zeros(100, 1);
    period_QPGP = zeros(100, 1);
    period_EPGP = zeros(100, 1);
    period_PGP = zeros(100, 1);
    period_MLPE = zeros(100, 1);
    period_FNLS = zeros(100, 1);
    
    time_NRC = zeros(100, 1);
    time_QPGP = zeros(100, 1);
    time_EPGP = zeros(100, 1);
    time_PGP = zeros(100, 1);
    time_MLPE = zeros(100, 1);
    time_FNLS = zeros(100, 1);
    for i = 1:100
        periodRecord(i) = signals(i).period;
        periodRange = signals(i).period;
        period_NRC(i) = results(i).modelNRCPE.period;
        time_NRC(i) = results(i).modelNRCPE.time;
        if  ismember (period_NRC(i),periodRange)
            count_NRC(j) = count_NRC(j) + 1;
        end
        period_QPGP(i) = results(i).modelQPGP.period;
        time_QPGP(i) = results(i).modelQPGP.time;
        if ismember (period_QPGP(i),periodRange)
            count_QPGP(j) = count_QPGP(j) + 1;
        end
           
        period_EPGP(i) = results(i).modelEPGP.period;
        time_EPGP(i) = results(i).modelEPGP.time;
        if ismember (period_EPGP(i),periodRange)
            count_EPGP(j) = count_EPGP(j) + 1;
        end
        
        period_PGP(i) = results(i).modelCPGP.period;
        time_PGP(i) = results(i).modelCPGP.time;
        if ismember (period_PGP(i),periodRange)
            count_PGP(j) = count_PGP(j) + 1;
        end
        period_MLPE(i) = results(i).modelMLPE.period;
        time_MLPE(i) = results(i).modelMLPE.time;
        if ismember (period_MLPE(i),periodRange)
            count_MLPE(j) = count_MLPE(j) + 1;
        end  
        period_FNLS(i) = results(i).modelFNLS.period;
        time_FNLS(i) = results(i).modelFNLS.time;
        if ismember (period_FNLS(i),periodRange)
            count_FNLS(j) = count_FNLS(j) + 1;
        end  
    end   
    T_NRC(j) = mean(time_NRC);
    T_QPGP(j) = mean(time_QPGP);
    T_EPGP(j) = mean(time_EPGP);
    T_PGP(j) = mean(time_PGP);
    T_MLPE(j) = mean(time_MLPE);
    T_FNLS(j) = mean(time_FNLS);
end

figure;
y=[count_PGP(1:p) count_QPGP(1:p) count_EPGP(1:p) count_FNLS(1:p) count_NRC(1:p) count_MLPE(1:p)];
b=bar(y);
grid on
set(gca,'XTickLabel',{N})
legend({'TPGP(CPGP)','QPGP','EPGP','FNLS','NRC','MLPE'},'Location','NorthWest')
xlabel('N');
ylabel('Accuracy (%)');