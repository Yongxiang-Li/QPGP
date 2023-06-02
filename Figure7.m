% Period estimation accuracy on sunspot data with varying SNR

close all; clear; clc;

addpath('routine');
load sunspot.dat
originSignal = sunspot(:,2);

S = 100;
for SNR = (-4:1:2)
    filename = ['Sunspot(SNR=',num2str(SNR),').mat'];
    if exist(filename, 'file')
        continue;
    end
    signals = [];    results = [];
    save(filename,'signals','results')
    parfor s = 1 : S
        % generate signal
        % rng(s) % For reproducibility
        noise = randn(size(originSignal));
        noise = noise*rms(originSignal)/rms(noise)/10^(SNR/20); % 20*log10(rms(originSignal)/rms(noise))
        noiseSignal = originSignal + noise;
        signal = struct('SNR',SNR,'noiseSignal',noiseSignal);
        % init analysis
        result = struct();
        period = (2:1:30)';
        result.Searchingperiod = period;
        % QPGP
        lob = [0.1 0.1 0.8];
        upb = [2 3 0.99];
        tic
        modelQPGP = fit_QPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
        modelQPGP.time = toc;
        % PGP
        lob = [0.1 0.1];
        upb = [2 3];
        tic
        modelCPGP = fit_CPGP(period, noiseSignal, @regpoly0, @period_sin_gauss_cov, lob, upb);
        modelCPGP.time = toc;
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
        modelFNLS.maxNoHarmonics = 10;
        modelFNLS.f0Bounds = [1/period(end), 1/period(1)];
        tic
        [modelFNLS] = FNLS( noiseSignal, length(noiseSignal), modelFNLS.maxNoHarmonics, modelFNLS.f0Bounds );
        modelFNLS.time = toc;
         % EPGP
        lob = [0.1 0.1 1];
        upb = [2 3 3];
        tic
        modelEPGP = fit_EPGP(period, noiseSignal, @regpoly0, @exponential_periodic_cov, lob, upb);
        EPGPtime = toc;
        modelEPGP.time = EPGPtime;

        result.modelEPGP = modelEPGP;
        result.modelNRCPE = modelNRCPE;
        result.modelQPGP = modelQPGP;
        result.modelCPGP = modelCPGP;
        result.modelMLPE = modelMLPE;
        result.modelFNLS = modelFNLS;
        results = [results; result];
        signals = [signals signal];
    end
    save(filename,'signals','results')
end

% Plot the results -- Figure 7 in the manuscript
count_NRC = zeros(20,1);
count_QPGP = zeros(20,1);
count_EPGP = zeros(20,1);
count_PGP = zeros(20,1);
count_MLPE = zeros(20,1);
count_FNLS = zeros(20,1);
SNR = -4:2;
periodRange = 11;
p = length(SNR);
period_NRC = zeros(100, p);
period_QPGP = zeros(100, p);
period_EPGP = zeros(100, p);
period_PGP = zeros(100, p);
period_MLPE = zeros(100, p);
period_FNLS = zeros(100, p);
for j = 1:p
    load(['Sunspot(SNR=',num2str(SNR(j)),').mat'])
    for i = 1:100
        period_NRC(i,j) = results(i).modelNRCPE.period;
        if  ismember (period_NRC(i,j),periodRange)
            count_NRC(j) = count_NRC(j) + 1;
        end
        period_QPGP(i,j) = results(i).modelQPGP.period;
        if ismember (period_QPGP(i,j),periodRange)
            count_QPGP(j) = count_QPGP(j) + 1;
        end
        period_EPGP(i,j) = results(i).modelEPGP.period;
        if ismember (period_EPGP(i,j),periodRange)
            count_EPGP(j) = count_EPGP(j) + 1;
        end 
        period_PGP(i,j) = results(i).modelCPGP.period;
        if ismember (period_PGP(i,j),periodRange)
            count_PGP(j) = count_PGP(j) + 1;
        end
        
        period_MLPE(i,j) = results(i).modelMLPE.period;
        if ismember (period_MLPE(i,j),periodRange)
            count_MLPE(j) = count_MLPE(j) + 1;
        end  
        period_FNLS(i,j) = results(i).modelFNLS.period;
        if ismember (period_FNLS(i,j),periodRange)
            count_FNLS(j) = count_FNLS(j) + 1;
        end  
    end
    
end

figure;
y=[count_PGP(1:p) count_QPGP(1:p) count_EPGP(1:p) count_FNLS(1:p) count_NRC(1:p) count_MLPE(1:p)];
b=bar(y);
grid on;
set(gca,'XTickLabel',{SNR})
ylim([0 100])
legend({'PGP','QPGP','EPGP','FNLS','NRC','MLPE'},'Location','NorthWest')
xlabel('SNR (dB)');
ylabel('Accuracy (%)');





