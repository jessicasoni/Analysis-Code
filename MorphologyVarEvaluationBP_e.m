function [ energyinbandcollectionNTWD, energyinbandcollectionTWD, energyinbandcollectionTWDon2on, energyinbandcollectionTWDon2onmedian, medianwavecol ] = MorphologyVarEvaluationBP( bpsignal, ann, figsavepath )
%MorphologyVarEvaluationBP, Evaluate morphological variability in blood
%pressure waveform signal.
%   Detailed explanation goes here

% Recompute annotations for signal

N = length(bpsignal);
% isolate bp wave pks
bpwavespks = NaN(1, size(ann,2)-1);
for bpwaveidx = 1:length(bpwavespks)
    % find local max
    data = bpsignal(ann(bpwaveidx):ann(bpwaveidx+1));
    [pks,locs] = findpeaks(data);
    maxpk = max(pks); maxloc = locs(maxpk == pks);
    if (~isempty(maxloc))
        bpwavespks(bpwaveidx) = ann(bpwaveidx) - 1 + maxloc(1); % in case of duplicate
    else    % unable to find pk, assign as the mid point
        bpwavespks(bpwaveidx) = mean(ann(bpwaveidx), ann(bpwaveidx+1));
    end
%     % plot data and annotations
%     figure; plot(data); hold on;
%     plot(maxloc, data(maxloc), 'or', 'Color', 'b');
%     hold off;
%     close(gcf);
end

% confirm peaks
% % plot signal and annotations
% figure; hold on;
% plot(bpsignal);
% plot(ann, bpsignal(ann), 'or');
% plot(bpwavespks, bpsignal(bpwavespks), 'or', 'Color', 'b');
% hold off;

% compute median onset-pk and offset-pk interval
oninterv = median(bpwavespks - ann(1:end-1));
offinterv = median(ann(2:end) - bpwavespks);
% isolate each wave in the median interval around the peak
medianbpinterval = offinterv + oninterv;
noofwaves = length(bpwavespks) - 2;
%bpwaves = NaN(noofwaves, medianbpinterval+1);

% computation over 5 minute segments, may want to compute median interval
% over each 5 min segment
N5 = 5*60*240;
noofsegments = floor(N/N5);
energyinbandcollectionNTWD = NaN(1,noofsegments);
energyinbandcollectionTWD = NaN(1,noofsegments);
energyinbandcollectionTWDon2on = NaN(1,noofsegments);
energyinbandcollectionTWDon2onmedian = NaN(1,noofsegments);

% also store median bp wave for each seg
medianwavecol = NaN(noofsegments, 1000);

for segmentidx = 1:noofsegments % 55
    %segmentidx
    try
        
    segstart = 1 + (segmentidx-1)*N5;
    segend = segmentidx*N5;
    % find first bpwave onset after start idx
    bpwavestartidx = bpwavespks > segstart; % bpwavespks-oninterv > segstart;
    bpwavestartidx = find(bpwavestartidx);
    bpwavestartidx = bpwavespks(bpwavestartidx(1));
    % find last bpwave offset before end idx
    bpwaveendidx = bpwavespks < segend; % bpwavespks+offinterv < segend;
    bpwaveendidx = find(bpwaveendidx);
    bpwaveendidx = bpwavespks(bpwaveendidx(end));
    % find pks for cur window
    fiveminpks = bpwavespks(bpwavespks <= bpwaveendidx); fiveminpks = fiveminpks(fiveminpks >= bpwavestartidx);
    fiveminpks = fiveminpks - bpwavestartidx(1) + 1;    % 
    
    % remove one of two adjacent detected peaks less than 2 samples wide
    pk2remove = find(diff(fiveminpks) <= 2) + 1;
    fiveminpks(pk2remove) = [];
    
    % five minute signal
    fiveminsignal = bpsignal(bpwavestartidx:bpwaveendidx);
    abp = fiveminsignal;
    % annotations consistent with five min data. peak locations optimized.
    % recompute valleys
    
    % find the min point between two pks and assign as valley
    fiveminann = [];    % NaN(1,length(fiveminpks)-1);
    for pkidx = 1:length(fiveminpks)-1
        data = fiveminsignal(fiveminpks(pkidx):fiveminpks(pkidx+1));
        
        %if (length(data) <= 2)
        %    wait = 1;
        %end
        
        [val,locs] = findpeaks(-1*data);
%         figure; plot(-1*data);
%         close(gcf);
        minpk = max(val); minloc = locs(minpk == val);
        if (~isempty(minloc))
            fiveminann = [fiveminann (fiveminpks(pkidx) - 1 + minloc(1))]; % in case of duplicate
        else    % unable to find pk, assign as the mid point
            fiveminann = [fiveminann mean(fiveminpks(pkidx),fiveminpks(pkidx+1))];
        end
    end
    
    % plot the five min segment
%     figure; hold on;
%     plot(abp);
%     plot(fiveminann, fiveminsignal(fiveminann), 'or');
%     plot(fiveminpks, fiveminsignal(fiveminpks), 'or', 'Color', 'b');
%     hold off;
    
    % determine features for jsqi
    features = abpfeature(fiveminsignal,fiveminann');
    %disp('reach here abpfeature');
    % determine median interval over 5 minute segment (ignore last beat as it may not have offset)
    onintervseg = round(median(fiveminpks(2:end-1) - fiveminann(1:end-1)));
    offintervseg = round(median(fiveminann(2:end) - fiveminpks(2:end-1)));
    % isolate each wave in the median interval around the peak
    medianbpintervalseg = offintervseg + onintervseg;
    fiveminnoofwaves = length(fiveminpks(2:end-1));
    fiveminbpwaves = NaN(fiveminnoofwaves, medianbpintervalseg+1);
    fiveminbpwavesonset2onset = NaN(fiveminnoofwaves, (max(diff(fiveminann))+1));
    fiveminonset2pkampl = NaN(1,fiveminnoofwaves);
    
    % segment bp waves
    for pkidx = 1:fiveminnoofwaves   % start from 1 to end - 1 (last wave may not have an offset)
        fiveminbpwaves(pkidx,:) = abp(fiveminpks(pkidx+1)-onintervseg:fiveminpks(pkidx+1)+offintervseg);
        fiveminbpwavesonset2onset(pkidx,1:(fiveminann(pkidx+1)-fiveminann(pkidx))+1) = abp(fiveminann(pkidx):fiveminann(pkidx+1));
        % (pkidx+1) because second pk corresponds to first valley
        fiveminonset2pkampl(pkidx) = abp(fiveminpks(pkidx+1)) - abp(fiveminann(pkidx));
    end

    % use jsqi to determine clean windows
    [BeatQ, rjSQI] = jSQI(features, fiveminann, abp); %rjSQI
    % remove bad beats
    for bpwaveidx = 1:size(fiveminbpwaves,1)
        if (BeatQ(bpwaveidx,1) == 1)    % set bad beats to NaN
            fiveminbpwaves(bpwaveidx,:) = NaN;  % may cause error, double check
            fiveminbpwavesonset2onset(bpwaveidx,:) = NaN;
            fiveminonset2pkampl(bpwaveidx) = NaN;
        end
    end
    
    % remove nans
    originalwaves = size(fiveminbpwaves,1);
    fiveminbpwaves(~any(~isnan(fiveminbpwaves), 2),:)=[];
    fiveminbpwavesonset2onset(~any(~isnan(fiveminbpwavesonset2onset), 2),:)=[];
    fiveminonset2pkampl = fiveminonset2pkampl(~isnan(fiveminonset2pkampl));
    cleanwaves = size(fiveminbpwaves,1);
    % fraction of clean beats
    cleanfraction = (cleanwaves/originalwaves);
    % compute normalization factor
    normalizationfactor = median(fiveminonset2pkampl);
    
    % median wave
    if (size(fiveminbpwaves,1) > 1)
        medianwave = median(fiveminbpwaves);
        medianwavecol(segmentidx,1:length(medianwave)) = medianwave;
    else
        % signal unclean
        continue;
    end

    if (rjSQI > 0.9)
        
        % normalize
        fiveminbpwaves = fiveminbpwaves - median(fiveminbpwavesonset2onset(:,1));
        fiveminbpwavesonset2onset = fiveminbpwavesonset2onset - median(fiveminbpwavesonset2onset(:,1));
        fiveminbpwavesnorm = fiveminbpwaves ./ normalizationfactor;
        fiveminbpwavesonset2onsetnorm = fiveminbpwavesonset2onset ./ normalizationfactor;
        
        %plot(fiveminbpwavesnorm(1:2,:)'); plot(fiveminbpwavesonset2onsetnorm(1:2,:)');

        %  get sum of squared difference series
        % [NTWDseries TWDseries TWDseriesonset2onset TWDseriesonset2onsetmedian] = DiffSeriesCalcBP(fiveminbpwavesnorm(1:end,:), fiveminbpwavesonset2onsetnorm(1:end,:), medianwave ./ normalizationfactor);%, figsavepath, segmentidx);
        %[~, ~, TWDseriesonset2onset, TWDseriesonset2onsetmedian] = DiffSeriesCalcBP(fiveminbpwavesnorm(1:end,:), fiveminbpwavesonset2onsetnorm(1:end,:), medianwave ./ normalizationfactor);%, figsavepath, segmentidx);
        [~,~,TWDseriesonset2onset,TWDseriesonset2onsetmedian] = DiffSeriesCalcBP_e(fiveminbpwavesnorm(1:end,:), fiveminbpwavesonset2onsetnorm(1:end,:), medianwave ./ normalizationfactor);%, figsavepath, segmentidx);

        %% Compute energy in power spectrum (TWDseriesonset2onset)
        nfft = 256; % standardize the length of the fft to 256
        pxx = pwelch(TWDseriesonset2onset,[],[],nfft);
        N = nfft;
        stepsize = 0.5/length(pxx);
        
        dl = 2; du = 7;
        indexeverytwobeats = round(N/dl, 0);    % corresponding to every 2 beats
        indexeverysevenbeats = floor(N/du); % corresponding to every 7 beats
        
%         figure; plot((1:129)*stepsize, log10(pxx)); title(['pwelch 5 min (TWDon2on)']); hold on;
%         scatter(indexeverysevenbeats*stepsize, log10(pxx(indexeverysevenbeats)));
%         scatter(indexeverytwobeats*stepsize, log10(pxx(indexeverytwobeats)));
%         hold off;
        
        % compute energy in every 2-7 beats region
        energyinband = trapz((indexeverysevenbeats:indexeverytwobeats) * stepsize, pxx(indexeverysevenbeats:indexeverytwobeats));
        energyinbandcollectionTWDon2on(segmentidx) = energyinband;
%         figure(gcf); xlabel(['energy ' num2str(energyinband)]);
        
        %% Compute energy in power spectrum (TWDseriesonset2onsetmedian)
        nfft = 256; % standardize the length of the fft to 256
        pxx = pwelch(TWDseriesonset2onsetmedian,[],[],nfft);
        N = nfft;
        stepsize = 0.5/length(pxx);
        
        dl = 2; du = 7;
        indexeverytwobeats = round(N/dl, 0);    % corresponding to every 2 beats
        indexeverysevenbeats = floor(N/du); % corresponding to every 7 beats
        
%         figure; plot((1:129)*stepsize, log10(pxx)); title(['pwelch 5 min (TWDon2onmedian)']); hold on;
%         scatter(indexeverysevenbeats*stepsize, log10(pxx(indexeverysevenbeats)));
%         scatter(indexeverytwobeats*stepsize, log10(pxx(indexeverytwobeats)));
%         hold off;
        
        % compute energy in every 2-7 beats region
        energyinband = trapz((indexeverysevenbeats:indexeverytwobeats) * stepsize, pxx(indexeverysevenbeats:indexeverytwobeats));
        energyinbandcollectionTWDon2onmedian(segmentidx) = energyinband;
%         figure(gcf); xlabel(['energy ' num2str(energyinband)]);

    else
        disp(['unclean segment ' num2str(rjSQI)]);
        %energyinbandcollectionNTWD = [energyinbandcollectionNTWD NaN];
        %energyinbandcollectionTWD = [energyinbandcollectionTWD NaN];
        energyinbandcollectionTWDon2on(segmentidx) = NaN;
        energyinbandcollectionTWDon2onmedian(segmentidx) = NaN;
    end
    close all;
    
    catch
        disp('error evaluating segment');
        % if segment evaluation results in error put in NaN
        %energyinbandcollectionNTWD = [energyinbandcollectionNTWD NaN];
        %energyinbandcollectionTWD = [energyinbandcollectionTWD NaN];
        energyinbandcollectionTWDon2on(segmentidx) = NaN;
        energyinbandcollectionTWDon2onmedian(segmentidx) = NaN;
    end
    
end

end

