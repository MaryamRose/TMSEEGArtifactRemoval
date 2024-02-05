% ===================================================================================================================
% The function to detect blink artifacts from TMS-EEG data
% 
% blinkPulse = unBlinkWithFDCollective(xTrds, M, clustNum, fdThreshold, blinkAmpThr, dsRatio,...
%                                     lagThreshold, reducedSR, conservativeLen, showClustersPic)
%
%       blinkPulse: a binary array (0/1) where 1 represents the existence
%                   of an artifact and 0 shows lack of it
%       xTrds: the 1-dimensional array of a trial in a channel (size: 1 x time_points)
%       M: Window length in embedding step 
%       clustNum: number of clusters
%       fdThreshold: the larger 'fdThreshold' value, more likely to remove slow waves (blinks) - 
%                   a very small value could miss blinks, a very large value over-corrects
%       blinkAmpThr: a blink should be larger than blinkAmpThr uV to be
%                    detected
%       dsRatio: down-sampling ratio
%       lagThreshold: attach two pulses together if the spacing is less than lagThreshold in ms
%       reducedSR: the new reduced sampling rate (e.g. 250 Hz)
%       conservativeLen: the number of samples to add to the strat and end of every piece of pulses 
%       showClustersPic: Flag to check if the right cluster number was chosen
%
% This function comes as part of the package for cleaning TMS-EEG datasets with the proposed method in: 
% "A Novel Approach to Artifact Removal in TMS-EEG Data Using Curve Fitting and Wavelet-Based Estimation",
% Maryam Rostami, Reza Zomorrodi, Reza Rostami, Gholam-Ali Hosseinzadeh
%
% Author of the codes: Maryam Rostami [mar.rostami@ut.ac.ir]
%
% The codes are shared publicly based on the requirements of the journal
% "Brain Topography" and are open-source with no waranty for the results.
% Please avoid using these codes for clinical and diagnostics applications.
% ====================================================================================================================

function blinkPulse = unBlinkWithFDCollective(xTrds, M, clustNum, fdThreshold, blinkAmpThr, dsRatio, lagThreshold, reducedSR, conservativeLen,...
    showClustersPic)
% The trial/channel-wise blink detection is adapted from the paper: 
% "Maddirala AK, Veluvolu KC. Eye ‑ blink artifact removal from single channel EEG with k ‑ means and SSA 2021:1–14."

    
    %%% embedding
    N = length(xTrds);
    K = N-M+1;
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    tags = {FigList.Tag};
    if showClustersPic
        if ismember('QC_fd', tags)
            figure(FigList(ismember('QC_fd', tags)))
        else
            figure('units','normal','position',[0 0 1 1],'tag','QC_fd'); 
        end
    end
    blinkPulse = zeros(1,N*dsRatio);
    blinkPulseTimes = zeros(size(xTrds,1),N*dsRatio);
    for iCh = 1:size(xTrds,1)
        
        xTrdsArr = xTrds(iCh,:); 
        xTrdsArr = detrend(xTrdsArr,1);
        
        % first filter the data to facilitate the slow blink detection
        [b,a] = fir1(20, [0.5 80]/(reducedSR/2), 'bandpass');
        [b1,a1] = butter(3, [45 55]/(reducedSR/2), 'stop');
        xTrf1 = filtfilt(b,a,xTrdsArr);
        xTrf2 = filtfilt(b1,a1,xTrf1);
        
        xTemp = arrayfun(@(i) xTrf2(i:i+M-1)', 1:K, 'un', 0);
        xEmbed = [xTemp{:}];

        %%% feature extraction
        f1 = sum(xEmbed.^2); % energy of the signal
        f2 = sqrt(var(diff(xEmbed,[],1)) ./ var(xEmbed)); % Hjorth mobility
        f3 = kurtosis(xEmbed); % kurtosis
        f4 = max(xEmbed,[],1)-abs(min(xEmbed,[],1)); % dynamic range
        featuresMat = vertcat(f1,f2,f3,f4);

        %%% k-means clustering  
        rng('default'); % for setting the seed of the random generator in kmeans as zero
        fIdx = kmeans(featuresMat', clustNum);
        xBar = zeros([size(xEmbed), clustNum]);
        for iClust = 1:clustNum
            xBar(:,:,iClust) = xEmbed;
            lenTemp = sum(fIdx~=iClust);
            xBar(:, fIdx~=iClust, iClust) = zeros(M,lenTemp);
        end

        %%% diagonal averaging
        sBars = zeros(clustNum,N);
        for iSbar=1:clustNum
            xTemp = squeeze(xBar(:,:,iSbar));  
            xDAvg1 = arrayfun(@(n) 1/n.*trace(xTemp(1:n, n-(1:n)+1)), 1:M, 'un', 0);
            xDAvg2 = arrayfun(@(n) 1/M.*trace(xTemp(1:M,n-(1:M)+1)), 1+M:K, 'un', 0);            
            xDAvg3 = arrayfun(@(n) 1/(N-n+1).*trace(xTemp(n-K+1:N-K+1,n-(n-K+1:N-K+1)+1)), K+1:N, 'un', 0);
            xDAvg = horzcat(xDAvg1, xDAvg2, xDAvg3);
            sBars(iSbar,:) = [xDAvg{:}];
        end

        KFD = zeros(1,clustNum);
        for iSbar = 1:clustNum                
            KFD(iSbar) = Katz_FD(sBars(iSbar,:));
        end

        if showClustersPic       
            subplot(clustNum+2,size(xTrds,1),iCh); plot(xTrdsArr); title('signal')
            subplot(clustNum+2,size(xTrds,1),iCh+size(xTrds,1)); plot(xTrf2);
            title('filtered signal')
            for i=1:clustNum
                subplot(clustNum+2,size(xTrds,1),iCh+(i+1)*size(xTrds,1)); plot(sBars(i,:));
                title(num2str(KFD(i)));
            end
        end
        
        sBarsMax = max(abs(sBars),[],2)';
        sBarIndsToRemove = KFD<=fdThreshold & sBarsMax>=blinkAmpThr; 
        aBarAr = sBars(sBarIndsToRemove,:);        
        if sum(sBarIndsToRemove)~=0 

            aBarAr = sum(aBarAr,1);
            aBarB = zeros(size(aBarAr)); % The 0-1 square pulse
            aBarB(abs(aBarAr)>0) = 1;   
            aBarB_us = interp1(1:length(aBarB), aBarB, linspace(1,length(aBarB),length(aBarB)*dsRatio));
            blinkPulseTimes(iCh,:) = aBarB_us;
            blinkPulseTimes(iCh,blinkPulseTimes(iCh,:)~=0) = 1;
            
            
        end
        
    end
    
    % Select the infomative channels and combine them
    voters = ~ismember(sum(blinkPulseTimes,2), [0 N*dsRatio]); % channels which offer a pulse not a full zeros or ones
    voters2Ones = ismember(sum(blinkPulseTimes,2), [N*dsRatio]); 
    if sum(voters2Ones)>=length(voters2Ones)/2
        blinkPulse = ones(1,N*dsRatio);
    elseif sum(voters)>=3 % then there is a blink
        pulsesTemp = blinkPulseTimes(voters,:);
        pulseDiff = diff([zeros(size(pulsesTemp,1),1) pulsesTemp zeros(size(pulsesTemp,1),1)],1,2);
        [tStarts, tEnds] = deal([]);
        for ip = 1:size(pulsesTemp,1)
            if length(find(pulseDiff(ip,:)==1)) > 1
                tStarts(ip,1:length(find(pulseDiff(ip,:)==1))) = find(pulseDiff(ip,:)==1);
                tEnds(ip,1:length(find(pulseDiff(ip,:)==1))) = find(pulseDiff(ip,:)==-1)-1;
            else
                tStarts(ip,1) = find(pulseDiff(ip,:)==1);
                tEnds(ip,1) = find(pulseDiff(ip,:)==-1)-1;
            end
        end
        
        numPulses = zeros(1,size(pulsesTemp,1));
        for iip = 1:size(pulsesTemp,1)
            [pulsesTemp(iip,:), tStartsNew, ~] = ...
                fillPulseGaps(tStarts(iip,1:sum(tStarts(iip,:)~=0)), tEnds(iip,1:sum(tStarts(iip,:)~=0)), pulsesTemp(iip,:), lagThreshold);
            numPulses(iip) = sum(tStartsNew~=0);
        end
                
        % the number of blink in final pulse is the mode of blink numbers in selected channel 
        pulses = pulsesTemp(numPulses == mode(numPulses),:);

        % Combine the pulses of selected channels by OR
        blinkPulse = pulses(1,:); %cell array containing arrays 
        for iand=2:size(pulses,1)
            blinkPulse = or(blinkPulse,pulses(iand,:));
        end
        bpdiff = diff([0 blinkPulse 0]);
        tStarts = find(bpdiff==1);
        tEnds = find(bpdiff==-1)-1;
        intervallength = tEnds - tStarts;
        
        % add the conservative length of ones before and after every rectangles that are larger than 100 ms (min length of a blink)
        for iiInd=1:length(tStarts)
            if intervallength(iiInd) > 100e-3*dsRatio*reducedSR
                if tStarts(iiInd) < conservativeLen || tEnds(iiInd) + conservativeLen > length(blinkPulse)
                    if tStarts(iiInd) < conservativeLen
                        blinkPulse(1:min([tEnds(iiInd)+conservativeLen, length(blinkPulse)])) = 1;
                    end
                    if tEnds(iiInd) + conservativeLen > length(blinkPulse)
                        blinkPulse(max([1,tStarts(iiInd)-conservativeLen]):length(blinkPulse)) = 1;
                    end

                else
                    blinkPulse(max([1,tStarts(iiInd)-conservativeLen]):min([tEnds(iiInd)+conservativeLen,length(blinkPulse)])) = 1;
                end
            end
        end
        bpdiff = diff([0 blinkPulse 0]);
        tStarts = find(bpdiff==1);
        tEnds = find(bpdiff==-1)-1;
        
        
        % fill the gaps that are shorter than 100 ms. between every two rectangle 
        [blinkPulse, ~, ~] = fillPulseGaps(tStarts, tEnds, blinkPulse, lagThreshold);
               
        
    end
  
end



function [blinkPulse, tStarts, tEnds] = fillPulseGaps(tStarts, tEnds, blinkPulse, lagThreshold)
% attach two pulses if they are so close
        emptyLags = tStarts(2:end) - tEnds(1:end-1);
        for ilags = 1:length(emptyLags)
            if emptyLags(ilags) < lagThreshold
                blinkPulse(tEnds(ilags):tStarts(ilags+1)) = 1;
            end
        end
        
        bpdiff = diff([0 blinkPulse 0]);
        tStarts = find(bpdiff==1);
        tEnds = find(bpdiff==-1)-1;
end