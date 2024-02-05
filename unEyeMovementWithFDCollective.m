% ===================================================================================================================
% The function to detect eye movement artifacts from TMS-EEG data
% 
% movePulse = unEyeMovementWithFDCollective(twoFrLatData,eyeMoveAmpThr,twinLen,minMoveLen,conservativeLen,lagThreshold,dsRatio)
%
%       movePulse: a binary array (0/1) where 1 represents the existence
%                   of an artifact and 0 shows lack of it
%       twoFrLatData: the matrix for two fronto-lateral channels (e.g. F7 and F8) (size: 2 x time points)
%       eyeMoveAmpThr: multiplier of F7/F8 channels should be larger than this to be detected as eye movement
%       twinLen: windows of this length are compared between F7/F8 to detect the eye movement pulse
%       minMoveLen: the minimum length of an eye movement artifact in samples
%       conservativeLen: the number of samples to add to the strat and end of every piece of pulses 
%       lagThreshold: attach two pulses together if the spacing is less than lagThreshold in ms
%       dsRatio: down-sampling ratio
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

function movePulse = unEyeMovementWithFDCollective(twoFrLatData,eyeMoveAmpThr,twinLen,minMoveLen,conservativeLen,lagThreshold,dsRatio)


    f7f8Vals = arrayfun(@(x) mean(twoFrLatData(:,x:x+twinLen),2), 1:length(twoFrLatData)-twinLen, 'un',0);
    f7f8Vals = [f7f8Vals{:}];
    f7f8Prods = prod(f7f8Vals,1);

    negFlags = f7f8Prods<eyeMoveAmpThr;
    movePulse_ds = zeros(1,length(twoFrLatData));
    movePulse_ds(negFlags) = 1;
    if f7f8Prods(end)<0
        movePulse_ds(end-twinLen+1:end) = ones(1,twinLen);
    else
        movePulse_ds(end-twinLen+1:end) = zeros(1,twinLen);
    end

    % attach two pulses together if the spacing is less than 100 ms
    pulseDiff = diff([0 movePulse_ds 0]);
    tStarts = find(pulseDiff==1);
    tEnds = find(pulseDiff==-1)-1;
    if sum(movePulse_ds)~=0 % skip if there is no eye movement 
        [movePulse_ds, ~, ~] = fillPulseGaps(tStarts, tEnds, movePulse_ds, lagThreshold);
    end

    movePulsedd = diff([0 movePulse_ds 0]);
    startAndEndFlags = [find(movePulsedd==1); find(movePulsedd==-1)-1];
    flagLens = startAndEndFlags(2,:) - startAndEndFlags(1,:);
    noisyInds = find(flagLens<minMoveLen);
    for iIndF = 1:length(noisyInds)
        movePulse_ds(startAndEndFlags(1,noisyInds(iIndF)):startAndEndFlags(2,noisyInds(iIndF))) = 0;
    end

    movePulse = interp1(1:length(movePulse_ds), movePulse_ds, linspace(1,length(movePulse_ds),length(movePulse_ds)*dsRatio));   
    movePulse(movePulse~=0) = 1;

    if sum(movePulse)~=0 % skip if there is no eye movement

        pulseDiff = diff([0 movePulse 0]);
        tStarts = find(pulseDiff==1);
        tEnds = find(pulseDiff==-1)-1;

        for iiInd=1:length(tStarts)
            if tStarts(iiInd) < conservativeLen || tEnds(iiInd) + conservativeLen > length(movePulse)
                if tStarts(iiInd) < conservativeLen
                    movePulse(1:min([tEnds(iiInd)+conservativeLen, length(movePulse)])) = 1;
                end
                if tEnds(iiInd) + conservativeLen > length(movePulse)
                    movePulse(max([1,tStarts(iiInd)-conservativeLen]):length(movePulse)) = 1;
                end

            else
                movePulse(max([1,tStarts(iiInd)-conservativeLen]):min([tEnds(iiInd)+conservativeLen,length(movePulse)])) = 1;
            end
        end

        pulseDiff = diff([0 movePulse 0]);
        tStarts = find(pulseDiff==1);
        tEnds = find(pulseDiff==-1)-1;
        [movePulse, tStarts, tEnds] = fillPulseGaps(tStarts, tEnds, movePulse, lagThreshold);

        if tStarts(1) < minMoveLen*dsRatio
            movePulse(1:tStarts(1)) = 1;
        end
        if length(movePulse) - tEnds(end) < minMoveLen*dsRatio
            movePulse(tEnds(end):end) = 1;
        end

        [movePulse, ~, ~] = fillPulseGaps(tStarts, tEnds, movePulse, lagThreshold);

    end

    
    
   
end




function [movePulse_us, tStarts, tEnds] = fillPulseGaps(tStarts, tEnds, movePulse_us, lagThreshold)
% attach two pulses if they are so close
        emptyLags = tStarts(2:end) - tEnds(1:end-1);
        for ilags = 1:length(emptyLags)
            if emptyLags(ilags) < lagThreshold
                movePulse_us(tEnds(ilags):tStarts(ilags+1)) = 1;
            end
        end
        
        bpdiff = diff([0 movePulse_us 0]);
        tStarts = find(bpdiff==1);
        tEnds = find(bpdiff==-1)-1;
end