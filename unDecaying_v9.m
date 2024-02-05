% ===================================================================================================================
% The function to fit and remove the decay (discharge) artifact from TMS-EEG data
% 
% [cleanSignal, iniExpInd, endExpInd, f0] = unDecaying_v9(signal, timeToZP, ...
%    baselineCorrTime, fsample, time, zeroSample, timeToFit, tolTime, nChngTime,...
%    lowerLim, upperLim,startPoint,focusMsecs)
%
%       cleanSignal: the signal which is cleaned from decay artifacts
%       iniExpInd: the index for the start of the decay fitting
%       endExpInd: the index for the end of the decay fitting
%       f0: the fit model 
%       signal: the 1-dimensional array of a trial in a channel (size: 1 x time_points)
%       timeToZP: the TMS pulse time range that was interpolated in seconds
%       baselineCorrTime: the time range to do baseline correction for
%                        'signal', in seconds
%       fsample: sampling rate (Hz)
%       time: time array in sec
%       zeroSample: the index of the sample representing 0 in time
%       timeToFit: the time range to fit the decay, in case of manual
%                   estimation of the start and ending time points of the decay;
%                   otherwise leave it as empty to find the optimal time
%       tolTime: tolerance time to estimate the starting time point of the
%                decay fit
%       nChngTime: No change time point, representing the time after which
%                   no TMS-evoked activity is expected
%       lowerLim: lower limitation for the fit function while estimating
%                  the exponential parameters
%       upperLim: upper limitation for the fit function while estimating
%                  the exponential parameters
%       startPoint: starting point for the fit function estimations
%       focusMsecs: the time length after the zero-padded range, in which
%                   large weights are assigned to capture the main part of the decay
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

function [cleanSignal, iniExpInd, endExpInd, f0] = unDecaying_v9(signal, timeToZP, ...
    baselineCorrTime, fsample, time, zeroSample, timeToFit, tolTime, nChngTime, lowerLim, upperLim,startPoint,focusMsecs)


tCut1 = dsearchn(time', timeToZP(1));
tCut2 = dsearchn(time', timeToZP(2));
ind500 = dsearchn(time', nChngTime);
indEnd = length(time);


tolTimeSample = tolTime * fsample;


signal(1:tCut1) = detrend(signal(1:tCut1),1); % detrend the baseline to avoid drops/jumps after tCut2
signal(tCut2:end) = detrend(signal(tCut2:end),1);
signal = signal - mean(signal(1+round(zeroSample+baselineCorrTime(1)*fsample):round(zeroSample+baselineCorrTime(2)*fsample)-1)); % baseline correction
signal(tCut1:tCut2) = signal(tCut2);


% >>>>>> exponential fit model <<<<<<<<<
g = fittype('a+b*exp(c*x)'); % for fit function
cleanSignal = signal;
f0 = [];

signalSmoothed = smooth(signal(tCut2:end),20e-3*fsample);
if isempty(timeToFit)
    
    % initial point
    x11 = tCut1; 
    y11 = signal(tCut1);
    % ending point
    x22 = tCut2+tolTimeSample; 
    y22 = signal(tCut2+tolTimeSample);
    
    % use triangle method to find the start of the exp function (the initial exp index):
    % define point to line distance : https://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    pointsToLineDist = arrayfun(@(i) abs((x22-x11)*(y11-signal(i)) - (x11-i)*(y22-y11)) / sqrt((x22-x11)^2+(y22-y11)^2) ,x11:x22,'un',0);
    pointsToLineDist = [pointsToLineDist{:}];
    pointsAboveLine = arrayfun(@(i) (signal(i)-y11)/(i-x11) > (y22-y11)/(x22-x11), x11:x22, 'un', 0);
    pointsAboveLine = [pointsAboveLine{:}];

    endExpIndTemp = tCut2+50e-3*fsample; % this is only used for estimating the sign of exponential (rising/falling) and it's the first 50 msec. after tCut2
    slopeEstAvg = mean(diff(signalSmoothed(1:endExpIndTemp-tCut2+1)));
    if slopeEstAvg <= 0
        [~, maxDistIndTemp] = max(pointsToLineDist(pointsAboveLine));
        tempInds = find(pointsAboveLine);
    else
        pointsBelowLine = ~pointsAboveLine;
        [~, maxDistIndTemp] = max(pointsToLineDist(pointsBelowLine));
        tempInds = find(pointsBelowLine);
    end
    maxDistInd = tempInds(maxDistIndTemp);
    iniExpInd = maxDistInd-1+x11;
    if iniExpInd < tCut2
        iniExpInd = tCut2;
    end

    endExpInd = length(time); 
    blLevel = mean(signal(ind500:indEnd));    
    
else
    iniExpInd = timeToFit(1)*fsample + zeroSample;
    endExpInd = timeToFit(2)*fsample + zeroSample;  
end

%%%%%%%%%%%%%%%

if  ~isempty(iniExpInd) && ~isempty(endExpInd)  % do the fit if you found indices for the start and end of the exponential decays
    
    % the "real" decay should be at least 5 msec. and fit start should be earlier than tCut2+tolTimeSample
    if (endExpInd-iniExpInd > 5e-3*fsample) && (iniExpInd <= tCut2+tolTimeSample) 

        yDiff_raw = abs(diff(signal([iniExpInd,endExpInd])));
        yDiff_smooth = abs(diff(signalSmoothed([iniExpInd-tCut2+1,endExpInd-tCut2+1])));
        if ~(yDiff_smooth<5) && ~(yDiff_raw<5) % if only the change in y axis is larger than 5 uV, then some exp. decay might be seen, otherwise it's clean
                lenTemp = ind500-1-iniExpInd+1; 
                if slopeEstAvg>0
                    f0 = fit([time(iniExpInd:ind500-1) time(ind500:indEnd)]', [signal(iniExpInd:ind500-1) blLevel*ones(1,(indEnd-ind500+1))]', g, ...
                        'TolFun',1e-4, 'TolX',1e-4,'StartPoint', [1 startPoint.b startPoint.c],...
                        'Weights',[1e10 1e4*ones(1,focusMsecs*fsample-1) ones(1,lenTemp-focusMsecs*fsample-1) 1e4 1e4*ones(1,indEnd-ind500+1)],...
                        'lower', lowerLim, 'upper', upperLim); 
                else
                    f0 = fit([time(iniExpInd:ind500-1) time(ind500:indEnd)]', [signal(iniExpInd:ind500-1) blLevel*ones(1,(indEnd-ind500+1))]', g, ...
                        'TolFun',1e-4, 'TolX',1e-4,'StartPoint', [1 startPoint.b startPoint.c],... 
                        'Weights',[1e10 1e4*ones(1,focusMsecs*fsample-1) ones(1,lenTemp-focusMsecs*fsample-1) 1e4 1e4*ones(1,indEnd-ind500+1)],...
                        'lower', lowerLim, 'upper', upperLim);           
                end

                zeroReach = endExpInd;

                if zeroReach > iniExpInd
                    yFitted = f0(time(iniExpInd:zeroReach))';
                    unDecayed = signal(iniExpInd:zeroReach) - yFitted;

                    lateSigblFixed = signal(zeroReach:end) - (signal(zeroReach)-unDecayed(end));
                    cleanSignal(zeroReach:end) = lateSigblFixed;

                    cleanSignal(iniExpInd:zeroReach) = unDecayed;    
                end
        end
    end


    cleanSignal(tCut1:iniExpInd) = 0;

elseif isempty(iniExpInd)
    iniExpInd = NaN;

elseif isempty(endExpInd)
    endExpInd = NaN;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT the raw/fit model/clean signals: For checking inside the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

figure; plot(1000*time, signal,'b','linewidth',2);
hold on, plot(1000*time, cleanSignal,'k','linewidth',2);
plot(1000*time(iniExpInd:endExpInd), f0(time(iniExpInd:endExpInd)),'m','linewidth',2);
plot(1000*time(tCut2:endExpInd), signalSmoothed, 'c')
% plot(1000*time([iniExpInd, endExpPossibleInd]), signalSmoothed([iniExpInd-tCut2+1, end]), 'r')
legend('Raw signal','Cleaned signal','Fit model','Smoothed signal')
set(gca,'XLim',[-10 100])


figure;
plot(time([iniExpInd, endExpPossibleInd]),signalSmoothed([iniExpInd-tCut2+1, end]), 'r');
hold on
plot(time([x0 x3]), [y0 y3], 'g.-')
axis square

%}











