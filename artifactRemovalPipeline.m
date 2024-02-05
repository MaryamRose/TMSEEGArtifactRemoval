% =========================================================================
% The main function to run removing the artifacts from TMS-EEG data using
% the proposed method in: 
% "A Novel Approach to Artifact Removal in TMS-EEG Data Using Curve Fitting and Wavelet-Based Estimation",
% Maryam Rostami, Reza Zomorrodi, Reza Rostami, Gholam-Ali Hosseinzadeh
%
% Author of the codes: Maryam Rostami [mar.rostami@ut.ac.ir]
%
% The codes are shared publicly based on the requirements of the journal
% "Brain Topography" and are open-source with no waranty for the results.
% Please avoid using these codes for clinical and diagnostics applications.
% =========================================================================


clear all
clc


%% Add Paths and load files
% savefilepath = '~\preprocessed data\'; % the path to save the files
% addpath '~\eeglab2022.1'; % add EEGLAB to the path
eeglab; 


%% parameters

% general parameters
tms_artifact_CLIP=[-0.002 0.014]; % [cut1 cut2]
baselineCorrTime = [-0.800 -0.010]; % baseline correction time interval
srate = 5000; % enter sampling rate of the data


% decay removal parameters
tolTime = 2e-3; % tolerance time (in seconds)
noChangeTimeStart = 0.5; % from 500 msec. it should not change anymore
focusWeightms = 40e-3; % to put larger weights (more role) on this initial time after TMS pulse cut2


% blink detectin parameters
channelsToRemoveBlink = {'fp1','fp2','fpz','af7','af8','af3','af4'};
reducedSR = 250; % e.g. from 5000 to 250
dsRatio = srate/reducedSR; % down-sampling ratio
M = 200*1e-3*reducedSR; % Window length in embedding step 
clustNum = 3; % number of clusters
fdThreshold = 1.1; % the larger 'fdThreshold' value, more likely to remove slow waves (blinks) - a very small value could miss blinks, a very large value over-corrects
blinkAmpThr = 30; % a blink should be larger than 30 uV
blinkAreaRatio = 0.1; % Assigns how large the peak amplitude of a blink should be to remove/keep a blink in a channel (X/Fpz blink peak); the larger the parameter, the less removing
showClustersPic = 0; % Flag to check if the right cluster number was chosen


% eye movement detection parameters
eyeMoveAmpThr = -500; % multiplier of F7/F8 channels should be larger than this to be detected as eye movement
twinLen = round(20e-3*reducedSR); % windows of this length are compared between F7/F8 to detect the eye movement pulse
minMoveLen = 0.2*reducedSR; % the windows should be 200 ms long at least 
conservativeLen = 100e-3*srate; % add these samples to the strat and end of every piece of pulses 
moveAreaAmpRatio = 0.5; % the ratio to choose the right channels to remove the artifacts from
lagThreshold = 100e-3*dsRatio*reducedSR; % attach two pulses together if the spacing is less than 100 ms


% eye artifact correction parameters
correctionMethod = 'DWT-PCA';
waveletType = 'db3'; 
waveletDegree = 7;
eigThr = 0.5; % eigen-value threshold


eyeArtifactParams.M = M;
eyeArtifactParams.clustNum = clustNum;
eyeArtifactParams.fdThreshold = fdThreshold;
eyeArtifactParams.blinkAmpThr = blinkAmpThr;
eyeArtifactParams.blinkAreaRatio = blinkAreaRatio;
eyeArtifactParams.eyeMoveAmpThr = eyeMoveAmpThr;
eyeArtifactParams.twinLen = twinLen;
eyeArtifactParams.minMoveLen = minMoveLen;
eyeArtifactParams.conservativeLen = conservativeLen;
eyeArtifactParams.moveAreaAmpRatio = moveAreaAmpRatio;
eyeArtifactParams.eigThr = eigThr;
eyeArtifactParams.waveletDegree = waveletDegree;



%% loading the data and channel file

[fileName, filePath] = uigetfile('.set','Select the raw data file in .set format.');
EEG = pop_loadset(fullfile(filePath, fileName));
[~, subjectName, ~] = fileparts(fileName);
fprintf('\n\n*** Dataset name: %s **\n\n',subjectName);  

[chanfilename, chanfilenamePath] = uigetfile({'.locs';'.ced';'.loc'},'Select the channels file.');
chanfile = fullfile(chanfilenamePath, chanfilename);
[allChanLocs, ~] = pop_readlocs(chanfile); 
allChanLabels = {allChanLocs.labels};


%% Primary Processing steps
          
EEG = pop_select( EEG,'nochannel',{'ECG','EOG','EMG'});  % omit unneccessary channels
EEG = pop_chanedit(EEG, 'lookup', chanfile);
EEG.allchan = EEG.chanlocs;
EEG = pop_epoch( EEG, {  'condition 2'  }, [-1  2], 'epochinfo', 'yes'); % tag for stimulation
EEG = pop_rmbase(EEG, 1e3*baselineCorrTime);
raw_epochs_eeg = double(EEG.data);


%% Remove the decay artifact 

% prelocating
incomleteFitInfo = []; incompleteFitInd = 0;    
cleanSignalUD = double(EEG.data); % in case the fit gets skipped, leave the signal the same as before 


% complementary variables for decay removal 
timeSec = EEG.times./1000;
tCut1 = dsearchn(timeSec', tms_artifact_CLIP(1));
tCut2 = dsearchn(timeSec', tms_artifact_CLIP(2));
zeroSample = dsearchn(timeSec',0);             
ind500 = dsearchn(timeSec', noChangeTimeStart); 
iniExpInd = dsearchn(timeSec', tms_artifact_CLIP(2)); % initial point of the exponential decay
lenTemp = ind500-1-iniExpInd+1; % from the start of the exponential to 500 ms 
indEnd = length(timeSec); % total length of an epoch

% perform a fit on average first, to get an estimate of the fit parameters and use them as initial points for trial-wise fit
rng('default')
signalAvg = double(mean(EEG.data,3));
for ich = 1:size(signalAvg,1)

    tt = timeSec(iniExpInd:indEnd)';
    xx = signalAvg(ich,iniExpInd:indEnd)';
    ww = [1e10 1e4*ones(1,focusWeightms*srate-1) ones(1,lenTemp-focusWeightms*srate-1) 1e4*ones(1,indEnd-ind500+2)]';

    f0 = fit(tt, xx, fittype('b*exp(c*x)'), ...
        'TolFun',1e-4, 'TolX',1e-4,'Weights',ww,...
        'lower', [], 'upper', [], 'StartPoint', randn(1,2)); 

    startPoint(ich).b = f0.b;
    startPoint(ich).c = f0.c;

    if ~isempty(f0)
        avgEstimate = [f0.b f0.c];
        lowerLim{ich} = -10*abs(avgEstimate);
        upperLim{ich} = 10*abs(avgEstimate);
    else
        lowerLim{ich} = [];
        upperLim{ich} = [];
    end
end


% Main block of channel and trial-wise decay removal
iniExpInd = zeros(size(cleanSignalUD));
endExpInd = iniExpInd;
fprintf('\nUnDecaying is being applied. \nTrial ::: ')
for trialNum = 1:EEG.trials
    fprintf('%d,', trialNum);
    for chanNum = 1:EEG.nbchan
        try
            signal = double(EEG.data(chanNum,:,trialNum));
            [cleanSignalUD(chanNum,:,trialNum), iniExpInd(chanNum,:,trialNum), endExpInd(chanNum,:,trialNum), ~] = ...
                unDecaying_v9(signal, tms_artifact_CLIP, baselineCorrTime, EEG.srate, EEG.times./1000, ...
                zeroSample, [], tolTime, noChangeTimeStart, [-1e4 lowerLim{chanNum}], [1e4 upperLim{chanNum}], startPoint(chanNum), focusWeightms);   
        catch ME
            incompleteFitInd = incompleteFitInd + 1;
            incomleteFitInfo(incompleteFitInd).message = ME.message;
            incomleteFitInfo(incompleteFitInd).subjSesTrCh = [trialNum chanNum]; % trial and channel that has caused error in fitting 
        end
     end
end 
fprintf('\n\n');

cleanSignalUD(:,tCut1:tCut2,:) = repmat(cleanSignalUD(:,tCut2,:), [1, tCut2-tCut1+1]);           
EEG.data = single(cleanSignalUD);

% save UnDecayed data 
if ~isfolder(savefilepath)
    mkdir(savefilepath);
end
pop_saveset(EEG, fullfile(savefilepath, [subjectName '_UD.set']));


%% detect the blinks in channel data

EEG = pop_rmbase(EEG, 1e3*baselineCorrTime);                                    
EEG = pop_tesa_removedata( EEG, [1e3*tms_artifact_CLIP(1), 1e3*tms_artifact_CLIP(2)] );
EEG = pop_tesa_interpdata( EEG, 'linear', [1,1]);   

EEGds = pop_resample( EEG, reducedSR); 

blinkChanInds = find(ismember(lower({EEG.allchan.labels}), lower(channelsToRemoveBlink)));

blinkPulseMat = zeros(EEG.trials,EEG.pnts,EEG.nbchan);
for trialNum = 1:EEG.trials
    fprintf('%d,', trialNum);
    xTr = double(EEG.data(blinkChanInds,:,trialNum));
    xTrds = double(EEGds.data(blinkChanInds,:,trialNum));
    blinkPulse = unBlinkWithFDCollective(xTrds, M, clustNum, fdThreshold, blinkAmpThr, dsRatio, lagThreshold, reducedSR, conservativeLen, showClustersPic);
    blinkPulseDS = downsample(blinkPulse, EEG.srate/reducedSR);  

    bpdiff = diff([0 blinkPulseDS 0]);
    tStarts = find(bpdiff==1);
    tEnds = find(bpdiff==-1)-1;

    signaldsfpz = double(squeeze(EEGds.data(find(ismember(lower({EEG.allchan.labels}), 'fpz')),:,trialNum))); % select a channel close to FPz, if it is not available in your data (e.g. Fz, AFz)
    for chanNum = 1:EEG.nbchan
        signalds = double(squeeze(EEGds.data(chanNum,:,trialNum)));

        if ~isempty(tStarts)
            [sig1, sig4, areaSigCurChan, areaSigFpz] = deal([]);
            for iSeg = 1:length(tStarts)                                
                sig1Temp = signalds(tStarts(iSeg):tEnds(iSeg)); sig1Temp = sig1Temp - sig1Temp(1); 
                sig1 = [sig1, sig1Temp];
                sig4Temp = signaldsfpz(tStarts(iSeg):tEnds(iSeg)); sig4Temp = sig4Temp - sig4Temp(1); 
                sig4 = [sig4, sig4Temp];
                areaSigCurChan = [areaSigCurChan trapz(1:length(sig1), abs(sig1))];
                areaSigFpz = [areaSigFpz trapz(1:length(sig4), abs(sig4))];                         
            end


            if sum(areaSigCurChan./areaSigFpz > blinkAreaRatio) == length(areaSigFpz)
                blinkPulseMat(trialNum, :, chanNum) = blinkPulse;
            end

        end                    

    end

end
fprintf('\nDone!\n')


%% Detect the time interval for eye movement 

movePulseMat = zeros(EEG.trials, EEG.pnts, EEG.nbchan);
fprintf('Trial number:');
for trialNum=1:EEG.trials
    fprintf('%d,',trialNum);

    f7f8DataTemp = double(EEGds.data(find(ismember(lower({EEG.allchan.labels}), {'f7','f8'})),:,trialNum));  % F7 and F8 (fronto-lateral) channels              
    f7f8Data = [mean(f7f8DataTemp(1:size(f7f8DataTemp,1)/2,:),1); mean(f7f8DataTemp(size(f7f8DataTemp,1)/2+1:end,:),1)];

    movePulse = unEyeMovementWithFDCollective(f7f8Data,eyeMoveAmpThr,twinLen,minMoveLen,conservativeLen,lagThreshold,dsRatio);
    movePulseDS = downsample(movePulse, dsRatio);     

    mvdiff = diff([0 movePulseDS 0]);
    tStarts = find(mvdiff==1);
    tEnds = find(mvdiff==-1)-1;

    signaldsf7 = double(squeeze(EEGds.data(find(ismember(lower({EEG.allchan.labels}), {'f7'})),:,trialNum)));
    signaldsf8 = double(squeeze(EEGds.data(find(ismember(lower({EEG.allchan.labels}), {'f8'})),:,trialNum)));
    for chanNum = 1:EEG.nbchan
        signalds = double(squeeze(EEGds.data(chanNum,:,trialNum)));

        if ~isempty(tStarts)
            [sig1, sig2, sig3, areaSigCurChan, areaSigFp1, areaSigFp2] = deal([]);
            for iSeg = 1:length(tStarts)                
                sig1Temp = signalds(tStarts(iSeg):tEnds(iSeg)); 
                sig1 = [sig1, sig1Temp];
                sig2Temp = signaldsf7(tStarts(iSeg):tEnds(iSeg)); 
                sig2 = [sig2, sig2Temp];
                sig3Temp = signaldsf8(tStarts(iSeg):tEnds(iSeg));  
                sig3 = [sig3, sig3Temp];     
                areaSigCurChan = [areaSigCurChan trapz(1:length(sig1), abs(sig1))];
                areaSigFp1 = [areaSigFp1 trapz(1:length(sig2), abs(sig2))];
                areaSigFp2 = [areaSigFp2 trapz(1:length(sig3), abs(sig3))];
            end

             if sum([areaSigCurChan areaSigCurChan]./[areaSigFp1 areaSigFp2] > moveAreaAmpRatio) == length([areaSigFp1 areaSigFp2])
                movePulseMat(trialNum, :, chanNum) = movePulse;
            end

        end

    end

end
fprintf('\nDone!\n')



%% Check fitting of the artifact pulses
   
eyePulse = blinkPulseMat | movePulseMat;

figure('unit','normalized','position',[0 0 1 1])
checkPulseFit(EEG, eyePulse, [0*3 10*3], 'epochs', [])
% set(gca,'XLim',[60*3 70*3]); % to update the view of the figure to see
% e.g. 60th to 70th trial (3 sec. is the length of the epoch)


%% Now correct using the pulse for both eye movement and blinks

EEG_NoEyeArtifact = EEG; % prelocate
cleanedSignal = DWTPCAcorrection(EEG, eyePulse, waveletDegree, waveletType, eigThr, 'epoch', []);
fprintf('\nDone!\n')

EEG_NoEyeArtifact.data = cleanedSignal;
EEG_NoEyeArtifact.eyeArtifactParams = eyeArtifactParams;

% To visually check and compare raw and eye-artifact cleaned
% eeglab redraw
% [ALLEEG, ~, ~] = pop_newset(ALLEEG,EEG_NoEyeArtifact,1);
% eeglab redraw

% Save the eye artifact cleaned data
if ~isfolder(savefilepath)
    mkdir(savefilepath);
end
pop_saveset(EEG_NoEyeArtifact, fullfile(savefilepath, [subjectName '_UD_noEyeArtifact.set']));

EEG = EEG_NoEyeArtifact;

%% Apply post-processing steps to remove the remainder of artifacts 

% mark bad trials and channels and remove them
eeglab redraw 
pop_eegplot(EEG);
input('\nHighlight bad trials and press REJECT, then press enter. ');
badChannels = input('\nWrite the name of bad channels: ');
EEG.preprocessing.badTrials = find(EEG.reject.rejmanual);
EEG.preprocessing.badChannels = badChannels;         
EEG = pop_rejepoch( EEG, EEG.preprocessing.badTrials ,0); 
EEG = pop_select( EEG,'nochannel',EEG.preprocessing.badChannels);

% Apply SSP-SIR to reject muscle artifacts
EEG = pop_tesa_sspsir( EEG, 'artScale','automatic','PC',[]); % artifact composes 10% variance of the data
fprintf('PCs removed by SSP-SIR: %d\n', EEG.sspsirPCsRemoved);

% Apply bandpass and notch filter 
EEG = pop_tesa_filtbutter( EEG, 0.5, 90, 4, 'bandpass' );
EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' );

% Re-reference to average, interpolate missing channels and remove baseline
EEG = pop_reref( EEG, []);
EEG = pop_interp(EEG, EEG.allchan, 'spherical');
EEG = pop_rmbase(EEG, baselineCorrTime);

%% check the final TEP results

eeglab redraw 
ChanIndToMark = find(ismember(allChanLabels,'F3'));
figure; plot(EEG.times, mean(EEG.data,3),'b');
hold on; plot(EEG.times, mean(EEG.data(ChanIndToMark,:,:),3), 'r', 'linewidth',2);
axis([-100 500 -20 20]); grid minor; grid on
            
figure('units','normalized','position',[0 0 1 1]);
tep_topoplot_compare_data(mean(raw_epochs_eeg,3),mean(EEG.data,3),EEG.times,chanfile,ChanIndToMark,tms_artifact_CLIP,[-100 500],[-50 50],'Raw','Clean')

