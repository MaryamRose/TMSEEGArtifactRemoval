% ==========================================================================================================
% The function to plot EEG channels signals over time and the detected
% eye-artifact pulses on top of it for visual check of the artifact
% detection
% 
% checkPulseFit(EEG, eyePulse, xLim, dataType, secLength)
%
%       EEG: input data in EEGLAB format
%       eyePulse: The binary (0/1) array representing the detected time of the
%                 eye artifacts
%       xLim: The x-axis range to show in seconds(e.g. [120 130])
%       dataType: 'epochs' for epoched data, or 'rest' for resting state
%                 continuous data
%       secLength: in case of 'rest' dataType, the length of the window to show
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
% ==========================================================================================================

function checkPulseFit(EEG, eyePulse, xLim, dataType, secLength)

    switch dataType
        case 'epochs'
            yDiv = 100 ;
            xDiv = EEG.pnts / EEG.srate;
            for i=1:EEG.nbchan
                sigToPlot = reshape(squeeze(EEG.data(i,:,:)),1,EEG.pnts*EEG.trials); 
                pulseToPlot = reshape(squeeze(eyePulse(:,:,i))',1,EEG.pnts*EEG.trials);
                timeToPlot = 0:1/EEG.srate:EEG.trials*(EEG.pnts/EEG.srate)-1/EEG.srate;
                dcLine = (EEG.nbchan-i+1) * yDiv;
                plot(timeToPlot, sigToPlot+dcLine, 'b', 'linewidth', 1);
                hold 'on'; 
                plot(timeToPlot, 50*pulseToPlot+dcLine, 'r');
                set(gca,'TickLength', [0 0],'FontSize',7.5,...
                    'YTick',[],'XLim',[0 5]);
            end
            set(gca,'YLim',[-100 6500]);
            set(gca,'XLim',xLim)
            xlabel('Time (sec.)')
            ylabel('Voltage (\muV)')


            for iL=1:EEG.trials-1
                plot([xDiv*iL xDiv*iL], [-100 6500], 'color', 0.5*[1 1 1])      
            end
            
        case 'rest'
            yDiv = 100 ;
            xDiv = secLength;
            for i=1:EEG.nbchan
                sigToPlot = detrend(EEG.data(i,:),1); 
                pulseToPlot = eyePulse(i,:);
                timeToPlot = 0:1/EEG.srate:EEG.trials*(EEG.pnts/EEG.srate);
                dcLine = (EEG.nbchan-i+1) * yDiv;
                plot(timeToPlot(1:length(pulseToPlot)), sigToPlot(1:length(pulseToPlot))+dcLine, 'b', 'linewidth', 1);
                hold 'on'; 
                plot(timeToPlot(1:length(pulseToPlot)), 50*pulseToPlot+dcLine, 'r');
                set(gca,'TickLength', [0 0],'FontSize',7.5,...
                    'YTick',[],'XLim',[0 5]);
            end
            set(gca,'YLim',[-100 6500]);
            set(gca,'XLim',xLim)
            xlabel('Time (sec.)')
            ylabel('Voltage (\muV)')


            for iL=1:size(eyePulse,2)/EEG.srate/secLength-1
                plot([xDiv*iL xDiv*iL], [-100 6500], 'color', 0.5*[1 1 1])      
            end
    end
            
   
end