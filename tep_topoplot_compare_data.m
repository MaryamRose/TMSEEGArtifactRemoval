% ===================================================================================================================
% The function to plot the TEPs and topographies of two datasets (e.g. raw
% data vs. cleaned data) for comparison purposes
% 
% tep_topoplot_compare_data(data1,data2,times,chanFileName,ChanIndToMark,tms_artifact_CLIP,Xlim,Ylim,title1,title2)
%       data1: the first 2-dimensional data to plot (size: channels x time_points)
%       data2: the second 2-dimensional data to plot (size: channels x time_points)
%       times: time array in ms.
%       chanFileName: the path to channels file (in .ced or .locs format)
%       ChanIndToMark: The channel index to mark (e.g. F3)
%       tms_artifact_CLIP: the TMS pulse time range that was interpolated
%                           in seconds
%       Xlim: x-axis time range in msec.
%       Ylim: y-axis range in micro-volts
%       title1: the title for the left plot (showing data1)
%       title2: the title for the right plot (showing data2)
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

function tep_topoplot_compare_data(data1,data2,times,chanFileName,ChanIndToMark,tms_artifact_CLIP,Xlim,Ylim,title1,title2)

    P30_times=find(times >=25 & times <=35); %Time interval for P30
    N45_times=find(times >=40 & times<=55); %Time interval for N45 (why -40???)
    P60_times=find(times >55 & times <=65); %Time interval for P60
    N100_times=find(times > 90 & times <=120); %Time interval for N100

    T{1}=P30_times; T{2}=N45_times; T{3}=P60_times; T{4}=N100_times;
    tl={'P30','N45','P60','N100'};
        
    
    if ~isempty(data2)
        subplot(2,8,1:4);
        plot(times, squeeze(mean(data1,3)),'b','LineWidth',1 ),hold on;
        plot(times, squeeze(mean(data1(ChanIndToMark,:,:),3)),'r','LineWidth',2 )
        xlim(Xlim),ylim(Ylim), grid on
        xlabel('Time(ms)'), ylabel('TEP (\muV)'), title(title1)
        patch(1e3*[tms_artifact_CLIP(1) tms_artifact_CLIP(2) tms_artifact_CLIP(2) tms_artifact_CLIP(1)], 10*[Ylim(1) Ylim(1) Ylim(2) Ylim(2)],0.7*[1 1 1],...
            'markerEdgeColor','none','LineStyle','none')


        subplot(2,8,5:8);
        plot(times, squeeze(mean(data2,3)),'b','LineWidth',1 ),hold on;
        plot(times, squeeze(mean(data2(ChanIndToMark,:,:),3)),'r','LineWidth',2 )
        xlim(Xlim),ylim(Ylim), grid on
        xlabel('Time(ms)'), title(title2)
        patch(1e3*[tms_artifact_CLIP(1) tms_artifact_CLIP(2) tms_artifact_CLIP(2) tms_artifact_CLIP(1)], 10*[Ylim(1) Ylim(1) Ylim(2) Ylim(2)],0.7*[1 1 1],...
            'markerEdgeColor','none','LineStyle','none')



        for i=1:4
            subplot(2,8,8+i);
            C= squeeze(  mean( mean(data1(:,T{i},:),3),2)   );
            topoplot(C,chanFileName,'plotrad',.55,'headrad',.55,'shading','interp','whitebk','on', 'electrodes','off','style','map');%,colorbar  %'maplimits',[0 1.5],   \
            title(tl(i));    
            cbar('horiz')
        end


        for i=1:4
            subplot(2,8,12+i);
            C= squeeze(  mean( mean(data2(:,T{i},:),3),2)   );
            topoplot(C,chanFileName,'plotrad',.55,'headrad',.55,'shading','interp','whitebk','on', 'electrodes','off','style','map');%,colorbar  %'maplimits',[0 1.5],   \
            title(tl(i));   
            cbar('horiz')
        end
        
    else
        
        subplot(2,4,1:4);
        plot(times, squeeze(mean(data1,3)),'b','LineWidth',1 ),hold on;
        plot(times, squeeze(mean(data1(ChanIndToMark,:,:),3)),'r','LineWidth',2 )
        xlim(Xlim),ylim(Ylim), grid on
        xlabel('Time(ms)'), ylabel('TEP (\muV)'), title(title1)
        patch(1e3*[tms_artifact_CLIP(1) tms_artifact_CLIP(2) tms_artifact_CLIP(2) tms_artifact_CLIP(1)], 10*[Ylim(1) Ylim(1) Ylim(2) Ylim(2)],0.7*[1 1 1],...
            'markerEdgeColor','none')
        
        for i=1:4
            subplot(2,4,4+i);
            C= squeeze(  mean( mean(data1(:,T{i},:),3),2)   );
            topoplot(C,chanFileName,'plotrad',.55,'headrad',.55,'shading','interp','whitebk','on', 'electrodes','off','style','map');%,colorbar  %'maplimits',[0 1.5],   \
            title(tl(i));    
        end

        
    end


