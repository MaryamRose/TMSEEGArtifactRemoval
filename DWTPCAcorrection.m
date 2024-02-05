% ==========================================================================================================
% The function to apply the correction of the eye artifacts that were
% detected in previous steps (using "unBlinkWithFDCollective" and
% "unEyeMovementWithFDCollective")
% 
% [cleanedSignal] = DWTPCAcorrection(EEG, eyePulse, waveletDegree, waveletType, eigThr, dataType, sectionLen)
%
%       cleanedSignal: the cleaned signal from detected artifacts (size: the same size of EEG.data)
%       EEG: input data in EEGLAB format
%       eyePulse: The binary (0/1) array representing the detected time of the
%                 eye artifacts
%       waveletDegree: The number of wavelet levels to decompose the signal
%       waveletType: the type of the wavelet (e.g. db3')
%       eigThr: The threshold to construct the artifact space using the eigenvectors of the signal within
%               the identified artifact pulse
%       dataType: 'epochs' for epoched data, or 'rest' for resting state
%                 continuous data
%       sectionLen: in case of 'rest' dataType, the length of the window to show
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

function [cleanedSignal] = DWTPCAcorrection(EEG, eyePulse, waveletDegree, waveletType, eigThr, dataType, sectionLen)
     
    switch dataType
        case 'epoch'
            cleanedSignal = EEG.data;
            fprintf('Trial number:');
            for trNum = 1:size(eyePulse,1)
                fprintf('%d,',trNum);   
                approxMat = [];
                [c,l] = deal(cell(1,EEG.nbchan));

                for iiCh = 1:EEG.nbchan
                        signalFull = double(squeeze(EEG.data(iiCh,:,trNum)));
                        signal = signalFull .* eyePulse(trNum,:,iiCh);

                        [c{iiCh},l{iiCh}] = wavedec(signal,waveletDegree,waveletType);
                        approxMat(iiCh,:) = appcoef(c{iiCh},l{iiCh},waveletType);
                end
                
                
                nonZeroChans = sum(approxMat,2)~=0;
                pulseDiff = diff([zeros(size(approxMat,1),1) approxMat zeros(size(approxMat,1),1)],[],2);
                tStartInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1), 1:size(pulseDiff,1), 'un',0);
                tStartInds = [tStartInds{:}];
                tEndInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1,'last')-1, 1:size(pulseDiff,1), 'un',0);
                tEndInds = [tEndInds{:}];
                
                tStartsIndsOne = min(tStartInds);
                tEndIndsOne = max(tEndInds);
                
                approxMatNonZero = zeros(size(approxMat));
                if ~isempty(tStartInds)
                    approxMatNonZero = zeros(size(approxMat,1), tEndIndsOne-tStartsIndsOne+1);
                    approxMatNonZero(nonZeroChans,:) = approxMat(nonZeroChans, tStartsIndsOne:tEndIndsOne);
                end
                
                % apply SVD on channels that have some artifact pattern unless all channels show zero artifact
                if sum(approxMatNonZero)~=0
                    [~, eigVal, eigVec] = svd(approxMatNonZero(nonZeroChans,:));
                else 
                    [~, eigVal, eigVec] = svd(approxMatNonZero);
                end
                
                D = diag(eigVal);

                eigValSum = cumsum(D)./sum(D);
                eigNumInd = find(eigValSum>=eigThr,1);

                for iiCh1  = 1:EEG.nbchan

                    signalFull = double(squeeze(EEG.data(iiCh1,:,trNum)));
                    b1 = glmfit(eigVec(:,1:eigNumInd), approxMatNonZero(iiCh1,:));

                    artifactModel = zeros(1,size(approxMat,2));
                    if nonZeroChans(iiCh1)
                        artifactModel(tStartsIndsOne:tEndIndsOne) = ...
                            b1'*[ones(size(eigVec,1),1) eigVec(:,1:eigNumInd)]'; % main components of the approximation
                    end
                    
                    cFixed = [artifactModel zeros(1,length(c{iiCh1})-l{iiCh1}(1))];
                    artifactRecTemp = waverec(cFixed,l{iiCh1},waveletType);
                    artifactRec = zeros(1,length(signalFull));
                    artifactRec( eyePulse(trNum,:,iiCh1) ) = artifactRecTemp( eyePulse(trNum,:,iiCh1) );
                    signalRec = signalFull - artifactRec;

        %             plotDWTPCAResults(signalFull, signalRec,artifactRec, eyePulse(trNum,:,iiCh1), EEG.srate)

                    pulseDiff = diff([0 eyePulse(trNum,:,iiCh1) 0]);
                    tStarts = find(pulseDiff==1);
                    tEnds = find(pulseDiff==-1)-1;

                    unqVals = setdiff(unique(diff(eyePulse(trNum,:,iiCh1))),0);
                    tStarts0 = []; 
                    tEnds0 = [];
                    for iii = 1:length(unqVals)
                        if unqVals(iii)==1
                            tEnds0 = find(diff(eyePulse(trNum,:,iiCh1))==unqVals(iii));
                        elseif unqVals(iii)==-1
                            tStarts0 = find(diff(eyePulse(trNum,:,iiCh1))==unqVals(iii))+1;
                        end
                    end


                    if sum(eyePulse(trNum,:,iiCh1))~=EEG.pnts && ~isempty(tStarts) && ~isempty(tEnds)
                        if ~isempty(tStarts0) && isempty(tEnds0)
                            tEnds0 = length(eyePulse(trNum,:,iiCh1));
                        end
                        if isempty(tStarts0) && ~isempty(tEnds0)
                            tStarts0 = 1;
                        end
                        if isempty(tStarts0(tStarts0 < tEnds0(1)))
                            tStarts0 = [1 tStarts0];
                        end
                        if isempty(tEnds0(tEnds0 > tStarts0(end)))
                            tEnds0 = [tEnds0 length(eyePulse(trNum,:,iiCh1))];
                        end
                    end


                    lineSig0 = zeros(1,length(signalFull));
                    for tInd = 1:length(tStarts)
                        t1 = tStarts(tInd); 
                        t2 = tEnds(tInd);
                        lineAround0 = signalRec(t1)+(signalRec(t2)-signalRec(t1))/(t2-t1)*((t1:t2)-t1); 
                        lineSig0(t1:t2) = lineAround0;
                    end


                    tempSig = signalRec - lineSig0; 
                    cleanedSignal(iiCh1,:,trNum) = tempSig;

                    if ~isempty(tStarts) && ~isempty(tEnds) && ~isempty(tStarts0) && ~isempty(tEnds0)
                        breakPoints = [tStarts tStarts0; tEnds tEnds0];

                        for bpi = 1:size(breakPoints,2)
                            interval = breakPoints(1,bpi):breakPoints(2,bpi);
                            cleanedSignal(iiCh1,interval,trNum) = detrend(tempSig(interval), 1);
                        end
                    end 
                end


            end
            
        case 'rest'
            
            sectionsNumber = ceil(EEG.xmax/sectionLen);
            cleanedSignal = EEG.data;
            fprintf('Section number:');
            for iSec = 1:sectionsNumber
                fprintf('%d,',iSec);   
                approxMat = [];
                [c,l] = deal(cell(1,EEG.nbchan));
                if iSec~=sectionsNumber
                    lenSecIndsCurrent = 1+(iSec-1)*sectionLen*EEG.srate : iSec*sectionLen*EEG.srate;  
                else
                    lenSecIndsCurrent = 1+(iSec-1)*sectionLen*EEG.srate : EEG.pnts;                      
                end
                                

                for iiCh = 1:EEG.nbchan
                    lenSecInds = lenSecIndsCurrent;                        
                    signalFull = double(squeeze(EEG.data(iiCh,lenSecInds)));
                    signalFull = detrend(signalFull,0); % demean
                    signal = signalFull .* eyePulse(iiCh,lenSecInds);

                    [c{iiCh},l{iiCh}] = wavedec(signal,waveletDegree,waveletType);
                    approxMat(iiCh,1:l{iiCh}) = appcoef(c{iiCh},l{iiCh},waveletType);
                end
                
                nonZeroChans = sum(approxMat,2)~=0;
                pulseDiff = diff([zeros(size(approxMat,1),1) approxMat zeros(size(approxMat,1),1)],[],2);
                tStartInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1), 1:size(pulseDiff,1), 'un',0);
                tStartInds = [tStartInds{:}];
                tEndInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1,'last')-1, 1:size(pulseDiff,1), 'un',0);
                tEndInds = [tEndInds{:}];
                
                tStartsIndsOne = min(tStartInds);
                tEndIndsOne = max(tEndInds);
                
                approxMatNonZero = zeros(size(approxMat));
                if ~isempty(tStartInds)
                    approxMatNonZero = zeros(size(approxMat,1), tEndIndsOne-tStartsIndsOne+1);
                    approxMatNonZero(nonZeroChans,:) = approxMat(nonZeroChans, tStartsIndsOne:tEndIndsOne);
                end
                
               
                % apply SVD on channels that have some artifact pattern unless all channels show zero artifact
                if sum(approxMatNonZero)~=0
                    [~, eigVal, eigVec] = svd(approxMatNonZero(nonZeroChans,:));
                else 
                    [~, eigVal, eigVec] = svd(approxMatNonZero);
                end
                
                D = diag(eigVal);
                eigValSum = cumsum(D)./sum(D);
                eigNumInd = find(eigValSum>=eigThr,1);
                

                for iiCh1  = 1:EEG.nbchan
                    

                    lenSecInds = lenSecIndsCurrent;                        

                    signalFull = double(squeeze(EEG.data(iiCh1,lenSecInds)));
                    b1 = glmfit(eigVec(:,1:eigNumInd), approxMatNonZero(iiCh1,:));

                    artifactModel = zeros(1,size(approxMat,2));
                    if nonZeroChans(iiCh1)
                        artifactModel(tStartsIndsOne:tEndIndsOne) = ...
                            b1'*[ones(size(eigVec,1),1) eigVec(:,1:eigNumInd)]'; % main components of the approximation
                    end
                    
                    cFixed = [artifactModel zeros(1,length(c{iiCh1})-l{iiCh1}(1))];
                    artifactRecTemp = waverec(cFixed,l{iiCh1},waveletType);
                    artifactRec = zeros(1,length(signalFull));
                    artifactRec( eyePulse(iiCh1,lenSecInds) ) = artifactRecTemp( eyePulse(iiCh1,lenSecInds) );
                    signalRec = signalFull - artifactRec;

%                     plotDWTPCAResults(signalFull, signalRec,artifactRec, eyePulse(iiCh1,lenSecInds), EEG.srate)

                    pulseDiff = diff([0 eyePulse(iiCh1,lenSecInds) 0]);
                    tStarts = find(pulseDiff==1);
                    tEnds = find(pulseDiff==-1)-1;

                    unqVals = setdiff(unique(diff(eyePulse(iiCh1,lenSecInds))),0);
                    tStarts0 = []; 
                    tEnds0 = [];
                    for iii = 1:length(unqVals)
                        if unqVals(iii)==1
                            tEnds0 = find(diff(eyePulse(iiCh1,lenSecInds))==unqVals(iii));
                        elseif unqVals(iii)==-1
                            tStarts0 = find(diff(eyePulse(iiCh1,lenSecInds))==unqVals(iii))+1;
                        end
                    end


                    if sum(eyePulse(iiCh1,lenSecInds))~=length(signalFull) && ~isempty(tStarts) && ~isempty(tEnds)
                        if ~isempty(tStarts0) && isempty(tEnds0)
                            tEnds0 = length(eyePulse(iiCh1,lenSecInds));
                        end
                        if isempty(tStarts0) && ~isempty(tEnds0)
                            tStarts0 = 1;
                        end
                        if isempty(tStarts0(tStarts0 < tEnds0(1)))
                            tStarts0 = [1 tStarts0];
                        end
                        if isempty(tEnds0(tEnds0 > tStarts0(end)))
                            tEnds0 = [tEnds0 length(eyePulse(iiCh1,lenSecInds))];
                        end
                    end
                    
                    
                    lineSig = zeros(1,length(signalRec));
                    lineSig0 = zeros(1,length(signalRec));
                    for tInd = 1:length(tStarts)
                        t1 = tStarts(tInd); 
                        t2 = tEnds(tInd);
                        lineSigPiece = signalFull(t1)+(signalFull(t2)-signalFull(t1))/(t2-t1)*((t1:t2)-t1);   
                        lineSig(t1:t2) = lineSigPiece;
                    end

                    tempSig = signalRec  - lineSig0 + lineSig; 
                    cleanedSignal(iiCh1,lenSecInds) = tempSig;


                end


            end
            
            
            %%%->>>>> Complementary processing <<<<<<-
            %%% Check for eyePulse picese on the intersection of two segments
            segmentsEdges = zeros(1,size(eyePulse,2));
            segmentsEdges(sectionLen*EEG.srate:sectionLen*EEG.srate:end) = 1;
            segmentsEdges(sectionLen*EEG.srate-1:sectionLen*EEG.srate:end) = 1;
            segmentsEdges(sectionLen*EEG.srate+1:sectionLen*EEG.srate:end) = 1;
            pulseOnEdge = eyePulse & segmentsEdges;
            [~,cp] = find(pulseOnEdge==1);
            edgeTimes = unique(cp);
            edgeTimesUnique = edgeTimes(ismember(edgeTimes,edgeTimes-1) & ismember(edgeTimes,edgeTimes+1));
            fprintf('Checking the artifacts on the edge of sections: ');
            for ij=1:length(edgeTimesUnique)
                    fprintf('%d,',ij);   
                    approxMat = [];
                    [c,l] = deal(cell(1,EEG.nbchan));
                    for iiCh = 1:EEG.nbchan

                        eyePulse1 = eyePulse(iiCh,1+edgeTimesUnique(ij)-sectionLen*EEG.srate:edgeTimesUnique(ij)); 
                        eyePulse2 = eyePulse(iiCh,1+edgeTimesUnique(ij):sectionLen*EEG.srate+edgeTimesUnique(ij)); 

                        pulseDiff = diff([0 eyePulse1 0]);
                        tStarts = find(pulseDiff==1,1,'last')+edgeTimesUnique(ij)-sectionLen*EEG.srate;
                        pulseDiff = diff([0 eyePulse2 0]);
                        tEnds = find(pulseDiff==-1,1,'first')-1+edgeTimesUnique(ij);
                        lenSecInds = tStarts:tEnds;
                        
                        if ~isempty(lenSecInds)
                            signalFull = double(squeeze(EEG.data(iiCh,lenSecInds)));
                            signalFull = detrend(signalFull,0); 
                            signal = signalFull .* eyePulse(iiCh,lenSecInds);

                            [c{iiCh},l{iiCh}] = wavedec(signal,waveletDegree,waveletType);
                            approxMat(iiCh,1:l{iiCh}) = appcoef(c{iiCh},l{iiCh},waveletType);
                        else
                            [c{iiCh},l{iiCh}] = deal([]);
                            approxMat(iiCh,1:l{iiCh}) = zeros(1,length(lenSecInds));
                        end
                    end

                    nonZeroChans = sum(approxMat,2)~=0;
                    pulseDiff = diff([zeros(size(approxMat,1),1) approxMat zeros(size(approxMat,1),1)],[],2);
                    tStartInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1), 1:size(pulseDiff,1), 'un',0);
                    tStartInds = [tStartInds{:}];
                    tEndInds = arrayfun(@(x) find(pulseDiff(x,:)~=0,1,'last')-1, 1:size(pulseDiff,1), 'un',0);
                    tEndInds = [tEndInds{:}];

                    tStartsIndsOne = min(tStartInds);
                    tEndIndsOne = max(tEndInds);

                    approxMatNonZero = zeros(size(approxMat));
                    if ~isempty(tStartInds)
                        approxMatNonZero = zeros(size(approxMat,1), tEndIndsOne-tStartsIndsOne+1);
                        approxMatNonZero(nonZeroChans,:) = approxMat(nonZeroChans, tStartsIndsOne:tEndIndsOne);
                    end

                    if sum(approxMatNonZero)~=0
                        [~, eigVal, eigVec] = svd(approxMatNonZero(nonZeroChans,:));
                    else 
                        [~, eigVal, eigVec] = svd(approxMatNonZero);
                    end

                    D = diag(eigVal);
                    eigValSum = cumsum(D)./sum(D);
                    eigNumInd = find(eigValSum>=eigThr,1);

                    for iiCh1  = 1:EEG.nbchan

                        eyePulse1 = eyePulse(iiCh1,1+edgeTimesUnique(ij)-sectionLen*EEG.srate:edgeTimesUnique(ij)); 
                        eyePulse2 = eyePulse(iiCh1,1+edgeTimesUnique(ij):sectionLen*EEG.srate+edgeTimesUnique(ij)); 

                        pulseDiff = diff([0 eyePulse1 0]);
                        tStarts = find(pulseDiff==1,1,'last')+edgeTimesUnique(ij)-sectionLen*EEG.srate;
                        pulseDiff = diff([0 eyePulse2 0]);
                        tEnds = find(pulseDiff==-1,1,'first')-1+edgeTimesUnique(ij);
                        lenSecInds = tStarts:tEnds;

                        if ~isempty(lenSecInds)
                            signalFull = double(squeeze(EEG.data(iiCh1,lenSecInds)));
                            b1 = glmfit(eigVec(:,1:eigNumInd), approxMatNonZero(iiCh1,:));
                            artifactModel = zeros(1,size(approxMat,2));
                            if nonZeroChans(iiCh1)
                                artifactModel(tStartsIndsOne:tEndIndsOne) = ...
                                    b1'*[ones(size(eigVec,1),1) eigVec(:,1:eigNumInd)]'; % main components of the approximation
                            end

                            cFixed = [artifactModel zeros(1,length(c{iiCh1})-length(artifactModel))];
                            artifactRecTemp = waverec(cFixed,l{iiCh1},waveletType);
                            artifactRec = zeros(1,length(signalFull));
                            artifactRec( eyePulse(iiCh1,lenSecInds) ) = artifactRecTemp( eyePulse(iiCh1,lenSecInds) );
                            signalRec = signalFull - artifactRec;

                            pulseDiff = diff([0 eyePulse(iiCh1,lenSecInds) 0]);
                            tStarts = find(pulseDiff==1);
                            tEnds = find(pulseDiff==-1)-1;

                            unqVals = setdiff(unique(diff(eyePulse(iiCh1,lenSecInds))),0);
                            tStarts0 = []; 
                            tEnds0 = [];
                            for iii = 1:length(unqVals)
                                if unqVals(iii)==1
                                    tEnds0 = find(diff(eyePulse(iiCh1,lenSecInds))==unqVals(iii));
                                elseif unqVals(iii)==-1
                                    tStarts0 = find(diff(eyePulse(iiCh1,lenSecInds))==unqVals(iii))+1;
                                end
                            end


                            if sum(eyePulse(iiCh1,lenSecInds))~=length(signalFull) && ~isempty(tStarts) && ~isempty(tEnds)
                                if ~isempty(tStarts0) && isempty(tEnds0)
                                    tEnds0 = length(eyePulse(iiCh1,lenSecInds));
                                end
                                if isempty(tStarts0) && ~isempty(tEnds0)
                                    tStarts0 = 1;
                                end
                                if isempty(tStarts0(tStarts0 < tEnds0(1)))
                                    tStarts0 = [1 tStarts0];
                                end
                                if isempty(tEnds0(tEnds0 > tStarts0(end)))
                                    tEnds0 = [tEnds0 length(eyePulse(iiCh1,lenSecInds))];
                                end
                            end


                            lineSig = zeros(1,length(signalRec));
                            lineSig0 = zeros(1,length(signalRec));
                            for tInd = 1:length(tStarts)
                                t1 = tStarts(tInd); 
                                t2 = tEnds(tInd);
                                if ismember(edgeTimesUnique(ij), lenSecInds)
                                    signalFullNew = double(squeeze(cleanedSignal(iiCh1,lenSecInds)));
                                    lineSig(t1:t2) = signalFullNew(t1)+(signalFullNew(t2)-signalFullNew(t1))/(t2-t1)*((t1:t2)-t1);                                       
                                else
                                    lineSig(t1:t2) = signalFull(t1)+(signalFull(t2)-signalFull(t1))/(t2-t1)*((t1:t2)-t1);   
                                end
                                lineSig0(t1:t2) = signalRec(t1)+(signalRec(t2)-signalRec(t1))/(t2-t1)*((t1:t2)-t1); % bring it around zero

                            end



                            tempSig = signalRec  - lineSig0 + lineSig; 
                            cleanedSignal(iiCh1,lenSecInds) = tempSig;
                        end

                    end

            end
                           
            
    end
            
                
end



function plotDWTPCAResults(signalFull, signalRec,artifactRec, eyePulsePiece, srate)

   
    figure; subplot(211); plot(signalFull); hold on; plot(artifactRec,'m'); plot(eyePulsePiece*100,'k')
    legend('Raw','artifact','Pulse')
    subplot(212);  plot(signalFull); hold on; plot(signalRec);
    legend('Raw','reconstrcuted')
    


    [pxx, f] = pwelch(signalFull(~eyePulsePiece),[],[],[],srate);
    [pxx0, f0] = pwelch(signalFull(eyePulsePiece),[],[],[],srate);
    [pxx1, f1] = pwelch(signalFull,[],[],[],srate);
    [pxx2, f2] = pwelch(signalRec,[],[],[],srate);
    [pxx3, f3] = pwelch(artifactRec,[],[],[],srate);

    figure; plot(f0,pxx0); hold on; plot(f,pxx); plot(f1,pxx1); plot(f2,pxx2); plot(f3,pxx3);                        
    set(gca,'XLim',[0 30])
    legend('blink part','non-blink part','full epoch','reconstrcuted','artifact')


    

end
