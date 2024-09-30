%% Martin White, Kleckner Lab August 2024

%% Function Description
%detect peaks in signals using the findpeaks function of MATLAB.  
%detect 'humps' in signal (indicating focal/domainal peak/Gaussian that has
%merged with a neighboring peak) by using the findpeaks function of MATLAB
%on the first derivative of the signal (humps are indicated by changes in
%gradient).

%Input
% inputCell - this is a cell (output from function formatTriadData).
%each row is a traced chromosome. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column 1 is the position on the traced
%chromosome, column 2 is the signal intensity at that position

%Output
% peakAndHumpPositions -  a cell.  first column is the output for Zip3
% signals, the second column is the output for Hop1 signals and the third
% column is the output for Zip1 signals.  The data for each signal is
% organized as a matrix.  Each measured chromosome is on a separate row.
% The first column is the measured chromosome length.  The susequent
% columns are the measured peak/hump positions

function peakAndHumpPositions = getSignalPeakAndHumpPositions(inputCell)

[totChroms,totSignals]                          = size(inputCell);
peakPositions{totChroms,totSignals}             = [];
gradient_peakPositions{totChroms,totSignals}    = [];
gradient_valleyPositions{totChroms,totSignals}  = [];
peakAndHumpPositionsCell{totChroms,totSignals}  = [];


for i = 1:totChroms
    for j = 1:totSignals

        %find peaks only
        [~,peakPositions{i,j}] = findpeaks(inputCell{i,j}(:,2),inputCell{i,j}(:,1));

        %find peaks on the gradient
        [~,gradient_peakPositions{i,j}] = findpeaks(gradient(inputCell{i,j}(:,2)),inputCell{i,j}(:,1));

        %find valleys on the gradient
        [~,gradient_valleyPositions{i,j}] = findpeaks(-gradient(inputCell{i,j}(:,2)),inputCell{i,j}(:,1));

        %if no evidence of humps, use the peak positions

        if length(peakPositions{i,j}) >= length(gradient_peakPositions{i,j}) && length(peakPositions{i,j}) >= length(gradient_valleyPositions{i,j})
            peakAndHumpPositionsCell{i,j} = peakPositions{i,j};

            %else do peak and hump detection
        else
            if length(gradient_peakPositions{i,j}) < length(gradient_valleyPositions{i,j})
                gradient_peakPositions{i,j} = vertcat(inputCell{i,j}(1,1),gradient_peakPositions{i,j}(:,1));
            end

            if length(gradient_valleyPositions{i,j}) < length(gradient_peakPositions{i,j})
                gradient_valleyPositions{i,j} = vertcat(gradient_valleyPositions{i,j},inputCell{i,j}(end,1));
            end

            if length(gradient_peakPositions{i,j}(:,1)) ==  length(gradient_valleyPositions{i,j}(:,1))
                gradient_peakPositions{i,j}(:,2)    = gradient_valleyPositions{i,j}(:,1);
                peakAndHumpPositionsCell{i,j}       = mean(gradient_peakPositions{i,j},2);
            end

        end

    end
end

% Format peak positions for downstream analyses

numbPeaks   = cellfun(@(x) numel(x), peakAndHumpPositionsCell);
max_Peaks   = max(numbPeaks(:));

peakAndHumpPositions{1,totSignals} = [];
for i = 1:totChroms
    for j = 1:totSignals
        peakAndHumpPositions{1,j}(i,1:max_Peaks+1)  = nan;
        peakAndHumpPositions{1,j}(i,1)              = inputCell{i,j}(end,1);           %bivalent length
        peakAndHumpPositions{1,j}(i,2:length(peakAndHumpPositionsCell{i,j})+1)...
                                                    = peakAndHumpPositionsCell{i,j};
    end
end

end
