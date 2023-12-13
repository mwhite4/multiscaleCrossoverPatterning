%% Martin White, Kleckner Lab January 2023

%% Function Description
% takes each measured chromosome, and normalizes its Zip3/Zip2, Hop1 and
% Zip1 signal intensity values by subtracting the mean and dividing by the
% standard deviation.

%Input
%intensityProfiles - this is a cell (output from function formatTriadData).
%each row is a traced bivalent. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column is the position on the traced
%bivalent, column 2 is the signal intensity at that position

%Output
%normSignals - the same format as the input, but the original signal has
%been normalized by subtracting the mean and dividing by the standard
%deviation of that signal
%%
function normSignals = normalizeSignalsBySTD(intensityProfiles)

numOfBivs                   = length(intensityProfiles(:,1));
normSignals{numOfBivs,3}    = [];

for i = 1:numOfBivs  
    for j = 1:3
        signal                  = intensityProfiles{i,j}(:,2);
        signalMinusMean         = signal - mean(signal);
        stdDev                  = std(signalMinusMean);
        normSignals{i,j}(:,1)   = intensityProfiles{i,j}(:,1);
        normSignals{i,j}(:,2)   = signalMinusMean./stdDev;
    end
end

end