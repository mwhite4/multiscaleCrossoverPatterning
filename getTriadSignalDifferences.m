%% Martin White, Kleckner Lab Harvard University August 2024

%% Function Description
%Takes Zip2/3, Hop1, and Zip1 normalized signal intensity profiles and
%calculates the difference for each pairwise combination of signals.

%Input
%input - this is a cell (output from function normalizeSignalsBySTD).
%each row is a traced bivalent. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column is the position on the traced
%bivalent, column 2 is the normalized signal intensity at that position

%output
%column1: Hop1 - Zip1 difference
%column2: Hop1 - Zip3 difference
%column3: Zip1 - Zip3 difference


function [signalDifferences,signalMeanDifferences] = getTriadSignalDifferences(input)

[numOfBivs,numOfSignals]                            = size(input);
signalDifferences{numOfBivs,numOfSignals}           = [];
signalMeanDifferences(1:numOfBivs,1:numOfSignals)   = nan;

for i = 1:numOfBivs
    
    %column1: Hop1 - Zip1
    signalDifferences{i,1}(:,1) = input{i,1}(:,1);
    signalDifferences{i,1}(:,2) = input{i,2}(:,2) - input{i,3}(:,2);
    
    %column2: Hop1 - Zip3
    signalDifferences{i,2}(:,1) = input{i,1}(:,1);
    signalDifferences{i,2}(:,2) = input{i,2}(:,2) - input{i,1}(:,2);
    
    %column3: Zip1 - Zip3
    signalDifferences{i,3}(:,1) = input{i,1}(:,1);
    signalDifferences{i,3}(:,2) = input{i,3}(:,2) - input{i,1}(:,2);
    
end

for i = 1:numOfBivs
    for j = 1:numOfSignals
        signalMeanDifferences(i,j) = sum(abs(signalDifferences{i,j}(:,2)))./length(signalDifferences{i,j}(:,2));
    end
end


end
