%% Martin White, Kleckner Lab, Harvard University, June 2023

%% Function Description
% - calculates a population average signal for each triad component

%Input
% inputSignal - this is a cell (output from function formatTriadData).
%each row is a traced chromosome. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column 1 is the position on the traced
%chromosome, column 2 is the signal intensity at that position

%Output
%normLength - as the input, except each positions on the traced chromosome
%have all been normalized to the chromosome length

%uniformSampling - same format as input, except the signal intensity
%profiles of each measured chromosome have been uniformly sampled in units
%of normalized chromosome length

%averageSignal - population average signal intensity profile of each triad
%component.  Derived from variable uniformSampling
%%


function [normLength,uniformSampling,averageSignal] = getPopnAverageSignal(inputSignal)

[numBiv,numSignals] = size(inputSignal);

normLength{numBiv,numSignals} = [];
uniformSampling{numBiv,numSignals} = [];
averageSignal{1,numSignals} = [];

%step 1: normalize bivalent to uniform length

for i = 1:numBiv
    for j = 1:numSignals
        bivLength = length(inputSignal{i,j}(:,1));
        normLength{i,j} = (0:1/(bivLength-1):1)';
        normLength{i,j}(:,2) = inputSignal{i,j}(:,2);
    end

end


%step 2: upsample or fit spline to get uniform sampling rate for each
%bivalent

%step 2.1: calculate the smallest sampling rate in the population

longestBiv = cellfun(@(x) length(x),inputSignal,'UniformOutput',false);
longestBiv = cellfun(@(x) x(1,1),longestBiv);
longestBiv = max(longestBiv(:,1));
% 
% highestSampleRate = 1/(longestBiv - 1);

%step 2: fit spline to interpolate for resampling

for i = 1:numBiv
    for j = 1:numSignals
        uniformSampling{i,j} = (0:1/(longestBiv-1):1)';
        uniformSampling{i,j}(:,2) = (spline(normLength{i,j}(:,1),...
            normLength{i,j}(:,2),0:1/(longestBiv-1):1))';
    end

end


%step 3: average signal for all bivalents

%step 3.1: for each signal, make a matrix of the uniformly sampled
%intensity values

signalMatrices{1,numSignals} = [];

for i = 1: numBiv
    for j = 1:numSignals
        signalMatrices{1,j}(i,:) = uniformSampling{i,j}(:,2)';
    end
end


for i = 1:numSignals
    averageSignal{1,i}(:,1) = (0:1/(longestBiv-1):1)';
    averageSignal{1,i}(:,2) = (mean(signalMatrices{1,i}))';
end



end
