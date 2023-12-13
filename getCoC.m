%% Martin White, Kleckner Lab January 2019

%% This is a program for generating 'high resolution' CoC curves for
%chromosomes within a user-defined range of lengths.  The program:

%1) takes the input, normalizes each chromosome to unit length and then
% seperates each chromosome into a user-defined number of equally sized
% intervals

%2) calculates the coefficient of coincidence for each pair of intervals
% 
%3) calculates the mean coefficient of coincidence for each inter-interval
% distance.

%input
%inputData: this should be a matrix, each row is a separate chromosome.
%First column is the chromosome length, Additional columns have the
%positions of detected crossovers (in the same units as chromosome length).
%empty cells should be filled with NaNs

%numOfintervals: the total number of equally sized intervals that each
%chromosome will be divided into

%output
%output: a matrix.  Column 1 are inter-interval distances, column 2
%contains associated CoC values.


%%

function output = getCoC(inputData,numOfIntervals)

%Step 1: format input data
inputData(inputData==0)=NaN;                                                %convert 0s to NaNs
inputData   = inputData(~all(isnan(inputData),2),:);                        %remove any rows and columns composed only of NaNs
inputData   = inputData(:,~all(isnan(inputData)));


%Step 3.1: Normalize CO positions relative to total chromosome length
normCO = inputData(:,2:end)./inputData(:,1);
[totChroms,~]   = size(normCO);


%Step 3.2: Generate CoC data for each total number of intervals

%Step 3.2.1: Calculate Interval Boundaries and labels
intSize     = 1/numOfIntervals;
edges       = 0:intSize:1;

%Step 3.2.2: Calculate the CO frequency for each interval
obsCOfreq  = histcounts(normCO,edges);
obsCOfreq  = obsCOfreq./totChroms;

%Step 3.2.3: For each pair of intervals, calculate the observed
%frequency of chromosomes with a(t least one) CO in both intervals
COcounts(1:totChroms,1:numOfIntervals) = nan;
for chrom = 1:totChroms
    COcounts(chrom,:) = histcounts(normCO(chrom,:),edges);
end
max1COcounts    = COcounts;
max1COcounts(max1COcounts>1) = 1;
for j=1:numOfIntervals-1
    for k=j+1:numOfIntervals
        obs2COfreq(j,k-j)=(sum(max1COcounts(:,j).*...                %each row is an interval, each column is an inter-interval distance
            max1COcounts(:,k)))./totChroms;
    end
end

%Step 3.2.4: For each pair of intervals, calculate the expected
%frequency of double COs, assuming indepedence
for j = 1:numOfIntervals-1
    for k = j+1:numOfIntervals
        exp2COfreq(j,k-j)  = obsCOfreq(j)*obsCOfreq(k);        %each row is an interval, each column is an inter-interval distance
        
    end
end

%Step 3.2.5: For each pair of intervals, calculate the Coefficient of
%Coincidence
CoC = obs2COfreq./exp2COfreq;

%Step 3.2.6: Calculate average CoC for all inter-interval distances
meanCoCperDistance   = mean(CoC,1,'omitnan');

%Step 3.2.7: Calculate inter-interval distances
interInterval           = 1:numOfIntervals-1;
relInterInterval        = 1/numOfIntervals*interInterval;
meanLength              = mean(inputData(:,1));
interIntervalDist       = relInterInterval*meanLength;

output(:,1) = interIntervalDist;
output(:,2) = meanCoCperDistance;

end





