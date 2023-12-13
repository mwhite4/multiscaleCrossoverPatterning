
%% Martin White, Kleckner Lab November 2023

%% Function Description
% A function for calculating the distances from each detected domainal peak
% to its nearest focal peaks.

%Input

%DomainalPeaks: A cell. Row 1 Column 1 contains a matrix of Zip3 peak
%positions.  Row 1 Column 2 contains a matrix of Hop1 peak positions. Row 1
%Column 3 contains a matrix of Zip1 peak positions.  
% 
%Each matrix takes the following format:
%Each row is a separate chromosome.
%First column is the chromosome length, Additional columns have the
%positions of detected crossovers (in the same units as chromosome length).
%empty cells should be filled with NaNs

%FocalPeaks: same format as DomainalPeaks, but data is detected focal peaks


%Output
%DomainaltoFocalPeakSpacing: A Cell.  

%first row of output is data for each traced bivalent.
%second row of output is cumulated data for all traced bivalents

%first column is closest distances between Zip3 domainal and focal peaks
%second column isclosest distances between Hop1 domainal and focal peaks
%third column isclosest distances between Zip1 domainal and focal peaks

function DomainaltoFocalPeakSpacing = getDomainalToFocalPeakDistances(Domainalpeaks,Focalpeaks)

%Step 1: set up output cell
DomainaltoFocalPeakSpacing{2,3} = [];

[numOfBivs,n] = size(Domainalpeaks{1,1});

DomainaltoFocalPeakSpacing{1,1}(1:numOfBivs,1:(n-1)*2) = nan;
DomainaltoFocalPeakSpacing{1,2}(1:numOfBivs,1:(n-1)*2) = nan;
DomainaltoFocalPeakSpacing{1,3}(1:numOfBivs,1:(n-1)*2) = nan;

%Step 2: calculate minimum distances between peaks of two types for each
%traced bivalent

for i = 1:numOfBivs
    Zip3SpacingMatrix = sort(abs(Domainalpeaks{1,1}(i,2:end) - Focalpeaks{1,1}(i,2:end).'));
    DomainaltoFocalPeakSpacing{1,1}(i,:) = horzcat(Zip3SpacingMatrix(1,:),Zip3SpacingMatrix(2,:));
    Hop1SpacingMatrix = sort(abs(Domainalpeaks{1,2}(i,2:end) - Focalpeaks{1,2}(i,2:end).'));
    DomainaltoFocalPeakSpacing{1,2}(i,:) = horzcat(Hop1SpacingMatrix(1,:),Hop1SpacingMatrix(2,:));
    Zip1SpacingMatrix = sort(abs(Domainalpeaks{1,3}(i,2:end) - Focalpeaks{1,3}(i,2:end).'));
    DomainaltoFocalPeakSpacing{1,3}(i,:) = horzcat(Zip1SpacingMatrix(1,:),Zip1SpacingMatrix(2,:));

end

 %Step 3: for each pairwise comparison, aggregate all the distances from
 %each traced bivalent and remove nans
 
 for i = 1:3
     DomainaltoFocalPeakSpacing{2,i} = sort(DomainaltoFocalPeakSpacing{1,i}(:));
     DomainaltoFocalPeakSpacing{2,i} = DomainaltoFocalPeakSpacing{2,i}(~isnan(DomainaltoFocalPeakSpacing{2,i}));
 end
 

end
