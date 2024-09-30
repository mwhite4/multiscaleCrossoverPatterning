
%% Martin White, Kleckner Lab November 2023
%% edited August 2024
    %edit deals with rare occurances when domainal peaks are not flanked by
    %two focal peaks

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
%second column is closest distances between Hop1 domainal and focal peaks
%third column is closest distances between Zip1 domainal and focal peaks



function DomtoFocPeakSpacing = getDomainalToFocalPeakDistances(domainalPeakPositions,focalPeakPositions)

%Step 1: set up output cell
DomtoFocPeakSpacing{2,3} = [];

[numOfBivs,n] = size(domainalPeakPositions{1,1});

DomtoFocPeakSpacing{1,1}(1:numOfBivs,1:(n-1)*2) = nan;
DomtoFocPeakSpacing{1,2}(1:numOfBivs,1:(n-1)*2) = nan;
DomtoFocPeakSpacing{1,3}(1:numOfBivs,1:(n-1)*2) = nan;

%Step 2: calculate minimum distances between peaks of two types for each
%traced bivalent

for i = 1:numOfBivs

    %Zip3
    Zip3SpacingMatrix = sort(abs(domainalPeakPositions{1,1}(i,2:end) - focalPeakPositions{1,1}(i,2:end).'));

    if domainalPeakPositions{1,1}(i,2) <= focalPeakPositions{1,1}(i,2)
        Zip3SpacingMatrix(2,1) = nan;
    end

    if domainalPeakPositions{1,1}(i,sum(~isnan(domainalPeakPositions{1,1}(i,1:end)))) >= focalPeakPositions{1,1}(i,sum(~isnan(focalPeakPositions{1,1}(i,1:end))))
        Zip3SpacingMatrix(2,sum(~isnan(Zip3SpacingMatrix(2,1:end)))) = nan;
    end

    DomtoFocPeakSpacing{1,1}(i,:) = horzcat(Zip3SpacingMatrix(1,:),Zip3SpacingMatrix(2,:));

    %Hop1
    Hop1SpacingMatrix = sort(abs(domainalPeakPositions{1,2}(i,2:end) - focalPeakPositions{1,2}(i,2:end).'));

    if domainalPeakPositions{1,2}(i,2) <= focalPeakPositions{1,2}(i,2)
        Hop1SpacingMatrix(2,1) = nan;
    end

    if domainalPeakPositions{1,2}(i,sum(~isnan(domainalPeakPositions{1,2}(i,1:end)))) >= focalPeakPositions{1,2}(i,sum(~isnan(focalPeakPositions{1,2}(i,1:end))))
        Hop1SpacingMatrix(2,sum(~isnan(Hop1SpacingMatrix(2,1:end)))) = nan;
    end

    DomtoFocPeakSpacing{1,2}(i,:) = horzcat(Hop1SpacingMatrix(1,:),Hop1SpacingMatrix(2,:));


    %Zip1
    Zip1SpacingMatrix = sort(abs(domainalPeakPositions{1,3}(i,2:end) - focalPeakPositions{1,3}(i,2:end).'));

    if domainalPeakPositions{1,3}(i,2) <= focalPeakPositions{1,3}(i,2)
        Zip1SpacingMatrix(2,1) = nan;
    end

    if domainalPeakPositions{1,3}(i,sum(~isnan(domainalPeakPositions{1,3}(i,1:end)))) >= focalPeakPositions{1,3}(i,sum(~isnan(focalPeakPositions{1,3}(i,1:end))))
        Zip1SpacingMatrix(2,sum(~isnan(Zip1SpacingMatrix(2,1:end)))) = nan;
    end

    DomtoFocPeakSpacing{1,3}(i,:) = horzcat(Zip1SpacingMatrix(1,:),Zip1SpacingMatrix(2,:));

end

%Step 3: for each pairwise comparison, aggregate all the distances from
%each traced bivalent and remove nans

for i = 1:3
    DomtoFocPeakSpacing{2,i} = sort(DomtoFocPeakSpacing{1,i}(:));
    DomtoFocPeakSpacing{2,i} = DomtoFocPeakSpacing{2,i}(~isnan(DomtoFocPeakSpacing{2,i}));
end




end
