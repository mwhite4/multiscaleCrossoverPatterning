%% Martin White, Kleckner Lab November 2023

%% Function Description
% A function for taking a list of detected peaks of one kind (Zip3/Zip2,
% Hop1, or Zip1) and measuring the distance to its nearest neighbor in a
% separate list containing detected peaks of a different kind (Zip3/Zip2,
% Hop1, or Zip1)
% 

%Input
%input: A cell. Row 1 Column 1 contains a matrix of Zip3 peak
%positions.  Row 1 Column 2 contains a matrix of Hop1 peak positions. Row 1
%Column 3 contains a matrix of Zip1 peak positions.  
% 
%Each matrix takes the following format:
%Each row is a separate chromosome.
%First column is the chromosome length, Additional columns have the
%positions of detected crossovers (in the same units as chromosome length).
%empty cells should be filled with NaNs


%Output
%closestPeakDist: A Cell.  

%first row of output is data for each traced bivalent.
%second row of output is cumulated data for all traced bivalents

%first column is closest distances between Zip3 domainal and focal peaks
%second column isclosest distances between Hop1 domainal and focal peaks
%third column isclosest distances between Zip1 domainal and focal peaks

%first row of output is data for each traced bivalent
%second row of output is cumulated data for all traced bivalent
%first column is closest distances between Zip3 and Hop1 peaks
%second column is closest distances between Zip3 and Zip1 peaks
%third column isclosest distances between Hop1 and Zip1 peaks

function closestPeakDist = getDistanceToClosestPeak(input)

%Step 1: set up output cell
closestPeakDist{2,3} = [];

[numOfBivs,n] = size(input{1,1});

closestPeakDist{1,1}(1:numOfBivs,1:n-1) = nan;
closestPeakDist{1,2}(1:numOfBivs,1:n-1) = nan;
closestPeakDist{1,3}(1:numOfBivs,1:n-1) = nan;

%Step 2: calculate minimum distances between peaks of two types for each
%traced bivalent

for i = 1:numOfBivs
    [closestPeakDist{1,1}(i,:),~] = min(abs(input{1,1}(i,2:end) - input{1,2}(i,2:end).'));
    [closestPeakDist{1,2}(i,:),~] = min(abs(input{1,1}(i,2:end) - input{1,3}(i,2:end).'));
    [closestPeakDist{1,3}(i,:),~] = min(abs(input{1,2}(i,2:end) - input{1,3}(i,2:end).'));
end

 %Step 3: for each pairwise comparison, aggregate all the distances from
 %each traced bivalent and remove nans
 
 for i = 1:3
     closestPeakDist{2,i} = sort(closestPeakDist{1,i}(:));
     closestPeakDist{2,i} = closestPeakDist{2,i}(~isnan(closestPeakDist{2,i}));
 end
 



end
