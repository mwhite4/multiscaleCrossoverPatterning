%% Martin White, Kleckner Lab January 2019

%% Function Description
% A function for calculating the distances between adjacent crossovers
% along a population of chromosomes


%input
%inputData: this should be a matrix, each row is a separate chromosome.
%First column is the chromosome length, Additional columns have the
%positions of detected crossovers (in the same units as chromosome length).
%empty cells should be filled with NaNs

%output
%output: a vector of distances between adjacent crossovers.  Units are the
%same as that used in the input.


function [COspacing]=getCOSpacing(input)

[rows,columns] = size(input);
COspacing(1:rows,1:columns-2) = nan;


for i = 1:rows
    COspacing(i,:) = input(i,3:end) - input(i,2:end-1);
end

COspacing = COspacing(:);
COspacing = COspacing(~isnan(COspacing));
    

end
