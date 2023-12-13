%% Martin White, Kleckner Lab Harvard University April 2023

%% Function Description
%Takes Zip2/3, Hop1, and Zip1 normalized signal intensity profiles and
%calculates the ratio

%Input
%input - this is a cell (output from function normalizeSignalsBySTD).
%each row is a traced bivalent. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column is the position on the traced
%bivalent, column 2 is the signal intensity at that position

%output
%column1: Zip3/Hop1 ratio
%column2: Zip3/Zip1 ratio
%column3: Hop1/Zip1 ratio


%NB
% I added an arbitrary value of 100 to normalized signal intensity values
% as an ad hoc solution to problems caused by the presence of negative
% values.

function signalRatios = getTriadSignalRatios(input)

numOfBivs                   = length(input(:,1));
signalRatios{numOfBivs,3}    = [];

for i = 1:numOfBivs
    
    %step 1: add the means back into the normalized data
    
    %column1: Zip3/Hop1
    signalRatios{i,1}(:,1) = input{i,1}(:,1);
    signalRatios{i,1}(:,2) = (input{i,1}(:,2) + 100)./...
                                (input{i,2}(:,2) + 100);
    
    %column2: Zip3/Zip1
    signalRatios{i,2}(:,1) = input{i,1}(:,1);
    signalRatios{i,2}(:,2) = (input{i,1}(:,2) + 100)./...
                                (input{i,3}(:,2) + 100);
    
    %column3: Hop1/Zip1
    signalRatios{i,3}(:,1) = input{i,1}(:,1);
    signalRatios{i,3}(:,2) = (input{i,2}(:,2) + 100)./...
                                (input{i,3}(:,2) + 100);
    
end


end
