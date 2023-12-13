function [intensityProfiles,zip3Traces,hop1Traces,zip1Traces] = formatTriadData(inputDirectory,extension)

%% Martin White, Kleckner Lab. January 2023

%% Function Description
% - Takes a folder containing intensity profiles of spread, S. cerevisiae
% pachytene chromosomes stained for Zip2, Hop1 and Zip1, import the data into
% matlab and formats as a cell.  

% Pixel size was set to 67 nm

%Input

%inputDirectory
%user defines the folder containing the raw data from Beth Weiner

%extension
% user needs to define the file extension

%Output

%triadIntensityProfiles
%This is a cell.  Each row is a traced bivalent. Column 1 is the Zip3 data,
%column 2 is the Hop1 data and column 3 is the Zip1 data.  Each data is a 
%matrix, column 1 is the position on the traced bivalent, column 2 is the 
%signal intensity at that position

%zip3Traces, hop1Traces and zip1Traces
%These are lists of imported data files. It can be useful if you need to 
%match the data in the Cell format to the raw data file

%%
% Step 1: Define Parameters

%1.1: File parameters : turn these into input variables
zip3Identifier      = 'zip3';                                               %This is the unique identifier in the file name that states the nature of the signal
hop1Identifier      = 'hop1';                                               %This is the unique identifier in the file name that states the nature of the signal 
zip1Identifier      = 'zip1';                                               %This is the unique identifier in the file name that states the nature of the signal

%Pixel to micrometer conversion factor
pixelSize           = 0.067;


%Step 2: get files
%Step 2.1 select folder from interactive window
% inputDirectory      = uigetdir()

%Step 2.2: get list of names of zip3 files
zip3Identifier      = strcat('*',zip3Identifier,'*',extension);
zip3Traces          = dir(fullfile(inputDirectory,zip3Identifier));

%Step 2.3: get list of names of hop1 files
hop1Identifier      = strcat('*',hop1Identifier,'*',extension);
hop1Traces          = dir(fullfile(inputDirectory,hop1Identifier));

%Step 2.4 get list of names of zip1 files
zip1Identifier      = strcat('*',zip1Identifier,'*',extension);
zip1Traces          = dir(fullfile(inputDirectory,zip1Identifier));


%Step 2.5: check that there is the same number of traces for all files.  If
%not end program and write message

if length(zip1Traces) ~= length(zip3Traces) || length(zip1Traces) ~= length(hop1Traces)
    disp('Error: One or more dataset is incomplete.  Program Ending')
    return
end
    
intensityProfiles{length(zip1Traces),3}         = [];


%Step 3: Analyze the Traces.  Process each traces, one at a time and save output

for dataset = 1:length(zip1Traces)
    
    %Step 3.1: load the zip3, hop1 and zip1 data for current bivalent
    
    %Step 3.1.1: load zip3 trace
    zip3TraceName       = zip3Traces(dataset).name;
    pathToZip3Trace     = fullfile(inputDirectory,zip3TraceName);
    zip3Trace           = readmatrix(pathToZip3Trace);
    zip3Trace(:,all(isnan(zip3Trace),1)) = [];                              %remove columns of NaNs
    zip3Trace(all(isnan(zip3Trace),2),:) = [];                              %remove rows of NaNs
    zip3Trace(:,1)      = zip3Trace(:,1).*pixelSize;
    
    intensityProfiles{dataset,1}(:,1) = zip3Trace(:,1);
    intensityProfiles{dataset,1}(:,2) = zip3Trace(:,2);
    

    %Step 3.1.2: load hop1 trace
    hop1TraceName       = hop1Traces(dataset).name;
    pathToHop1Trace     = fullfile(inputDirectory,hop1TraceName);
    hop1Trace           = readmatrix(pathToHop1Trace);
    hop1Trace(:,all(isnan(hop1Trace),1)) = [];                              %remove columns of NaNs
    hop1Trace(all(isnan(hop1Trace),2),:) = [];                              %remove rows of NaNs
    hop1Trace(:,1)      = hop1Trace(:,1).*pixelSize;
    
    intensityProfiles{dataset,2}(:,1) = hop1Trace(:,1);
    intensityProfiles{dataset,2}(:,2) = hop1Trace(:,2);
    
    
    %Step 3.1.3: load zip1 trace
    zip1TraceName       = zip1Traces(dataset).name;
    pathToZip1Trace     = fullfile(inputDirectory,zip1TraceName);
    zip1Trace           = readmatrix(pathToZip1Trace);
    zip1Trace(:,all(isnan(zip1Trace),1)) = [];                              %remove columns of NaNs
    zip1Trace(all(isnan(zip1Trace),2),:) = [];                              %remove rows of NaNs
    zip1Trace(:,1)      = zip1Trace(:,1).*pixelSize;
  
    intensityProfiles{dataset,3}(:,1) = zip1Trace(:,1);
    intensityProfiles{dataset,3}(:,2) = zip1Trace(:,2);
    
    
end


end


