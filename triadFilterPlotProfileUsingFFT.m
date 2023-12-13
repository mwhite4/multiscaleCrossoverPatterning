
%% Martin White, Kleckner Lab January 2023

%% Function Description
% for triad analysis. A Fourier approach for separating intensity profile
% signals into low and high frequency components.  Creates two new
% intensity profiles by removing the DC component and either filtering out
% frequencies above (lfIntensityProfiles), or below (hfIntensityProfiles),
% a user defined threshold (maxPeriod)


%Input
% intensityProfiles - this is a cell.  each row is a traced bivalent.
%column1 is the Zip3 data, column 2 is the Hop1 data and column 3 is the
%Zip1 data.  Each data is a matrix, column is the position on the traced
%bivalent, column 2 is the signal intensity at that position.  This is the
%output of function formatTriadData

% Fhats - the same format as intensityProfiles, but for each intensity 
%profile of each traced bivalent, there is a single column containing the 
%fourier coefficients of that intensity profile.  This is one of the
%outputs of function triadFrequencyDomain

% PSDs - the same format as intensityProfiles, but for each signal, of 
%each traced bivalent, column 1 is the frequency in units of cycles per 
%micron and column 2 is the power of that frequency.  This is one of the
%outputs of function triadFrequencyDomain

% maxPeriod - the threshold frequency/periodicity expressed as micrometers
% per cycle (i.e. periodicity, or wavelength, rather than frequency)
% for separating the signal into low/high frequency components

% minPower - the minimum power that a frequency must have in order to not
% be filtered out.  This is intended as a noise filter


%Output
% lfIntensityProfiles - the intensity profiles with every frequency higher
% than maxPeriod and with less power than minPower, filtered out

% hfIntensityProfiles - the intensity profiles with every frequency lower
% or equal to maxPeriod and with less power than minPower, filtered out

%%

function [lfIntensityProfiles,hfIntensityProfiles] = triadFilterPlotProfileUsingFFT (intensityProfiles,inputFhat,inputPSD,maxPeriodicity,minPower)

maxFreq = 1/maxPeriodicity;
lfIntensityProfiles{length(inputFhat(:,1)),3} = [];
hfIntensityProfiles{length(inputFhat(:,1)),3} = [];


for i = 1:length(inputFhat(:,1))
    
    lfIntensityProfiles{i,1}(:,1) = intensityProfiles{i,1}(:,1);
    lfIntensityProfiles{i,2}(:,1) = intensityProfiles{i,1}(:,1);
    lfIntensityProfiles{i,3}(:,1) = intensityProfiles{i,1}(:,1);
    
    hfIntensityProfiles{i,1}(:,1) = intensityProfiles{i,1}(:,1);
    hfIntensityProfiles{i,2}(:,1) = intensityProfiles{i,1}(:,1);
    hfIntensityProfiles{i,3}(:,1) = intensityProfiles{i,1}(:,1);
    
    
    %get the indices of the frequencies that you want to keep
    zip3Indices_hf = inputPSD{i,1}(:,1) >= maxFreq;
    hop1Indices_hf = inputPSD{i,2}(:,1) >= maxFreq;
    zip1Indices_hf = inputPSD{i,3}(:,1) >= maxFreq;
    
    %this will only do the first half
    %get the length of zeros then remove this (-1 for DC current) from the
    %end
    
    NegFreq                             = sum(zip3Indices_hf(:)==0) - 1;
    zip3Indices_hf(end-NegFreq+1:end,1) = 0;
    hop1Indices_hf(end-NegFreq+1:end,1) = 0;
    zip1Indices_hf(end-NegFreq+1:end,1) = 0;  
   
    %Invert the 'hf' indices to get the 'lf' indices
    zip3Indices_lf = 1 - zip3Indices_hf;
    hop1Indices_lf = 1 - hop1Indices_hf;
    zip1Indices_lf = 1 - zip1Indices_hf;
    
    %remove the DC current
    zip3Indices_lf(1) = 0;
    hop1Indices_lf(1) = 0;
    zip1Indices_lf(1) = 0;
    
    zip3Indices_hf(1) = 0;
    hop1Indices_hf(1) = 0;
    zip1Indices_hf(1) = 0;
    
    %Now remove the indices of any frequency whose power is below the set
    %threshold
    zip3LowPowerIndices = inputPSD{i,1}(:,2) >= minPower(1);
    hop1LowPowerIndices = inputPSD{i,2}(:,2) >= minPower(2);
    zip1LowPowerIndices = inputPSD{i,3}(:,2) >= minPower(3);
    
    %combine the two thresholds
    zip3Indices_hf = zip3Indices_hf.*zip3LowPowerIndices;
    zip3Indices_lf = zip3Indices_lf.*zip3LowPowerIndices;
    
    hop1Indices_hf = hop1Indices_hf.*hop1LowPowerIndices;
    hop1Indices_lf = hop1Indices_lf.*hop1LowPowerIndices;
    
    zip1Indices_hf = zip1Indices_hf.*zip1LowPowerIndices;
    zip1Indices_lf = zip1Indices_lf.*zip1LowPowerIndices;
    
    
    %zero out the power of the frequencies that you want to remove
    zip3Fhat_filt_hf        = zip3Indices_hf.*inputFhat{i,1}(:,1);
    zip3Fhat_filt_lf        = zip3Indices_lf.*inputFhat{i,1}(:,1);
    
    hop1Fhat_filt_hf        = hop1Indices_hf.*inputFhat{i,2}(:,1);
    hop1Fhat_filt_lf        = hop1Indices_lf.*inputFhat{i,2}(:,1);
    
    zip1Fhat_filt_hf        = zip1Indices_hf.*inputFhat{i,3}(:,1);
    zip1Fhat_filt_lf        = zip1Indices_lf.*inputFhat{i,3}(:,1);
    
    
    %perform an inverse fast fourier transform on the filtered frequencies
    hfIntensityProfiles{i,1}(:,2) = ifft(zip3Fhat_filt_hf);
    lfIntensityProfiles{i,1}(:,2) = ifft(zip3Fhat_filt_lf);
    
    hfIntensityProfiles{i,2}(:,2) = ifft(hop1Fhat_filt_hf);
    lfIntensityProfiles{i,2}(:,2) = ifft(hop1Fhat_filt_lf);
    
    hfIntensityProfiles{i,3}(:,2) = ifft(zip1Fhat_filt_hf);
    lfIntensityProfiles{i,3}(:,2) = ifft(zip1Fhat_filt_lf);
      
end


end

