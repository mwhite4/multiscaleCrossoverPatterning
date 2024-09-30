%% Martin White, Kleckner Lab January 2023

%revised August 2024, amplitudes are no longered squared to give the PSD

%% Function Description
% for triad analysis. Pixel size was set to 67 nm

%Input
%intensityProfiles - this is a cell (output from function formatTriadData).
%each row is a traced bivalent. column1 is the Zip3 data, column 2 is the 
%Hop1 data and column 3 is the Zip1 data.  
%Each data is a matrix, column is the position on the traced
%bivalent, column 2 is the signal intensity at that position

%Output
%Fhats - the same format as the input, but for each intensity profile of
%each traced bivalent, there is a single column containing the fourier 
%coefficients of that intensity profile


%Amplitudes - the same format as the input, but for each signal, of each traced
%bivalent, column 1 is the frequency in units of cycles per micron and
%column 2 is the amplitdue of the (co-)sine wave corresponding to that frequency

%oneSidedAmplitude - the signal is real valued so take the positive frequencies

%%


function [Fhats,Amplitudes,oneSidedAmplitudes] = triadFrequencyDomain(intensityProfiles)

pixelSize                       = 0.067;
nyquistFreq                     = 1/(2*pixelSize);
numOfBivalents                  = length(intensityProfiles(:,1));

Fhats{numOfBivalents,3}         = [];
Amplitudes{numOfBivalents,3}          = [];
oneSidedAmplitudes{numOfBivalents,3}  = [];

for i = 1:numOfBivalents
    n                           = length(intensityProfiles{i,1}(:,1));
    dX                          = pixelSize;
    fs                          = 1/dX;
    frequencies                 = (0:n-1)*(fs/n);
        
    for j = 1:3
        Fhats{i,j} = fft(intensityProfiles{i,j}(:,2));
        
        Amplitudes{i,j}(:,1) = frequencies;
        % PSDs{i,j}(:,2) = abs(Fhats{i,j}).^2/n;
        Amplitudes{i,j}(:,2) = abs(Fhats{i,j}./n);
        
        frequencies_oneSided            = Amplitudes{i,j}(1:floor(n/2),1)';
        nyqRange                        = frequencies_oneSided < nyquistFreq;
        frequencies_oneSided(~any(nyqRange,2)) = [];
        
        PSD_oneSided                    = Amplitudes{i,j}(1:floor(n/2),2);
        PSD_oneSided(2:end-1)           = 2*PSD_oneSided(2:end-1);
        PSD_oneSided(~any(nyqRange,2))  = [];
        
        oneSidedAmplitudes{i,j}(:,1)          = frequencies_oneSided;
        oneSidedAmplitudes{i,j}(:,2)          = PSD_oneSided;
        
    end
    
    
end



end
