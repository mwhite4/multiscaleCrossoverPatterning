%% Martin White, Kleckner Lab Harvard University December 2023

%% Function Description
%Takes a vector of measured distances between adjacent peaks and then, for
%each measured distance, picks a number between 0 and half that distance
%at random from a uniform distribution.  It repeats this process 100 times.

%Input
%peakSpacing: a vector of measured distances between adjacent peaks (output
%from getCOSpacing)

%Output
%A vector of distances

function nullModel = getNullModelIndependentPhasing(peakSpacing)

n = length(peakSpacing);
nullModel(1:n,1) = nan;

for i = 1:n
    for j = 1:100
        nullModel(i,j) = rand(1)*(0.5*peakSpacing(i));
    end
end

nullModel = nullModel(:);

end