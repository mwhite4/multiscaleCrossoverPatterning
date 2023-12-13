%% Martin White, Kleckner Lab. November 2023

%% Function Description
% accessory function for Beam-Film simulation software that allows for a
% second  patterning process using the precursors (and their associated
% sensitivity values) that remain after an initial patterning process

%source of function inspiration:
% from https://www.mathworks.com/matlabcentral/answers/386380-move-all-nan-to-end-of-matrix-columns

function sortedprecursorSensitivites = sortMinorityCO_precSensitiv_retained(precursor_sensitivity_rnd2_retain)

unsortedMatrix = precursor_sensitivity_rnd2_retain;
unsortedMatrix =unsortedMatrix';

[~,idr] = sort(isnan(unsortedMatrix),1);

S = size(unsortedMatrix);
[~,id2] = ndgrid(1:S(1),1:S(2));
idx = sub2ind(S,idr,id2);


sortedprecursorSensitivites = unsortedMatrix(idx);

sortedprecursorSensitivites = sortedprecursorSensitivites';


end
