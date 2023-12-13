%% Martin White. Kleckner Lab. Nov 2023

%% Script Description
%
% Script for simulating two-tiered crossover patterning by the Beam-Film
% Model.
% 
% Logic of simulation is based on that used to simulate SC nucleations in
% Sordaria.
% 
% Do 'canonical' BF simulation, but with a low Smax to get minority
% Crossovers.

% Repeat the simulation, using the same precursor array (with its
% associated sensitivity values) but use a higher Smax to get minority +
% majority crossovers.

% Subtract minority crossovers (first simulation) from minority + majority
% crossovers (second simulation) to get majority crossovers.

n               = 10000;
N               = 13;
E               = 0.75;
B               = 1;
M               = 1;
A               = 1;
L               = 0.1;
Bs              = 1;
Be              = 1;
Bd              = 1;
Y               = 1;
T2prob          = 0;
Bsmax           = 1;
cL              = 0.85;
cR              = 0.85;
SmaxMinority    = 1.5;
SmaxMajority    = 5;
chrXVlength     = 3.2;

%get fixed precursor array
[totalPrecursor]                        = event_per_object_from_population_mean(n,N,B);
[precursor_positions,totalPrecursor]    = generate_precursor_array(n,totalPrecursor,E,Bs,Be,Bd);
[precursor_sensitivity_matrix]          = generate_precursor_sensitivities(n,totalPrecursor,A);


%start at low Smax and get minority crossovers
maximum_stress_per_object_minority      = event_per_object_from_population_mean(n,SmaxMinority,Bsmax);
minority_COs                            = pattern_event_designations_according_to_beam_film_model(n,totalPrecursor,precursor_positions,precursor_sensitivity_matrix,maximum_stress_per_object_minority,L,cL,cR);

minority_COs                        = minority_COs- precursor_positions;
minority_COs                        = 1 - minority_COs;
minority_COs(minority_COs==1)       = 0;
minority_COs(minority_COs==0)       = NaN;
minority_COs                        = sort(minority_COs,2);
object_length(1:n,1)                = 1;
simCOs_minority                     = horzcat(object_length,minority_COs);
simCOs_minority                     = simCOs_minority.*chrXVlength;

simCOs_classII_CoC30                = getCoC(simCOs_minority,30);
simCOs_classII_spacing              = getCOSpacing(simCOs_minority);

%increase Smax to add majority crossovers
maximum_stress_per_object_majority  = event_per_object_from_population_mean(n,SmaxMajority,Bsmax);
majorityAndMinority_COs             = pattern_event_designations_according_to_beam_film_model(n,totalPrecursor,precursor_positions,precursor_sensitivity_matrix,maximum_stress_per_object_majority,L,cL,cR);

majorityAndMinority_COs                             = majorityAndMinority_COs- precursor_positions;
majorityAndMinority_COs                             = 1 - majorityAndMinority_COs;
majorityAndMinority_COs(majorityAndMinority_COs==1) = 0;
majorityAndMinority_COs(majorityAndMinority_COs==0) = NaN;
majorityAndMinority_COs                             = sort(majorityAndMinority_COs,2);
object_length(1:n,1)                                = 1;
simCOs_majorityAndMinority                          = horzcat(object_length,majorityAndMinority_COs);
simCOs_majorityAndMinority                          = simCOs_majorityAndMinority.*chrXVlength;

simCOs_majorityAndMinority_CoC30                    = getCoC(simCOs_majorityAndMinority,30);
simCOs_majorityAndMinority_spacing                  = getCOSpacing(simCOs_majorityAndMinority);

%remove minority crossvoers to get majority crossovers

simCOs_majority(1:n,1:length(simCOs_majorityAndMinority(1,:)))  = nan;
total_minorityCOs                                               = countTotalCOs(simCOs_minority);

simCOs_majority(:,1) = chrXVlength;

for i = 1:n
    simCOs_majority(i,2:end-total_minorityCOs(i)) = setdiff(simCOs_majorityAndMinority(i,2:end),simCOs_minority(i,2:end));
end


simCOs_majority_CoC30        = getCoC(simCOs_majority,30);
simCOs_majority_spacing      = getCOSpacing(simCOs_majority);
