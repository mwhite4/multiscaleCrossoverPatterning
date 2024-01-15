%% Martin White. Kleckner Lab. Nov 2023

%% Script Description
% do standard modelling of majority canonical and minority crossovers (aka
% 'typeII events')

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
T2prob          = 0.2;
Bsmax           = 1;
cL              = 0.85;
cR              = 0.85;
Smax            = 3.8;
chrXVlength     = 3.2;

%get fixed precursor array
[totalPrecursor]                        = event_per_object_from_population_mean(n,N,B);
[precursor_positions,totalPrecursor]    = generate_precursor_array(n,totalPrecursor,E,Bs,Be,Bd);
[precursor_sensitivity_matrix]          = generate_precursor_sensitivities(n,totalPrecursor,A);
maximum_stress_per_object               = event_per_object_from_population_mean(n,Smax,Bsmax);

%get canonical majority crossovers
canonical_COs                       = pattern_event_designations_according_to_beam_film_model(n,totalPrecursor,precursor_positions,precursor_sensitivity_matrix,maximum_stress_per_object,L,cL,cR);
canonicalAndMinority_COs            = add_typeII_events(n,totalPrecursor,canonical_COs,T2prob);

canonical_COs                       = canonical_COs - precursor_positions;
canonical_COs                       = 1 - canonical_COs;
canonical_COs(canonical_COs==1)     = 0;
canonical_COs(canonical_COs==0)     = NaN;
canonical_COs                       = sort(canonical_COs,2);
object_length(1:n,1)                = 1;
simCOs_canonical                    = horzcat(object_length,canonical_COs);
simCOs_canonical                    = simCOs_canonical.*chrXVlength;

simCOs_canonical_CoC30              = getCoC(simCOs_canonical,30);
simCOs_canonical_spacing            = getCOSpacing(simCOs_canonical);


%get canonical and majority crossovoers
canonicalAndMinority_COs                            = canonicalAndMinority_COs - precursor_positions;
canonicalAndMinority_COs                            = 1 - canonicalAndMinority_COs;
canonicalAndMinority_COs(canonicalAndMinority_COs==1)   = 0;
canonicalAndMinority_COs(canonicalAndMinority_COs==0)   = NaN;
canonicalAndMinority_COs                            = sort(canonicalAndMinority_COs,2);
object_length(1:n,1)                                = 1;
simCOs_canonicalAndMinority                         = horzcat(object_length,canonicalAndMinority_COs);
simCOs_canonicalAndMinority                         = simCOs_canonicalAndMinority.*chrXVlength;

simCOs_canonicalAndMajority_CoC30                   = getCoC(simCOs_canonicalAndMinority,30);
simCOs_canonicalAndMajority_spacing                 = getCOSpacing(simCOs_canonicalAndMinority);


% %remove canonical crossovers to get minority crossovers
simCOs_minority(1:n,1:length(simCOs_canonicalAndMinority(1,:)))  = nan;
total_canonical                                              = countTotalCOs(simCOs_canonical);

simCOs_minority(:,1) = chrXVlength;

for i = 1:n
    simCOs_minority(i,2:end-total_classIs(i)) = setdiff(simCOs_canonicalAndMinority(i,2:end),simCOs_canonical(i,2:end));
end

simCOs_minority_CoC30                   = getCoC(simCOs_minority,30);
simCOs_minority_spacing                 = getCOSpacing(simCOs_minority);
