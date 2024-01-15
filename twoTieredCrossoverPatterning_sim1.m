%% Martin White. Kleckner Lab. Nov 2023

%% Script Description
%
% Script for simulating two-tiered crossover patterning by the Beam-Film
% Model.
% 
%Logic of simulation: first do 'canonical' BF simulation for majority 
%crossovers, and then take the remaining precursors, with their determined
%sensitivites and do a second patterning process (with larger L) to get 
%minority crossovers

n               = 10000;
N               = 13;
E               = 0.75;
B               = 1;
M               = 1;
A               = 1;
L1              = 0.1;
L2              = 0.2;
Bs              = 1;
Be              = 1;
Bd              = 1;
Y               = 1;
T2prob          = 0;
Bsmax           = 1;
cL              = 0.85;
cR              = 0.85;
Smax            = 3.8;
chrXVlength     = 3.2;

%get fixed precursor array
totalPrecursor_p1                           = event_per_object_from_population_mean(n,N,B);
[precursor_positions_p1,totalPrecursor_p1]  = generate_precursor_array(n,totalPrecursor_p1,E,Bs,Be,Bd);
precursor_sensitivity_matrix_p1             = generate_precursor_sensitivities(n,totalPrecursor_p1,A);

maximum_stress_per_object                   = event_per_object_from_population_mean(n,Smax,Bsmax);

%get majority crossovers
majorityCO_Designations         = pattern_event_designations_according_to_beam_film_model(n,totalPrecursor_p1,precursor_positions_p1,precursor_sensitivity_matrix_p1,maximum_stress_per_object,L1,cL,cR);
majorityCO_Designations         = mature_designated_precursors(n,totalPrecursor_p1,precursor_positions_p1,majorityCO_Designations,M);

majority_COs                    = add_typeII_events(n,totalPrecursor_p1,majorityCO_Designations,T2prob);

majority_COs                    = majority_COs- precursor_positions_p1;
majority_COs                    = 1 - majority_COs;
majority_COs(majority_COs==1)   = 0;
majority_COs(majority_COs==0)   = NaN;
majority_COs                    = sort(majority_COs,2);
object_length(1:n,1)            = 1;
simCOs_majority                 = horzcat(object_length,majority_COs);
simCOs_majority                 = simCOs_majority.*chrXVlength;

simCOs_majority_CoC30           = getCoC(simCOs_majority,30);
simCOs_majority_spacing         = getCOSpacing(simCOs_majority);


% take remaining precursors, reset the stress landscape, and do a second
% patterning process using a larger L to get minority crossovers

%Get the positions of the remaining precursors
precursor_positions_p2                              = majorityCO_Designations;
precursor_positions_p2(precursor_positions_p2==1)   = 0;
precursor_positions_p2(precursor_positions_p2==0)   = NaN;
precursor_positions_p2                              = sort(precursor_positions_p2,2);

totalPrecursor_rnd2(1:n,1)      = nan;
for i = 1:n
    totalPrecursor_rnd2(i,1)    = sum(~isnan(precursor_positions_p2(i,:)));
end

%Get the sensitivity values of the remaining precursors
precursor_sensitivity_matrix_p2(1:n,1:max(totalPrecursor_p1)) = nan;
for i = 1:n
    for j = 1:max(totalPrecursor_p1)
        if majorityCO_Designations(i,j) ~= 1
            precursor_sensitivity_matrix_p2(i,j) = precursor_sensitivity_matrix_p1(i,j);
        end
    end
end

precursor_sensitivity_matrix_p2 = sortMinorityCO_precSensitiv_retained(precursor_sensitivity_matrix_p2);

minorityCO_designations_2L    = pattern_event_designations_according_to_beam_film_model(n,totalPrecursor_rnd2,precursor_positions_p2,precursor_sensitivity_matrix_p2,maximum_stress_per_object,L2,cL,cR);

minority_COs_2L                     = minorityCO_designations_2L - precursor_positions_p2;
minority_COs_2L                     = 1 - minority_COs_2L;
minority_COs_2L(minority_COs_2L==1) = 0;
minority_COs_2L(minority_COs_2L==0) = NaN;
minority_COs_2L                     = sort(minority_COs_2L,2);
object_length(1:n,1)                = 1;
simCOs_minority_2L                  = horzcat(object_length,minority_COs_2L);
simCOs_minority_2L                  = simCOs_minority_2L.*chrXVlength;

simCOs_minority_2L_CoC30            = getCoC(simCOs_minority_2L,30);
simCOs_minority_2L_spacing          = getCOSpacing(simCOs_minority_2L);

% combine canonical majority and minority crossovers
simCOs_majorityAndMinority_2L       = horzcat(simCOs_majority(:,2:end),simCOs_minority_2L(:,2:end));
simCOs_majorityAndMinority_2L       = sort(simCOs_majorityAndMinority_2L,2);
simCOs_majorityAndMinority_2L       = horzcat(simCOs_majority(:,1),simCOs_majorityAndMinority_2L);

simCOs_majorityandMinority_2L_CoC30 = getCoC(simCOs_majorityAndMinority_2L,30);
simCOs_majorityandMinorit_2L_spacing  = getCOSpacing(simCOs_majorityAndMinority_2L);


