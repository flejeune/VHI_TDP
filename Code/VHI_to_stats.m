%% Authors : 
%   - Marco D'Alonzo, PhD. Senior research associate.
%       marco.dalonzo@unicampus.it
%   - François Le Jeune, PhD. Post-doctoral fellow.
%       francois.le-jeune@hotmail.fr
%
%
% Affiliation of both authors at time of editing : 
%   - NeXT Lab, Università Campus Bio-Medico di Roma (UCBM), Roma, Italy.

%%
function VHI_to_stats( subj, folder_path )

s = length(subj);

load([folder_path 'subjects' num2str(s) '.mat'])
load([folder_path 'TDPTcounts' num2str(s) '.mat']);
load([folder_path 'TDPTmeans' num2str(s) '.mat']);
load([folder_path 'fitResults' num2str(s) '_createFit.mat']);
load([folder_path 'illusionResults' num2str(s) '.mat']);

% TDPT0
TDPT0(:,:) = numb_avams_m(4,:,:); 
TDPT0 = TDPT0';
TDPT0 = 100*TDPT0/8;
xlswrite([folder_path 'TDPT0_' num2str(s) '.xlsx'],TDPT0(:,1), ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'TDPT0_' num2str(s) '.xlsx'],TDPT0(:,4), ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'TDPT0_' num2str(s) '.xlsx'],TDPT0(:,2:3), ['C2:D' int2str(s+1)]);
entete = {'TDPT0 Pre', 'TDPT0 20A', 'TDPT0 20S', 'TDPT0 40S'};
xlswrite([folder_path 'TDPT0_' num2str(s) '.xlsx'], entete, 'A1:D1');

% Correct Answers
cor_ans(:,:) = sum(number_c_m,1);
cor_ans = cor_ans';
cor_ans = 100*cor_ans/48;
xlswrite([folder_path 'NCorrectAnswers' num2str(s) '.xlsx'],cor_ans(:,1), ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'NCorrectAnswers' num2str(s) '.xlsx'],cor_ans(:,4), ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'NCorrectAnswers' num2str(s) '.xlsx'],cor_ans(:,2:3), ['C2:D' int2str(s+1)]);
entete = {'NCorrectAns Pre', 'NCorrectAns 20A', 'NCorrectAns 20S', 'NCorrectAns 40S'};
xlswrite([folder_path 'NCorrectAnswers' num2str(s) '.xlsx'], entete, 'A1:D1');

% Forearm Answers
forearm(:,:) = sum(numb_avams_m,1);
forearm = forearm';
forearm = 100*forearm/48;
xlswrite([folder_path 'ForearmAnswers' num2str(s) '.xlsx'],forearm(:,1), ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'ForearmAnswers' num2str(s) '.xlsx'],forearm(:,4), ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'ForearmAnswers' num2str(s) '.xlsx'],forearm(:,2:3), ['C2:D' int2str(s+1)]);
entete = {'ForeAns Pre', 'ForeAns 20A', 'ForeAns 20S', 'ForeAns 40S'};
xlswrite([folder_path 'ForearmAnswers' num2str(s) '.xlsx'], entete, 'A1:D1');

% Point of Subjective Equality
PSE(:,:) = coefficients(:,1,:);
PSE = PSE';
xlswrite([folder_path 'PSE' num2str(s) '.xlsx'],PSE(:,1), ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'PSE' num2str(s) '.xlsx'],PSE(:,4), ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'PSE' num2str(s) '.xlsx'],PSE(:,2:3), ['C2:D' int2str(s+1)]);
entete = {'PSE Pre', 'PSE 20A', 'PSE 20S', 'PSE 40S'};
xlswrite([folder_path 'PSE' num2str(s) '.xlsx'], entete, 'A1:D1');

% Point of Subjective Equality - Delta (Post-Pre)
PSE_delta(:,1) = PSE(:,4) - PSE(:,1);
PSE_delta(:,2) = PSE(:,2) - PSE(:,1);
PSE_delta(:,3) = PSE(:,3) - PSE(:,1);
xlswrite([folder_path 'PSEdelta' num2str(s) '.xlsx'],PSE_delta, ['A2:C' int2str(s+1)]);
entete = {'PSEdelta 20A', 'PSEdelta 20S', 'PSEdelta 40S'};
xlswrite([folder_path 'PSEdelta' num2str(s) '.xlsx'], entete, 'A1:C1');

% Esteem Accuracy
EA(:,:) = coefficients(:,2,:);
EA = EA';
xlswrite([folder_path 'EA' num2str(s) '.xlsx'],EA(:,1), ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'EA' num2str(s) '.xlsx'],EA(:,4), ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'EA' num2str(s) '.xlsx'],EA(:,2:3), ['C2:D' int2str(s+1)]);
entete = {'EA Pre', 'EA 20A', 'EA 20S', 'EA 40S'};
xlswrite([folder_path 'EA' num2str(s) '.xlsx'], entete, 'A1:D1');

% Illusion
entete = {'RHI 20A', 'RHI 20S', 'RHI 40S'};
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'A1:C1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmA.RHI_index, ['A2:A' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmS.RHI_index, ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_40cmS.RHI_index, ['C2:C' int2str(s+1)]);
entete = {'Vividness 20A', 'Vividness 20S', 'Vividness 40S'};
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'D1:F1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmA.vividness, ['D2:D' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmS.vividness, ['E2:E' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_40cmS.vividness, ['F2:F' int2str(s+1)]);
entete = {'Prevalence 20A', 'Prevalence 20S', 'Prevalence 40S'};
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'G1:I1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmA.prevalence, ['G2:G' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmS.prevalence, ['H2:H' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_40cmS.prevalence, ['I2:I' int2str(s+1)]);
entete = {'PD 20A', 'PD 20S', 'PD 40S'};
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'J1:L1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmA.pd, ['J2:J' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmS.pd, ['K2:K' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_40cmS.pd, ['L2:L' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'J1:L1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmA.pd, ['J2:J' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_20cmS.pd, ['K2:K' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], illusion.cdt_40cmS.pd, ['L2:L' int2str(s+1)]);
entete = {'PD Pre 20A', 'PD Post 20A', 'PD Pre 20S', 'PD Post 20S', 'PD Pre 40S', 'PD Post 40S'};
PD20APre = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['G3:G' int2str(s+2)]);
PD20APost = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['H3:H' int2str(s+2)]);
PD20SPre = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['O3:O' int2str(s+2)]);
PD20SPost = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['P3:P' int2str(s+2)]);
PD40SPre = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['W3:W' int2str(s+2)]);
PD40SPost = xlsread([folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'], ['X3:X' int2str(s+2)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], entete, 'M1:R1');
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD20APre, ['M2:M' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD20APost, ['N2:N' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD20SPre, ['O2:O' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD20SPost, ['P2:P' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD40SPre, ['Q2:Q' int2str(s+1)]);
xlswrite([folder_path 'illusionResults' num2str(s) '.xlsx'], PD40SPost, ['R2:R' int2str(s+1)]);

end

