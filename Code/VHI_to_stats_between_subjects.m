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
function VHI_to_stats_between_subjects( subj, folder_path)

s = length(subj);

% load([folder_path 'subjects' num2str(s) '.mat'])
load([folder_path 'BetweenSubjects_TDPT' num2str(s) '.mat']);
load([folder_path 'BetweenSubjects_TDPTmean' num2str(s) '.mat']);
load([folder_path 'BetweenSubjects_fitResults' num2str(s) '_createFit.mat']);
load([folder_path 'BetweenSubjects_Illusion' num2str(s) '.mat']);

% TDPT0
% BetweenSubjects_TDPT0(:,:) = numb_avams_m(4,:,:);
% BetweenSubjects_TDPT0 = BetweenSubjects_TDPT0';
% BetweenSubjects_TDPT0 = 100*BetweenSubjects_TDPT0/8;
% xlswrite([folder_path 'BetweenSubjects_TDPT0' num2str(s) '.xlsx'], BetweenSubjects_TDPT0, ['A2:B' int2str(s+1)]);
% entete = {'TDPT0 Pre', 'TDPT0 Post', 'Condition'};
% xlswrite([folder_path 'BetweenSubjects_TDPT0' num2str(s) '.xlsx'], entete, 'A1:C1');
% BetweenSubjects_TDPT0delta = BetweenSubjects_TDPT0(:,2) - BetweenSubjects_TDPT0(:,1);
% xlswrite([folder_path 'BetweenSubjects_TDPT0delta' num2str(s) '.xlsx'], BetweenSubjects_TDPT0delta, ['A2:A' int2str(s+1)]);
% entete = {'TDPT0delta', 'Condition'};
% xlswrite([folder_path 'BetweenSubjects_TDPT0delta' num2str(s) '.xlsx'], entete, 'A1:B1');

% Correct Answers
% BetweenSubjects_CA(:,:) = sum(ca_m,1);
% BetweenSubjects_CA = BetweenSubjects_CA';
% BetweenSubjects_CA = 100*BetweenSubjects_CA/48;
% xlswrite([folder_path 'BetweenSubjects_CA' num2str(s) '.xlsx'], BetweenSubjects_CA, ['A2:B' int2str(s+1)]);
% entete = {'CA Pre', 'CA Post', 'Condition'};
% xlswrite([folder_path 'BetweenSubjects_CA' num2str(s) '.xlsx'], entete, 'A1:C1');
% BetweenSubjects_CAdelta = BetweenSubjects_CA(:,2) - BetweenSubjects_CA(:,1);
% xlswrite([folder_path 'BetweenSubjects_CAdelta' num2str(s) '.xlsx'], BetweenSubjects_CAdelta, ['A2:A' int2str(s+1)]);
% entete = {'CAdelta', 'Condition'};
% xlswrite([folder_path 'BetweenSubjects_CAdelta' num2str(s) '.xlsx'], entete, 'A1:B1');

% Forearm Answers
BetweenSubjects_forearm(:,:) = sum(numb_avams_m,1);
BetweenSubjects_forearm = BetweenSubjects_forearm';
BetweenSubjects_forearm = 100*BetweenSubjects_forearm/48;
xlswrite([folder_path 'BetweenSubjects_ForeAns' num2str(s) '.xlsx'], BetweenSubjects_forearm, ['A2:B' int2str(s+1)]);
entete = {'ForeAns Pre', 'ForeAns Post', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_ForeAns' num2str(s) '.xlsx'], entete, 'A1:C1');
BetweenSubjects_forearmdelta = BetweenSubjects_forearm(:,2) - BetweenSubjects_forearm(:,1);
xlswrite([folder_path 'BetweenSubjects_ForeAnsdelta' num2str(s) '.xlsx'], BetweenSubjects_forearmdelta, ['A2:A' int2str(s+1)]);
entete = {'ForeAnsdelta', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_ForeAnsdelta' num2str(s) '.xlsx'], entete, 'A1:B1');

% Point of Subjective Equality
BetweenSubjects_PSE(:,:) =  coefficients(:,1,:);
BetweenSubjects_PSE = BetweenSubjects_PSE';
xlswrite([folder_path 'BetweenSubjects_PSE' num2str(s) '.xlsx'], BetweenSubjects_PSE, ['A2:B' int2str(s+1)]);
entete = {'PSE Pre', 'PSE Post', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_PSE' num2str(s) '.xlsx'], entete, 'A1:C1');
BetweenSubjects_PSEdelta = BetweenSubjects_PSE(:,2) - BetweenSubjects_PSE(:,1);
xlswrite([folder_path 'BetweenSubjects_PSEdelta' num2str(s) '.xlsx'],BetweenSubjects_PSEdelta, ['A2:A' int2str(s+1)]);
entete = {'PSEdelta', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_PSEdelta' num2str(s) '.xlsx'], entete, 'A1:B1');

% Esteem Accuracy
BetweenSubjects_EA(:,:) =  coefficients(:,2,:);
BetweenSubjects_EA = BetweenSubjects_EA';
xlswrite([folder_path 'BetweenSubjects_EA' num2str(s) '.xlsx'], BetweenSubjects_EA, ['A2:B' int2str(s+1)]);
entete = {'EA Pre', 'EA Post', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_EA' num2str(s) '.xlsx'], entete, 'A1:C1');
BetweenSubjects_EAdelta = BetweenSubjects_EA(:,2) - BetweenSubjects_EA(:,1);
xlswrite([folder_path 'BetweenSubjects_EAdelta' num2str(s) '.xlsx'],BetweenSubjects_EAdelta, ['A2:A' int2str(s+1)]);
entete = {'EAdelta', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_EAdelta' num2str(s) '.xlsx'], entete, 'A1:B1');

% Illusion
illuBetween = illuBetween';
xlswrite([folder_path 'BetweenSubjects_Illusion' num2str(s) '.xlsx'],illuBetween, ['A2:D' int2str(s+1)]);
entete = {'RHI Index', 'Vividness', 'Prevalence', 'PD', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_Illusion' num2str(s) '.xlsx'], entete, 'A1:E1');

condition = cell(length(subj),1);
i_20A = 1; i_20S = 1; i_40S = 1;
for s = 1:length(subj)
    cd(subj{s});
    D = dir('*.xlsx');
    disp(subj{s});
    condition{s} = D(2).name([8 9 12]);
    switch condition{s}
        case '20A'
            condition20A(:,i_20A) = coefficients(:,1,s);
            condition20Abis(:,i_20A) = coefficients(:,2,s);
            condition20Ater(i_20A,:) = BetweenSubjects_forearm(s,:);
            i_20A = i_20A + 1;
        case '20S'
            condition20S(:,i_20S) = coefficients(:,1,s);
            condition20Sbis(:,i_20S) = coefficients(:,2,s);
            condition20Ster(i_20S,:) = BetweenSubjects_forearm(s,:);
            i_20S = i_20S + 1;
        case '40S'
            condition40S(:,i_40S) = coefficients(:,1,s);
            condition40Sbis(:,i_40S) = coefficients(:,2,s);
            condition40Ster(i_40S,:) = BetweenSubjects_forearm(s,:);
            i_40S = i_40S + 1;
    end
    clear D
    cd ..
end

% PSE for t-test
condition20A = condition20A';
entete = {'PSE 20A Pre','PSE 20A Post'};
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], condition20A, ['A2:B' int2str(length(condition20A)+1)]);
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], entete, 'A1:B1');
condition20S = condition20S';
entete = {'PSE 20S Pre','PSE 20S Post'};
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], condition20S, ['C2:D' int2str(length(condition20S)+1)]);
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], entete, 'C1:D1');
condition40S = condition40S';
entete = {'PSE 40S Pre','PSE 40S Post'};
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], condition40S, ['E2:F' int2str(length(condition40S)+1)]);
xlswrite([folder_path 'BetweenSubjects_PSE_ttest' num2str(s) '.xlsx'], entete, 'E1:F1');

% EA for t-test
condition20Abis = condition20Abis';
entete = {'EA 20A Pre','EA 20A Post'};
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], condition20Abis, ['A2:B' int2str(length(condition20Abis)+1)]);
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], entete, 'A1:B1');
condition20Sbis = condition20Sbis';
entete = {'EA 20S Pre','EA 20S Post'};
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], condition20Sbis, ['C2:D' int2str(length(condition20Sbis)+1)]);
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], entete, 'C1:D1');
condition40Sbis = condition40Sbis';
entete = {'EA 40S Pre','EA 40S Post'};
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], condition40Sbis, ['E2:F' int2str(length(condition40Sbis)+1)]);
xlswrite([folder_path 'BetweenSubjects_EA_ttest' num2str(s) '.xlsx'], entete, 'E1:F1');

% Forearm for t-test
entete = {'ForeAns 20A Pre','ForeAns 20A Post'};
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], condition20Ater, ['A2:B' int2str(length(condition20Ater)+1)]);
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], entete, 'A1:B1');
entete = {'ForeAns 20S Pre','ForeAns 20S Post'};
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], condition20Ster, ['C2:D' int2str(length(condition20Ster)+1)]);
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], entete, 'C1:D1');
entete = {'ForeAns 40S Pre','ForeAns 40S Post'};
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], condition40Ster, ['E2:F' int2str(length(condition40Ster)+1)]);
xlswrite([folder_path 'BetweenSubjects_ForeAns_ttest' num2str(s) '.xlsx'], entete, 'E1:F1');


% Proprioceptive Drift
pds = xlsread([folder_path 'illusionResults' num2str(s) '.xlsx'],['M2:R' int2str(s+1)]);
condition = cell(length(subj),1);
for k = 1:length(subj)
    cd(subj{k});
    D = dir('*.xlsx');
    disp(subj{k});
    condition{k} = D(2).name([8 9 12]); % Number of the trial
    switch condition{k}
        case '20S'
            pd(k,1:2) = pds(k,3:4);
        case '40S'
            pd(k,1:2) = pds(k,5:6);
        case '20A'
            pd(k,1:2) = pds(k,1:2);
    end
    cd ..
end
entete = {'PD Pre', 'PD Post', 'Condition'};
xlswrite([folder_path 'BetweenSubjects_PD' num2str(s) '.xlsx'], pd, ['A2:B' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_PD' num2str(s) '.xlsx'], condition, ['C2:C' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_PD' num2str(s) '.xlsx'], entete, 'A1:C1');

% PD for t-test
i_20A = 1; i_20S = 1; i_40S = 1;
A = readcell([folder_path 'BetweenSubjects_PD' num2str(s) '.xlsx'],'Range', ['A2:C'  num2str(s+1)]);
for i = 1:s
    switch A{i,3}
        case '20A'
            B(i_20A,1) = A{i,1};
            B(i_20A,2) = A{i,2};
            i_20A = i_20A + 1;
        case '20S'
            B(i_20S,3) = A{i,1};
            B(i_20S,4) = A{i,2};
            i_20S = i_20S + 1;
        case '40S'
            B(i_40S,5) = A{i,1};
            B(i_40S,6) = A{i,2};
            i_40S = i_40S + 1;
        
    end
end
entete = {'PD 20A Pre', 'PD 20A Post', 'PD 20S Pre', 'PD 20S Post', 'PD 40S Pre', 'PD 40S Post'};
xlswrite([folder_path 'BetweenSubjects_PD_ttest' num2str(s) '.xlsx'], B, ['A2:F' int2str(length(B))]);
xlswrite([folder_path 'BetweenSubjects_PD_ttest' num2str(s) '.xlsx'], entete, 'A1:F1');

% Write condition column
% xlswrite([folder_path 'BetweenSubjects_TDPT0delta' num2str(s) '.xlsx'],condition, ['B2:B' int2str(s+1)]);
% xlswrite([folder_path 'BetweenSubjects_TDPT0' num2str(s) '.xlsx'],condition, ['C2:C' int2str(s+1)]);
% xlswrite([folder_path 'BetweenSubjects_CAdelta' num2str(s) '.xlsx'],condition, ['B2:B' int2str(s+1)]);
% xlswrite([folder_path 'BetweenSubjects_CA' num2str(s) '.xlsx'],condition, ['C2:C' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_PSEdelta' num2str(s) '.xlsx'],condition, ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_PSE' num2str(s) '.xlsx'],condition, ['C2:C' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_EAdelta' num2str(s) '.xlsx'],condition, ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_EA' num2str(s) '.xlsx'],condition, ['C2:C' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_Illusion' num2str(s) '.xlsx'],condition, ['E2:E' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_ForeAnsdelta' num2str(s) '.xlsx'],condition, ['B2:B' int2str(s+1)]);
xlswrite([folder_path 'BetweenSubjects_ForeAns' num2str(s) '.xlsx'],condition, ['C2:C' int2str(s+1)]);

end