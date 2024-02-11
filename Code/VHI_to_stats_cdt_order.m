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
function VHI_to_stats_cdt_order( subj, folder_path )

s = length(subj);

load([folder_path 'subjects' num2str(s) '.mat'])
load([folder_path 'CdtOrder_TDPTcount' num2str(s) '.mat']);
load([folder_path 'CdtOrder_TDPTmean' num2str(s) '.mat'])
load([folder_path 'CdtOrder_fitResults' num2str(s) '_createFit.mat'])
load([folder_path 'CdtOrder_Illusion' num2str(s) '.mat'])

% Point of Subject Equality
CdtOrder_PSE(:,:) =  coefficients(:,1,:);
CdtOrder_PSE = CdtOrder_PSE';
xlswrite([folder_path 'CdtOrder_PSE' num2str(s) '.xlsx'],CdtOrder_PSE, ['A2:D' int2str(s+1)]);
entete = {'PSE 1st - Pre', 'PSE 2nd', 'PSE 3rd', 'PSE 4th'};
xlswrite([folder_path 'CdtOrder_PSE' num2str(s) '.xlsx'], entete, 'A1:D1');

% Esteem Accuracy
CdtOrder_EA(:,:) =  coefficients(:,2,:);
CdtOrder_EA = CdtOrder_EA';
xlswrite([folder_path 'CdtOrder_EA' num2str(s) '.xlsx'],CdtOrder_EA, ['A2:D' int2str(s+1)]);
entete = {'EA 1st - Pre', 'EA 2nd', 'EA 3rd', 'EA 4th'};
xlswrite([folder_path 'CdtOrder_EA' num2str(s) '.xlsx'], entete, 'A1:D1');

% TDPT0
CdtOrder_tdpt0(:,:) = numb_avams_m(4,:,:);
CdtOrder_tdpt0 = CdtOrder_tdpt0';
CdtOrder_tdpt0 = 100*CdtOrder_tdpt0/8;
xlswrite([folder_path 'CdtOrder_TDPT0_' num2str(s) '.xlsx'],CdtOrder_tdpt0, ['A2:D' int2str(s+1)]);
entete = {'TDPT0 1st - Pre', 'TDPT0 2nd', 'TDPT0 3rd', 'TDPT0 4th'};
xlswrite([folder_path 'CdtOrder_TDPT0_' num2str(s) '.xlsx'], entete, 'A1:D1');

% Correct Answers
CdtOrder_cor_ans(:,:) = sum(number_c_m,1);
CdtOrder_cor_ans = CdtOrder_cor_ans';
CdtOrder_cor_ans = 100*CdtOrder_cor_ans/48;
xlswrite([folder_path 'CdtOrder_CA' num2str(s) '.xlsx'],CdtOrder_cor_ans, ['A2:D' int2str(s+1)]);
entete = {'CA 1st - Pre', 'CA 2nd', 'CA 3rd', 'CA 4th'};
xlswrite([folder_path 'CdtOrder_CA' num2str(s) '.xlsx'], entete, 'A1:D1');

% Forearm Answers
CdtOrder_forearm(:,:) = sum(numb_avams_m,1);
CdtOrder_forearm = CdtOrder_forearm';
CdtOrder_forearm = 100*CdtOrder_forearm/48;
xlswrite([folder_path 'CdtOrder_ForeAns' num2str(s) '.xlsx'], CdtOrder_forearm, ['A2:D' int2str(s+1)]);
entete = {'ForeAns 1st - Pre', 'ForeAns 2nd', 'ForeAns 3rd', 'ForeAns 4th'};
xlswrite([folder_path 'CdtOrder_ForeAns' num2str(s) '.xlsx'], entete, 'A1:D1');

% Illusion
entete = {'RHI Index 1st','RHI Index 2nd','RHI Index 3rd'};
RHI_index(:,:) = illuOrder(1,:,:);
RHI_index = RHI_index';
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], RHI_index, ['A2:C' int2str(s+1)]);
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], entete, 'A1:C1');
entete = {'Vividness 1st','Vividness 2nd','Vividness 3rd'};
vividness(:,:) = illuOrder(2,:,:);
vividness = vividness';
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], vividness, ['D2:F' int2str(s+1)]);
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], entete, 'D1:F1');
entete = {'Prevalence 1st','Prevalence 2nd','Prevalence 3rd'};
prevalence(:,:) = illuOrder(3,:,:);
prevalence = prevalence';
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], prevalence, ['G2:I' int2str(s+1)]);
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], entete, 'G1:I1');

pds = xlsread([folder_path 'illusionResults' num2str(s) '.xlsx'],['M2:R' int2str(s+1)]);
count = zeros(3);
for k = 1:length(subj)
    cd(subj{k});
    D = dir('*.xlsx');
    disp(subj{k});
    for t = 1:3 % t for Trial
        condition = D(t+1).name(8:end-5); % Number of the trial
        switch condition
            case '20cmS'
                pd(k,(t-1)*2+1:(t-1)*2+2) = pds(k,3:4);
                count(1,t) = count(1,t) + 1;
            case '40cmS'
                pd(k,(t-1)*2+1:(t-1)*2+2) = pds(k,5:6);
                count(2,t) = count(2,t) + 1;
            case '20cmA'
                pd(k,(t-1)*2+1:(t-1)*2+2) = pds(k,1:2);
                count(3,t) = count(3,t) + 1;
        end
    end
    cd ..
end
entete = {'PD Pre 1st', 'PD Post 1st', 'PD Pre 2nd', 'PD Post 2nd',  'PD Pre 3rd', 'PD Post 3rd'};
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], pd, ['J2:O' int2str(s+1)]);
xlswrite([folder_path 'CdtOrder_Illusion' num2str(s) '.xlsx'], entete, 'J1:O1');

end

