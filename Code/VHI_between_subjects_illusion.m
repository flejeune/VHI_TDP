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
function i_fig = VHI_between_subjects_illusion( subj, folder_path, i_fig)

%% QUESTIONNAIRES RESULTS (RHI INDEX, VIVIDNESS, PREVALENCE, PROPRIOCEPTIVE DRIFT)

s=length(subj);

quest_filename = [folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'];

illusion.cdt_20cmA.RHI_index = xlsread(quest_filename, ['D3:D' num2str(s+2)]);
illusion.cdt_20cmA.vividness = xlsread(quest_filename, ['E3:E' num2str(s+2)]);
illusion.cdt_20cmA.prevalence = xlsread(quest_filename, ['F3:F' num2str(s+2)]);
illusion.cdt_20cmA.pd = xlsread(quest_filename, ['I3:I' num2str(s+2)]);

illusion.cdt_20cmS.RHI_index = xlsread(quest_filename, ['L3:L' num2str(s+2)]);
illusion.cdt_20cmS.vividness = xlsread(quest_filename, ['M3:M' num2str(s+2)]);
illusion.cdt_20cmS.prevalence = xlsread(quest_filename, ['N3:N' num2str(s+2)]);
illusion.cdt_20cmS.pd = xlsread(quest_filename, ['Q3:Q' num2str(s+2)]);

illusion.cdt_40cmS.RHI_index = xlsread(quest_filename, ['T3:T' num2str(s+2)]);
illusion.cdt_40cmS.vividness = xlsread(quest_filename, ['U3:U' num2str(s+2)]);
illusion.cdt_40cmS.prevalence = xlsread(quest_filename, ['V3:V' num2str(s+2)]);
illusion.cdt_40cmS.pd = xlsread(quest_filename, ['Y3:Y' num2str(s+2)]);
clear s

count = zeros(3);
for s = 1:length(subj)
    cd(subj{s});
    D = dir('*.xlsx');
    disp(subj{s});
    condition = D(2).name(8:end-5);
    switch condition
        case '20cmS'
            illuBetween(1,s) = illusion.cdt_20cmS.RHI_index(s);
            illuBetween(2,s) = illusion.cdt_20cmS.vividness(s);
            illuBetween(3,s) = illusion.cdt_20cmS.prevalence(s);
            illuBetween(4,s) = illusion.cdt_20cmS.pd(s);
            count(1) = count(1) + 1;
        case '40cmS'
            illuBetween(1,s) = illusion.cdt_40cmS.RHI_index(s);
            illuBetween(2,s) = illusion.cdt_40cmS.vividness(s);
            illuBetween(3,s) = illusion.cdt_40cmS.prevalence(s);
            illuBetween(4,s) = illusion.cdt_40cmS.pd(s);
            count(2) = count(2) + 1;
        case '20cmA'
            illuBetween(1,s) = illusion.cdt_20cmA.RHI_index(s);
            illuBetween(2,s) = illusion.cdt_20cmA.vividness(s);
            illuBetween(3,s) = illusion.cdt_20cmA.prevalence(s);
            illuBetween(4,s) = illusion.cdt_20cmA.pd(s);
            count(3) = count(3) + 1;
    end
    cd ..
end

mean_illuBetween = mean(illuBetween,3);
stderr_illuBetween = std(illuBetween,0,3)/sqrt(s);

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'BetweenSubjects_Illusion' num2str(s) '.mat'],'illuBetween', 'mean_illuBetween','stderr_illuBetween','count');
%---------------------------------PLOT------------------------------------%

                    %----------------------------------%
                    % Illusion - Questionnaire results %
                    %----------------------------------%
%{
X = categorical({'1st';'2nd';'3rd'});
colors = 'rbg';

% RHI scores
% Vividness scores
% Prevalence scores
% Proprioceptive drift
f = figure(20);
subplot(1,4,1)
hold on
for k = 1:3
    bar(X(k), mean_illuBetween(1,k), colors(k));
    errorbar(X(k),  mean_illuBetween(1,k),  stderr_illuBetween(1,k), 'k')
end
hold off
ylabel('RHI Index')
title('RHI Index')
subplot(1,4,2)
hold on
for k = 1:3
    bar(X(k), mean_illuBetween(2,k), colors(k));
    errorbar(X(k),  mean_illuBetween(2,k),  stderr_illuBetween(2,k), 'k')
end
hold off
ylabel('Vividness')
title('Vividness')
subplot(1,4,3)
hold on
for k = 1:3
    bar(X(k), mean_illuBetween(3,k), colors(k));
    errorbar(X(k),  mean_illuBetween(3,k),  stderr_illuBetween(3,k), 'k')
end
hold off
ylabel('Prevalence')
title('Prevalence')
subplot(1,4,4)
hold on
for k = 1:3
    bar(X(k), mean_illuBetween(4,k), colors(k));
    errorbar(X(k),  mean_illuBetween(4,k),  stderr_illuBetween(4,k), 'k')
end
hold off
ylabel('Proprioceptive Drift')
title('Proprioceptive Drift')
savefig(f, [folder_path '20-CdtOrder_IllusionResults.fig']);
%}
                    
end

