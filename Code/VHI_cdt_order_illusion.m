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
function i_fig = VHI_cdt_order_illusion( subj, folder_path, i_fig)
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
    for t = 2:4 % t for Trial
        condition = D(t).name(8:end-5); % Number of the trial
        switch condition
            case '20cmS'
                illuOrder(1,t-1,s) = illusion.cdt_20cmS.RHI_index(s);
                illuOrder(2,t-1,s) = illusion.cdt_20cmS.vividness(s);
                illuOrder(3,t-1,s) = illusion.cdt_20cmS.prevalence(s);
                illuOrder(4,t-1,s) = illusion.cdt_20cmS.pd(s);
                count(1,t-1) = count(1,t-1) + 1;
            case '40cmS'
                illuOrder(1,t-1,s) = illusion.cdt_40cmS.RHI_index(s);
                illuOrder(2,t-1,s) = illusion.cdt_40cmS.vividness(s);
                illuOrder(3,t-1,s) = illusion.cdt_40cmS.prevalence(s);
                illuOrder(4,t-1,s) = illusion.cdt_40cmS.pd(s);
                count(2,t-1) = count(2,t-1) + 1;
            case '20cmA'
                illuOrder(1,t-1,s) = illusion.cdt_20cmA.RHI_index(s);
                illuOrder(2,t-1,s) = illusion.cdt_20cmA.vividness(s);
                illuOrder(3,t-1,s) = illusion.cdt_20cmA.prevalence(s);
                illuOrder(4,t-1,s) = illusion.cdt_20cmA.pd(s);
                count(3,t-1) = count(3,t-1) + 1;
        end
    end
    cd ..
end

mean_illuOrder = mean(illuOrder,3);
stderr_illuOrder = std(illuOrder,0,3)/sqrt(s);

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'CdtOrder_Illusion' num2str(s) '.mat'],'illuOrder', 'mean_illuOrder','stderr_illuOrder','count');
%---------------------------------PLOT------------------------------------%
                    %----------------------------------%
                    % Illusion - Questionnaire results %
                    %----------------------------------%
                    
X = categorical({'1st';'2nd';'3rd'});
colors = 'rbg';

% RHI scores
% Vividness scores
% Prevalence scores
% Proprioceptive drift
f = figure(i_fig); i_fig = i_fig + 1;
subplot(1,4,1)
hold on
for k = 1:3
    bar(X(k), mean_illuOrder(1,k), colors(k));
    errorbar(X(k),  mean_illuOrder(1,k),  stderr_illuOrder(1,k), 'k')
end
hold off
ylabel('RHI Index')
title('RHI Index')
subplot(1,4,2)
hold on
for k = 1:3
    bar(X(k), mean_illuOrder(2,k), colors(k));
    errorbar(X(k),  mean_illuOrder(2,k),  stderr_illuOrder(2,k), 'k')
end
hold off
ylabel('Vividness')
title('Vividness')
subplot(1,4,3)
hold on
for k = 1:3
    bar(X(k), mean_illuOrder(3,k), colors(k));
    errorbar(X(k),  mean_illuOrder(3,k),  stderr_illuOrder(3,k), 'k')
end
hold off
ylabel('Prevalence')
title('Prevalence')
subplot(1,4,4)
hold on
for k = 1:3
    bar(X(k), mean_illuOrder(4,k), colors(k));
    errorbar(X(k),  mean_illuOrder(4,k),  stderr_illuOrder(4,k), 'k')
end
hold off
ylabel('Proprioceptive Drift')
title('Proprioceptive Drift')
savefig(f, [folder_path '23-CdtOrder_IllusionResults.fig']);

end

