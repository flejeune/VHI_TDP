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
function i_fig = VHI_analysis(subj, column_name, folder_path, i_fig)
%VHI_analysis - Analyse TDPT & Questionnaire data
%   Compute metrics of interest from TDPT and save them
%   Compute metrics of interest from questionnaire and save them

%   Inputs:
%       - subj:        cell array containing participants' folder name
%       - column_name: cell array containing column name for Excel file 'All_Data.xls'
%       - folder_path: folder path of the folder where data is saved


mAvam = zeros(56,4,length(subj));
% Compute variables of interest
for s = 1:length(subj) % Nb of subjects
    eval(sprintf(['cd ' subj{s}]));
    D = dir('*.xlsx');
    disp(subj{s});
    for t = 1:length(D) % 1 to 6
        if (size(D(t).name,2)>6 && strcmp(D(t).name(1:5),'Trial')) % If correct name of variable in field "name"

            switch D(t).name(8:end-5)
                case 'Pre'  ,cond = 1;
                case '20cmS',cond = 2;
                case '40cmS',cond = 3;
                case '20cmA',cond = 4;
            end

            Trial = str2num(D(t).name(6)); % Number of the trial

            % Extract results from excel file into vectors (length: 56)
            Stimolo1   = xlsread(D(t).name, 'J2:J57'); % Size of 1st Stimulus
            [~, Area1] = xlsread(D(t).name, 'K2:K57'); % Area of 1st Stimulus
            Stimolo2   = xlsread(D(t).name, 'L2:L57'); % Size of 2nd Stimulus
            [~, Area2] = xlsread(D(t).name, 'M2:M57'); % Area of 2nd Stimulus
            risposta   = xlsread(D(t).name, 'N2:N57'); % Subject's answer

            % Initialization of variables of interest
            number_c  = zeros(6,1);  % Number of correct responses for each pair of stimulus. Excludes the case where stimuli are equal. To be divided by 8.
            numb_avam = zeros(7,1);  % Number of responses where the stimulus on the forearm is felt larger for each pair of stimulus. To be divided by eight. It is grouped with respect to the order of the pair of stimulus.
            avam_c = zeros(6,1);     % Number of responses where the stimulus on the forearm is correctly felt longer for each pair of stimulus. It is grouped with respect to the order of the pair of stimulus. Excludes the case where the stimuli are equal. To be divided by 4.
            avam_c_s = zeros(6,1);   % Number of responses in which the response is correct for each pair of stimulus. It is grouped with respect to the reference stimulation site. Excludes the case in which the stimuli are equal. To be divided by 8.
            risp_1 = zeros(7,1);     % Number of responses in which the response was 1 for each pair of stimulus. To be divided by 8.
            numb_avams = zeros(7,1); % Number of responses where the stimulus on the forearm is felt larger for each pair of stimulus. It is grouped with respect to the reference stimulation site. To be divided by 8.
            avam = 0;

            for i = 1:length(Stimolo1) % 56
                Dif = Stimolo1(i) - Stimolo2(i);
                n = ((Dif > 0)&(risposta(i) == 1))|((Dif < 0)&(risposta(i) == 2)); % Correct n=1, wrong n=0
                % Count number of CORRECT answers & number of answers "1"
                switch Dif
                    case -15
                        number_c(1) = number_c(1)+n;
                        if (risposta(i) == 1); risp_1(1) = risp_1(1)+1; end
                    case -10
                        number_c(2) = number_c(2)+n;
                        if (risposta(i) == 1); risp_1(2) = risp_1(2)+1; end
                    case -5
                        number_c(3) = number_c(3)+n;
                        if (risposta(i) == 1); risp_1(3) = risp_1(3)+1; end
                    case 0
                        if (risposta(i) == 1); risp_1(4) = risp_1(4)+1; end
                    case 5
                        number_c(4) = number_c(4)+n;
                        if (risposta(i) == 1); risp_1(5) = risp_1(5)+1; end
                    case 10
                        number_c(5) = number_c(5)+n;
                        if (risposta(i) == 1); risp_1(6) = risp_1(6)+1; end
                    case 15
                        number_c(6) = number_c(6)+n;
                        if (risposta(i) == 1); risp_1(7) = risp_1(7)+1; end
                end

                % Count number of times the LARGER stimulus is felt on the FOREARM
                if (strcmp(Area1{i}, 'AVAMBRACCIO')) % If first stimulus on forearm
                    Dif_a = Dif; % Dif_a = stimAvambraccio - stimTesta
                    if (risposta(i) == 1) % If answer is "one"
                        mAvam(i,t,s) = 1;
                        switch Dif_a
                            case -15
                                numb_avams(1) = numb_avams(1)+1;
                            case -10
                                numb_avams(2) = numb_avams(2)+1;
                            case -5
                                numb_avams(3) = numb_avams(3)+1;
                            case 0
                                numb_avams(4) = numb_avams(4)+1;
                            case 5
                                numb_avams(5) = numb_avams(5)+1;
                            case 10
                                numb_avams(6) = numb_avams(6)+1;
                            case 15
                                numb_avams(7) = numb_avams(7)+1;
                        end
                    end
                else % If first stimulus on head, second on forearm
                    Dif_a = -Dif;
                    if (risposta(i) == 2) % If answer is two
                        mAvam(i,t,s) = 1;
                        switch Dif_a
                            case -15
                                numb_avams(1) = numb_avams(1)+1;
                            case -10
                                numb_avams(2) = numb_avams(2)+1;
                            case -5
                                numb_avams(3) = numb_avams(3)+1;
                            case 0
                                numb_avams(4) = numb_avams(4)+1;
                            case 5
                                numb_avams(5) = numb_avams(5)+1;
                            case 10
                                numb_avams(6) = numb_avams(6)+1;
                            case 15
                                numb_avams(7) = numb_avams(7)+1;
                        end
                    end
                end              
                % Count number of CORRECT answers - SPOT OF STIMULUS
                if (strcmp(Area1{i}, 'AVAMBRACCIO')) % If first stimulus on forearm
                    Dif_a = Dif; % Dif_a = stimAvam - stimTesta
                    switch Dif_a
                        case -15
                            if(risposta(i) == 2); avam_c_s(1) = avam_c_s(1) + 1; end
                        case -10
                            if(risposta(i) == 2); avam_c_s(2) = avam_c_s(2) + 1; end
                        case -5
                            if(risposta(i) == 2); avam_c_s(3) = avam_c_s(3) + 1; end
                        case 5
                            if(risposta(i) == 1); avam_c_s(4) = avam_c_s(4) + 1; end
                        case 10
                            if(risposta(i) == 1); avam_c_s(5) = avam_c_s(5) + 1; end
                        case 15
                            if(risposta(i) == 1); avam_c_s(6) = avam_c_s(6) + 1; end
                    end
                else % If first stimulus on head, second on forearm
                    Dif_a = -Dif; % Dif_a = StimAvam - StimTesta
                    switch Dif_a
                        case -15
                            if(risposta(i) == 1); avam_c_s(1) = avam_c_s(1) + 1; end
                        case -10
                            if(risposta(i) == 1); avam_c_s(2) = avam_c_s(2) + 1; end
                        case -5
                            if(risposta(i) == 1); avam_c_s(3) = avam_c_s(3) + 1; end
                        case 5
                            if(risposta(i) == 2); avam_c_s(4) = avam_c_s(4) + 1; end
                        case 10
                            if(risposta(i) == 2); avam_c_s(5) = avam_c_s(5) + 1; end
                        case 15
                            if(risposta(i) == 2); avam_c_s(6) = avam_c_s(1) + 1; end
                    end
                end          
                % If first stimulus on forearm and answer is "one"
                % OR
                % If first stimulus on head and answer is "two"
                % Count number of times the LARGER stimulus is felt on the FOREARM(numb_avam)
                % Count number of times the LARGER stimulus is CORRECTLY felt on the FOREARM(avam_c)
                if (strcmp(Area1{i}, 'AVAMBRACCIO')&&(risposta(i) == 1))||(strcmp(Area1{i}, 'TESTA')&&(risposta(i) == 2))
                    avam = avam+1;
                    switch Dif
                        case -15
                            numb_avam(1) = numb_avam(1)+1;
                            if (risposta(i) == 2); avam_c(1) = avam_c(1)+1; end
                        case -10
                            numb_avam(2) = numb_avam(2)+1;
                            if (risposta(i) == 2); avam_c(2) = avam_c(2)+1; end
                        case -5
                            numb_avam(3) = numb_avam(3)+1;
                            if (risposta(i) == 2); avam_c(3) = avam_c(3)+1; end
                        case 0
                            numb_avam(4) = numb_avam(4)+1;
                        case 5
                            numb_avam(5) = numb_avam(5)+1;
                            if (risposta(i) == 1); avam_c(4) = avam_c(4)+1; end
                        case 10
                            numb_avam(6) = numb_avam(6)+1;
                            if (risposta(i) == 1); avam_c(5) = avam_c(5)+1; end
                        case 15
                            numb_avam(7) = numb_avam(7)+1;
                            if (risposta(i) == 1); avam_c(6) = avam_c(6)+1; end
                    end
                end

            end

            % Save current variable (current subject, current condition) in the general matrix
            number_c_m(:,cond,s) = number_c;     % size : 6*nbConditions*nbSubjects
            numb_avam_m(:,cond,s) = numb_avam;   % size : 7*nbConditions*nbSubjects
            avam_c_m(:,cond,s) = avam_c;         % size : 6*nbConditions*nbSubjects
            avam_c_s_m(:,cond,s) = avam_c_s;     % size : 6*nbConditions*nbSubjects
            risp_1_m(:,cond,s) = risp_1;         % size : 7*nbConditions*nbSubjects
            numb_avams_m(:,cond,s) = numb_avams;
            avam_v(cond,s) = avam;
            trial_v(cond,s) = Trial;
            clear number_c numb_avam avam_c avam_c_s risp_1 numb_avams avam Trial
        end
    end    
    cd ..
end

% Write all data in an Excel file called "All_Data.xls"
for cnd = 1:size(risp_1_m,2)
    xlswrite([folder_path 'All_Data'],squeeze(100*risp_1_m(:,cnd,:)/8),'Resp1_percentage',[column_name{2} num2str(3+(cnd-1)*(1+size(risp_1_m,1))) ':' column_name{1+size(risp_1_m,3)} num2str(9+(cnd-1)*(1+size(risp_1_m,1)))]);
    xlswrite([folder_path 'All_Data'],squeeze(100*numb_avams_m(:,cnd,:)/8),'Respforearm_percentage',[column_name{2} num2str(3+(cnd-1)*(1+size(risp_1_m,1))) ':' column_name{1+size(risp_1_m,3)} num2str(9+(cnd-1)*(1+size(risp_1_m,1)))]);   
end
titleExcel = {'Pre' '-15' '-10' '-5' '0' '5' '10' '15' '20cmS' '-15' '-10' '-5' '0' '5' '10' '15' '40cmS' '-15' '-10' '-5' '0' '5' '10' '15' '20cmA' '-15' '-10' '-5' '0' '5' '10' '15'};
xlswrite([folder_path 'All_Data'],subj,'Resp1_percentage',         [column_name{2} '1:' column_name{1+size(risp_1_m,3)} '1']);
xlswrite([folder_path 'All_Data'],subj,'Respforearm_percentage',   [column_name{2} '1:' column_name{1+size(risp_1_m,3)} '1']);
xlswrite([folder_path 'All_Data'],subj,'Respforearmtot_percentage',[column_name{2} '1:' column_name{1+size(risp_1_m,3)} '1']);
xlswrite([folder_path 'All_Data'],titleExcel','Resp1_percentage',      ['A2:A' num2str(1+size(titleExcel,2))]);
xlswrite([folder_path 'All_Data'],titleExcel','Respforearm_percentage',['A2:A' num2str(1+size(titleExcel,2))]);
xlswrite([folder_path 'All_Data'],{titleExcel{1} titleExcel{9} titleExcel{17} titleExcel{25}}','Respforearmtot_percentage',['A2:A' num2str(1+size(avam_v,1))]);
xlswrite([folder_path 'All_Data'],avam_v/56,'Respforearmtot_percentage',[column_name{2} '2:' column_name{1+size(avam_v,2)} num2str(1+size(avam_v,1))]);
trial_name = {'Trial1' 'Trial2' 'Trial3' 'Trial4'};
%xlswrite('All_Data',trial_name,'Trial_Order',[column_name{2} '1:' column_name{1+size(trial_name,2)} '1']);
xlswrite([folder_path 'All_Data'],{titleExcel{1} titleExcel{9} titleExcel{17} titleExcel{25}}','Trial_order',['A2:A' num2str(1+size(avam_v,1))]);
xlswrite([folder_path 'All_Data'],subj,'Trial_Order',[column_name{2} '1:' column_name{1+size(risp_1_m,3)} '1']);
xlswrite([folder_path 'All_Data'],trial_v,'Trial_Order',[column_name{2} '2:' column_name{1+size(trial_v,2)} num2str(1+size(trial_v,1))]);

% Compute means
mean_number_c   = 100*mean(number_c_m,3)/8;
mean_numb_avam  = 100*mean(numb_avam_m,3)/8;
mean_numb_avams = 100*mean(numb_avams_m,3)/8;
mean_avam_c     = 100*mean(avam_c_m,3)/4;
mean_avam_c_s   = 100*mean(avam_c_s_m,3)/8;
mean_risp_1     = 100*mean(risp_1_m,3)/8;
mean_avam_v     = 100*mean(avam_v,2)/56;
mean_c          = 100*mean(sum(number_c_m,1),3)/48;

% Compute standard errros
stderr_number_c   = 100*(std(number_c_m,0,3)/sqrt(s))/8;
stderr_numb_avam  = 100*(std(numb_avam_m,0,3)/sqrt(s))/8;
stderr_numb_avams = 100*(std(numb_avams_m,0,3)/sqrt(s))/8;
stderr_avam_c     = 100*(std(avam_c_m,0,3)/sqrt(s))/8;
stderr_avam_c_s   = std(100*avam_c_s_m/8,0,3)/sqrt(s);
stderr_risp_1     = 100*(std(risp_1_m,0,3)/sqrt(s))/8;
stderr_avam_v     = std(100*avam_v/56,0,2)/sqrt(s);
stderr_c          = std(100*sum(number_c_m,1)/48,0,3)/sqrt(s);
% stderr_avam_v     = 100*(std(avam_v,0,2)/sqrt(s))/56;
% stderr_c          = 100*(std(sum(number_c_m,1),0,3)/sqrt(s))/48;

% Compute medians
median_number_c   = 100*median(number_c_m,3)/8;
median_numb_avam  = 100*median(numb_avam_m,3)/8;
median_numb_avams = 100*median(numb_avams_m,3)/8;
median_avam_c     = 100*median(avam_c_m,3)/4;
median_risp_1     = 100*median(risp_1_m,3)/8;

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'subjects' num2str(s) '.mat'],'subj','s');

save([folder_path 'TDPTcounts' num2str(s) '.mat'],'number_c_m','numb_avam_m','avam_c_m','avam_c_s_m', 'risp_1_m',...
    'numb_avams_m', 'avam_v', 'trial_v');

save([folder_path 'TDPTmeans' num2str(s) '.mat'],'mean_number_c','mean_numb_avam','mean_numb_avams',...
    'mean_avam_c','mean_avam_c_s','mean_risp_1','mean_avam_v','mean_c');

save([folder_path 'TDPTstderrs' num2str(s) '.mat'],'stderr_number_c','stderr_numb_avam','stderr_numb_avams',...
    'stderr_avam_c','stderr_avam_c_s','stderr_risp_1','stderr_avam_v','stderr_c');

save([folder_path 'TDPTmAvam' num2str(s) '.mat'], 'mAvam');

%% QUESTIONNAIRES RESULTS (RHI INDEX, VIVIDNESS, PREVALENCE, PROPRIOCEPTIVE DRIFT)

quest_filename = [folder_path 'QuestionnairesAnswers' num2str(s) '.xlsx'];

% 20A
illusion.cdt_20cmA.RHI_index = xlsread(quest_filename, ['D3:D' num2str(s+2)]);
illusion.cdt_20cmA.RHI_index_mean = mean(illusion.cdt_20cmA.RHI_index);
illusion.cdt_20cmA.RHI_index_std_dev = std(illusion.cdt_20cmA.RHI_index);
illusion.cdt_20cmA.RHI_index_std_err = std(illusion.cdt_20cmA.RHI_index)/sqrt(s);

illusion.cdt_20cmA.vividness = xlsread(quest_filename, ['E3:E' num2str(s+2)]);
illusion.cdt_20cmA.vividness_mean = mean(illusion.cdt_20cmA.vividness);
illusion.cdt_20cmA.vividness_std_dev = std(illusion.cdt_20cmA.vividness);
illusion.cdt_20cmA.vividness_std_err = std(illusion.cdt_20cmA.vividness)/sqrt(s);

illusion.cdt_20cmA.prevalence = xlsread(quest_filename, ['F3:F' num2str(s+2)]);
illusion.cdt_20cmA.prevalence_mean = mean(illusion.cdt_20cmA.prevalence);
illusion.cdt_20cmA.prevalence_std_dev = std(illusion.cdt_20cmA.prevalence);
illusion.cdt_20cmA.prevalence_std_err = std(illusion.cdt_20cmA.prevalence)/sqrt(s);

illusion.cdt_20cmA.pd = xlsread(quest_filename, ['I3:I' num2str(s+2)]);
illusion.cdt_20cmA.pd_mean = mean(illusion.cdt_20cmA.pd);
illusion.cdt_20cmA.pd_std_dev = std(illusion.cdt_20cmA.pd);
illusion.cdt_20cmA.pd_std_err = std(illusion.cdt_20cmA.pd)/sqrt(s);

% 20S
illusion.cdt_20cmS.RHI_index = xlsread(quest_filename, ['L3:L' num2str(s+2)]);
illusion.cdt_20cmS.RHI_index_mean = mean(illusion.cdt_20cmS.RHI_index);
illusion.cdt_20cmS.RHI_index_std_dev = std(illusion.cdt_20cmS.RHI_index);
illusion.cdt_20cmS.RHI_index_std_err = std(illusion.cdt_20cmS.RHI_index)/sqrt(s);

illusion.cdt_20cmS.vividness = xlsread(quest_filename, ['M3:M' num2str(s+2)]);
illusion.cdt_20cmS.vividness_mean = mean(illusion.cdt_20cmS.vividness);
illusion.cdt_20cmS.vividness_std_dev = std(illusion.cdt_20cmS.vividness);
illusion.cdt_20cmS.vividness_std_err = std(illusion.cdt_20cmS.vividness)/sqrt(s);

illusion.cdt_20cmS.prevalence = xlsread(quest_filename, ['N3:N' num2str(s+2)]);
illusion.cdt_20cmS.prevalence_mean = mean(illusion.cdt_20cmS.prevalence);
illusion.cdt_20cmS.prevalence_std_dev = std(illusion.cdt_20cmS.prevalence);
illusion.cdt_20cmS.prevalence_std_err = std(illusion.cdt_20cmS.prevalence)/sqrt(s);

illusion.cdt_20cmS.pd = xlsread(quest_filename, ['Q3:Q' num2str(s+2)]);
illusion.cdt_20cmS.pd_mean = mean(illusion.cdt_20cmS.pd);
illusion.cdt_20cmS.pd_std_dev = std(illusion.cdt_20cmS.pd);
illusion.cdt_20cmS.pd_std_err = std(illusion.cdt_20cmS.pd)/sqrt(s);

% 40S
illusion.cdt_40cmS.RHI_index = xlsread(quest_filename, ['T3:T' num2str(s+2)]);
illusion.cdt_40cmS.RHI_index_mean = mean(illusion.cdt_40cmS.RHI_index);
illusion.cdt_40cmS.RHI_index_std_dev = std(illusion.cdt_40cmS.RHI_index);
illusion.cdt_40cmS.RHI_index_std_err = std(illusion.cdt_40cmS.RHI_index)/sqrt(s);

illusion.cdt_40cmS.vividness = xlsread(quest_filename, ['U3:U' num2str(s+2)]);
illusion.cdt_40cmS.vividness_mean = mean(illusion.cdt_40cmS.vividness);
illusion.cdt_40cmS.vividness_std_dev = std(illusion.cdt_40cmS.vividness);
illusion.cdt_40cmS.vividness_std_err = std(illusion.cdt_40cmS.vividness)/sqrt(s);

illusion.cdt_40cmS.prevalence = xlsread(quest_filename, ['V3:V' num2str(s+2)]);
illusion.cdt_40cmS.prevalence_mean = mean(illusion.cdt_40cmS.prevalence);
illusion.cdt_40cmS.prevalence_std_dev = std(illusion.cdt_40cmS.prevalence);
illusion.cdt_40cmS.prevalence_std_err = std(illusion.cdt_40cmS.prevalence)/sqrt(s);

illusion.cdt_40cmS.pd = xlsread(quest_filename, ['Y3:Y' num2str(s+2)]);
illusion.cdt_40cmS.pd_mean = mean(illusion.cdt_40cmS.pd);
illusion.cdt_40cmS.pd_std_dev = std(illusion.cdt_40cmS.pd);
illusion.cdt_40cmS.pd_std_err = std(illusion.cdt_40cmS.pd)/sqrt(s);

% RHI_indices
illusion.RHI_index_means = [illusion.cdt_20cmA.RHI_index_mean;
    illusion.cdt_20cmS.RHI_index_mean;
    illusion.cdt_40cmS.RHI_index_mean];
illusion.RHI_index_std_devs = [illusion.cdt_20cmA.RHI_index_std_dev
    illusion.cdt_20cmS.RHI_index_std_dev
    illusion.cdt_40cmS.RHI_index_std_dev];
illusion.RHI_index_std_errs = [illusion.cdt_20cmA.RHI_index_std_err
    illusion.cdt_20cmS.RHI_index_std_err
    illusion.cdt_40cmS.RHI_index_std_err];

% Vividness
illusion.vividness_means = [illusion.cdt_20cmA.vividness_mean
    illusion.cdt_20cmS.vividness_mean
    illusion.cdt_40cmS.vividness_mean];
illusion.vividness_std_devs = [illusion.cdt_20cmA.vividness_std_dev 
    illusion.cdt_20cmS.vividness_std_dev
    illusion.cdt_40cmS.vividness_std_dev];
illusion.vividness_std_errs = [illusion.cdt_20cmA.vividness_std_err 
    illusion.cdt_20cmS.vividness_std_err
    illusion.cdt_40cmS.vividness_std_err];

% Prevalence
illusion.prevalence_means = [illusion.cdt_20cmA.prevalence_mean
    illusion.cdt_20cmS.prevalence_mean
    illusion.cdt_40cmS.prevalence_mean];
illusion.prevalence_std_devs = [illusion.cdt_20cmA.prevalence_std_dev 
    illusion.cdt_20cmS.prevalence_std_dev
    illusion.cdt_40cmS.prevalence_std_dev];
illusion.prevalence_std_errs = [illusion.cdt_20cmA.prevalence_std_err
    illusion.cdt_20cmS.prevalence_std_err
    illusion.cdt_40cmS.prevalence_std_err];

% Proprioceptive Drift
illusion.pd_means = [illusion.cdt_20cmA.pd_mean
    illusion.cdt_20cmS.pd_mean
    illusion.cdt_40cmS.pd_mean];
illusion.pd_std_devs = [illusion.cdt_20cmA.pd_std_dev 
    illusion.cdt_20cmS.pd_std_dev
    illusion.cdt_40cmS.pd_std_dev];
illusion.pd_std_errs = [illusion.cdt_20cmA.pd_std_err 
    illusion.cdt_20cmS.pd_std_err
    illusion.cdt_40cmS.pd_std_err];

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'illusionResults' num2str(s) '.mat'], 'illusion');

% END OF FUNCTION
end