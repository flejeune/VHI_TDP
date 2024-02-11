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
function i_fig = VHI_cdt_order_tdpt(subj, folder_path, i_fig)

for s = 1:length(subj) 
    cd(subj{s});
    D = dir('*.xlsx');
    disp(subj{s});
    for t = 1:4 % t for Trial            
        
        Trial = str2num(D(t).name(6)); % Number of the trial

        % Extract results from excel file into vectors (length: 56)
        Stimolo1   = xlsread(D(t).name, 'J2:J57'); % Size of 1st Stimulus
        [~, Area1] = xlsread(D(t).name, 'K2:K57'); % Area of 1st Stimulus
        Stimolo2   = xlsread(D(t).name, 'L2:L57'); % Size of 2nd Stimulus
        [~, Area2] = xlsread(D(t).name, 'M2:M57'); % Area of 2nd Stimulus
        risposta   = xlsread(D(t).name, 'N2:N57'); % Subject's answer

        % Initialization of variables of interest
        %  Number of times the stimulus on the forearm is felt larger, 
        %  for each difference of stimuli length (total : 7), 
        %  regrouped according to the stimulation spot (forearm).
        numb_avams = zeros(7,1);
        %  Number of correct answers for each difference of stimuli length.
        %  Does not consider the case in which stimuli lengths are equal.
        number_c  = zeros(6,1);

        for i = 1:length(Stimolo1) % 56
            % Count number of times the LARGER stimulus is felt on the FOREARM
            if (strcmp(Area1{i}, 'AVAMBRACCIO'))
                Dif = Stimolo1(i) - Stimolo2(i);
                if (risposta(i) == 1)
                    switch Dif
                        case -15; numb_avams(1) = numb_avams(1)+1;
                        case -10; numb_avams(2) = numb_avams(2)+1;
                        case -5;  numb_avams(3) = numb_avams(3)+1;
                        case 0;   numb_avams(4) = numb_avams(4)+1;
                        case 5
                            numb_avams(5) = numb_avams(5)+1;
                            number_c(4) = number_c(4) + 1;
                        case 10
                            numb_avams(6) = numb_avams(6)+1;
                            number_c(5) = number_c(5) + 1;
                        case 15
                            numb_avams(7) = numb_avams(7)+1;
                            number_c(6) = number_c(6) + 1;
                    end
                else
                    switch Dif
                        case -15; number_c(1) = number_c(1) + 1;
                        case -10; number_c(2) = number_c(2) + 1;
                        case -5;  number_c(3) = number_c(3) + 1;
                    end
                end
            else
                Dif = Stimolo2(i) - Stimolo1(i);
                if (risposta(i) == 2) % If answer is two
                    switch Dif
                        case -15; numb_avams(1) = numb_avams(1)+1;
                        case -10; numb_avams(2) = numb_avams(2)+1;
                        case -5;  numb_avams(3) = numb_avams(3)+1;
                        case 0;   numb_avams(4) = numb_avams(4)+1;
                        case 5
                            numb_avams(5) = numb_avams(5)+1;
                            number_c(4) = number_c(4) + 1;
                        case 10
                            numb_avams(6) = numb_avams(6)+1;
                            number_c(5) = number_c(5) + 1;
                        case 15
                            numb_avams(7) = numb_avams(7)+1;
                            number_c(6) = number_c(6) + 1;
                    end
                else
                    switch Dif
                        case -15; number_c(1) = number_c(1) + 1;
                        case -10; number_c(2) = number_c(2) + 1;
                        case -5;  number_c(3) = number_c(3) + 1;
                    end
                end
            end              
        end

        % Save current variable (current subject, current condition) in the general matrix
        numb_avams_m(:,t,s) = numb_avams;
        number_c_m(:,t,s) = number_c;
        clear numb_avams number_c
    end    
    cd ..
end

% Convert in percentage
numb_avams_p = 100*numb_avams_m/8;

% Compute mean
mean_numb_avams = mean(numb_avams_p,3);

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'CdtOrder_TDPTcount' num2str(s) '.mat'],'numb_avams_m', 'number_c_m');

save([folder_path 'CdtOrder_TDPTmean' num2str(s) '.mat'],'mean_numb_avams');

%---------------------------------PLOTS-----------------------------------%
                    %------------------------------%
                    % TDPT - Intra-subject results %
                    %------------------------------%
                    
x = [-15 -10 -5 0 5 10 15];

for s = 1:length(subj)
    % PERCENTAGE OF "FOREARM" RESPONSE - numb_avams_m (SPOT OF STIMULUS)
    if s<13
        f1 = figure(i_fig);
        subplot(3,4,s)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        if s<6
            title(['Subj ' num2str(s) ' - ' subj{s}(3:end)])
        else
            title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        end
        if s==12
            leg = legend('1st - Pre','2nd','3rd','4th');
            leg.Position(1) = 0.87;
            leg.Position(2) = 0.46;
        end
    elseif s<25
        f2 = figure(i_fig + 1);
        subplot(3,4,s-12)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        if s==length(subj)
            leg_due = legend('1st - Pre','2nd','3rd','4th');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<37
        f3 = figure(i_fig + 2);
        subplot(3,4,s-24)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        if s==length(subj)
            leg_due = legend('1st - Pre','2nd','3rd','4th');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<49
        f4 = figure(i_fig + 3);
        subplot(3,4,s-36)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        if s==length(subj)
            leg_due = legend('1st - Pre','2nd','3rd','4th');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<61
        f5 = figure(i_fig + 4);
        subplot(3,4,s-48)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        if s==length(subj)
            leg_due = legend('1st - Pre','2nd','3rd','4th');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    else
        f6 = figure(i_fig + 5);
        subplot(3,4,s-60)
        plot(x,numb_avams_p(:,1,s),'-ok'); hold on;
        plot(x,numb_avams_p(:,2,s),'-*b');
        plot(x,numb_avams_p(:,3,s),'-+g');
        plot(x,numb_avams_p(:,4,s),'-xr');
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj ' num2str(s) ' - ' subj{s}(4:end)])
        if s==length(subj)
            leg_due = legend('1st - Pre','2nd','3rd','4th');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    end
end
savefig(f1, [folder_path '17-CdtOrder_individualTDPT1.fig']);
savefig(f2, [folder_path '18-CdtOrder_individualTDPT2.fig']);
savefig(f3, [folder_path '19-CdtOrder_individualTDPT3.fig']);
savefig(f4, [folder_path '20-CdtOrder_individualTDPT4.fig']);
savefig(f5, [folder_path '21-CdtOrder_individualTDPT5.fig']);
savefig(f6, [folder_path '22-CdtOrder_individualTDPT6.fig']);

i_fig = i_fig + 6;

                    %-------------------------------%
                    % TDPT - Group analysis results %
                    %-------------------------------%
                    
% colori = 'kbgr';
colori = [128 128 128; 65 105 225; 127 255 0; 255 53 53]/255;
% Sigmoïd
f = figure(i_fig); i_fig = i_fig + 1;
plot(x,mean_numb_avams(:,1),'-ok');
hold on;
plot(x,mean_numb_avams(:,2),'-*b');
plot(x,mean_numb_avams(:,3),'-+g');
plot(x,mean_numb_avams(:,4),'-xr');
plot(x,ones(1,length(x))*50);
% for i=1:13
%     for j=1:4
%         errorbar(x(i), mean_numb_avams(i,j), stderr_numb_avams(i,j), colori(j));
%     end
% end
legend('1st - Pre','2nd','3rd','4th');
xlabel('Diff. of the stimulation distance between the forearm and forehead [mm]')
ylabel('Perc. of "forearm" responses')
title(['Mean over ' num2str(length(subj)) ' subjects - Forearm answers - SPOT OF STIMULUS'])
savefig(f, [folder_path '22-CdtOrder_averageTDPT.fig']);

% Forearm answers
f = figure(i_fig); i_fig = i_fig + 1;
temp(:,:) = 100*sum(numb_avams_m, 1)/48;
avamsMean = mean(temp,2);
avamsStd = std(temp,0,2)/sqrt(s);
X = categorical({'1st - Pre';'2nd';'3rd';'4th'});
colori = 'kbgr';
hold on
for k = 1:4
    bar(X(k), avamsMean(k,1), colori(k));
    errorbar(X(k),  avamsMean(k,1),  avamsStd(k,1), 'k')
end
hold off
ylim([0 100])
ylabel('% of forearm answers')
title(['Forearm answers over ' num2str(s) ' subjects'])
savefig(f, [folder_path '22-CdtOrder_forearm.fig']);

end