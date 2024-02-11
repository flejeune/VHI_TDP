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
function i_fig = VHI_first_plots(s, folder_path, i_fig)
%VHI_first_plots Plot first data
%   Plots TDPT & Questionnaires data

load([folder_path 'subjects' num2str(s) '.mat'])
load([folder_path 'TDPTcounts' num2str(s) '.mat'])
load([folder_path 'TDPTmeans' num2str(s) '.mat'])
load([folder_path 'TDPTstderrs' num2str(s) '.mat'])
load([folder_path 'illusionResults' num2str(s) '.mat'])

%---------------------------------PLOTS-----------------------------------%
                    %------------------------------%
                    % TDPT - Intra-subject results %
                    %------------------------------%
                    
for s = 1:size(risp_1_m,3)
    % PERCENTAGE OF "FOREARM" RESPONSE - numb_avams_m (SPOT OF STIMULUS)
    if s<13
        f1 = figure(i_fig);
        subplot(3,4,s)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        if s<6
            title(['Subj - ' subj{s}(3:end)])
        else
            title(['Subj - ' subj{s}(4:end)])
        end
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<25
        f2 = figure(i_fig + 1);
        subplot(3,4,s-12)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj - ' subj{s}(4:end)])
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<37
        f3 = figure(i_fig + 2);
        subplot(3,4,s-24)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj - ' subj{s}(4:end)])
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<49
        f4 = figure(i_fig + 3);
        subplot(3,4,s-36)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj - ' subj{s}(4:end)])
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif s<61
        f5 = figure(i_fig + 4);
        subplot(3,4,s-48)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj - ' subj{s}(4:end)])
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    else
        f6 = figure(i_fig + 5);
        subplot(3,4,s-60)
        plot(100*numb_avams_m(:,1,s)/8,'-ok'); hold on;
        plot(100*numb_avams_m(:,2,s)/8,'-*b');
        plot(100*numb_avams_m(:,3,s)/8,'-+g');
        plot(100*numb_avams_m(:,4,s)/8,'-xr');
        set(gca,'XTick', [1:7], 'XTickLabel', [-15:5:15]);
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        title(['Subj - ' subj{s}(4:end)])
        if s==size(risp_1_m,3)
            leg_due = legend('Pre','20cmS','40cmS','20cmA');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    end
end

i_fig = i_fig + 6;

savefig(f1, [folder_path '1-individualTDPT1.fig']);
savefig(f2, [folder_path '2-individualTDPT2.fig']);
savefig(f3, [folder_path '3-individualTDPT3.fig']);
savefig(f4, [folder_path '4-individualTDPT4.fig']);
savefig(f5, [folder_path '5-individualTDPT5.fig']);
savefig(f6, [folder_path '6-individualTDPT6.fig']);

                    %-------------------------------%
                    % TDPT - Group analysis results %
                    %-------------------------------%
                    
colori = 'kbgr';
x = -15:5:15;
colori = [[75 75 75]/255; 0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];
% colori = [0 0 0; [0 191 255]/255; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];
% colori = [0 0 0; [0 191 255]/255; [124 252 0]/255; 0.8500 0.3250 0.0980];

% Sigmoïd
f = figure(i_fig); i_fig = i_fig + 1;
% plot(x,mean_numb_avams(:,1),'-ok');
% hold on;
% plot(x,mean_numb_avams(:,2),'-*b');
% plot(x,mean_numb_avams(:,3),'-+g');
% plot(x,mean_numb_avams(:,4),'-xr');
ylim([0 100])
% plot(x,mean_numb_avams(:,1),'-ok', 'Color', [0 0 0],'LineWidth', 1, 'MarkerSize', 8);
% hold on;
% plot(x,mean_numb_avams(:,2),'-db', 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'MarkerSize', 8);
% plot(x,mean_numb_avams(:,3),'-sg', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1, 'MarkerSize', 8);
% plot(x,mean_numb_avams(:,4),'-^r', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1, 'MarkerSize', 8);
plot(x,mean_numb_avams(:,1),'-k', 'Color', colori(1,:),'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(x,mean_numb_avams(:,4),'-r', 'Color', colori(4,:), 'LineWidth', 1, 'MarkerSize', 8);
plot(x,mean_numb_avams(:,2),'-b', 'Color', colori(2,:), 'LineWidth', 1, 'MarkerSize', 8);
plot(x,mean_numb_avams(:,3),'-g', 'Color', colori(3,:), 'LineWidth', 1, 'MarkerSize', 8);
plot(linspace(-15,15,7),ones(1,7)*50,'--', 'Color', [155 155 155]/255)
for i=1:7
    for j=[1 4 2 3]
         e = errorbar(x(i), mean_numb_avams(i,j), stderr_numb_avams(i,j));
         e.Color = colori(j,:);
         e.LineWidth = 1;
    end
end
legend('Pre', '20A', '20S','40S', 'Location', 'Best');
xlabel('\DeltaL [mm]', 'FontWeight', 'bold')
ylabel('Perc. of "forearm" Answers [%]', 'FontWeight', 'bold')
title(['Precentage of Forearm Answers by Stimulation Distance Difference and by Condition'])
ylim([0 100])
savefig(f, [folder_path '6-averageTDPT1.fig']);

% Same plot than above but as bar graph
f = figure(i_fig); i_fig = i_fig + 1;
x = [-15 -10 -5 0 5 10 15];
y = [mean_numb_avams(1,1) mean_numb_avams(1,4) mean_numb_avams(1,2) mean_numb_avams(1,3);
    mean_numb_avams(2,1) mean_numb_avams(2,4) mean_numb_avams(2,2) mean_numb_avams(2,3);
    mean_numb_avams(3,1) mean_numb_avams(3,4) mean_numb_avams(3,2) mean_numb_avams(3,3);
    mean_numb_avams(4,1) mean_numb_avams(4,4) mean_numb_avams(4,2) mean_numb_avams(4,3);
    mean_numb_avams(5,1) mean_numb_avams(5,4) mean_numb_avams(5,2) mean_numb_avams(5,3);
    mean_numb_avams(6,1) mean_numb_avams(6,4) mean_numb_avams(6,2) mean_numb_avams(6,3);
    mean_numb_avams(7,1) mean_numb_avams(7,4) mean_numb_avams(7,2) mean_numb_avams(7,3)];
z = [stderr_numb_avams(1,1) stderr_numb_avams(1,4) stderr_numb_avams(1,2) stderr_numb_avams(1,3);
    stderr_numb_avams(2,1) stderr_numb_avams(2,4) stderr_numb_avams(2,2) stderr_numb_avams(2,3);
    stderr_numb_avams(3,1) stderr_numb_avams(3,4) stderr_numb_avams(3,2) stderr_numb_avams(3,3);
    stderr_numb_avams(4,1) stderr_numb_avams(4,4) stderr_numb_avams(4,2) stderr_numb_avams(4,3);
    stderr_numb_avams(5,1) stderr_numb_avams(5,4) stderr_numb_avams(5,2) stderr_numb_avams(5,3);
    stderr_numb_avams(6,1) stderr_numb_avams(6,4) stderr_numb_avams(6,2) stderr_numb_avams(6,3);
    stderr_numb_avams(7,1) stderr_numb_avams(7,4) stderr_numb_avams(7,2) stderr_numb_avams(7,3)];
b = bar(x, y);
b(1).FaceColor = colori(1,:);
b(2).FaceColor = colori(4,:);
b(3).FaceColor = colori(2,:);
b(4).FaceColor = colori(3,:);
hold on
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = 3.7;
for i=1:7
    for j=1:4
        xx = x(i) - (groupwidth/2) + (2*j-1) * (groupwidth / (2*nbars));
        errorbar(xx, y(i,j), z(i,j), 'k');%colori(j));
    end
end
hold off
legend('Pre','20cmA','20cmS','40cmS');
xlabel('Diff. of the stimulation distance between the forearm and forehead [mm]')
ylabel('Perc. of "forearm" responses')
title(['Mean over ' num2str(length(subj)) ' subjects - Forearm answers - SPOT OF STIMULUS'])
ylim([0 100]);
savefig(f, [folder_path '7-averageTDPT2.fig']);

f = figure(i_fig); i_fig = i_fig + 1;
temp(:,:) = 100*sum(numb_avams_m, 1)/48;
avamsMean = mean(temp,2);
avamsStd = std(temp,0,2)/sqrt(s);
X = categorical({'Pre';'20cmS';'40cmS';'20cmA'});
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
savefig(f, [folder_path '7-forearm.fig']);


                    %----------------------------------%
                    % Illusion - Questionnaire results %
                    %----------------------------------%
                    
X = categorical({'20cmA';'20cmS';'40cmS'});
Xbis = categorical({'Pre';'20cmS';'40cmS';'20cmA'});
colors = 'rbg';

% RHI scores
% Vividness scores
% Prevalence scores
% Proprioceptive drift
f = figure(i_fig); i_fig = i_fig + 1;
subplot(1,4,1)
hold on
for k1 = 1:3
    bar(X(k1), illusion.RHI_index_means(k1), colors(k1));
    errorbar(X(k1),  illusion.RHI_index_means(k1),  illusion.RHI_index_std_errs(k1), 'k')
end
hold off
ylabel('RHI Index')
title('RHI Index')
subplot(1,4,2)
hold on
for k1 = 1:3
    bar(X(k1), illusion.vividness_means(k1), colors(k1));
    errorbar(X(k1),  illusion.vividness_means(k1),  illusion.vividness_std_errs(k1), 'k')
end
hold off
ylabel('Vividness')
title('Vividness')
subplot(1,4,3)
hold on
for k1 = 1:3
    bar(X(k1), illusion.prevalence_means(k1), colors(k1));
    errorbar(X(k1),  illusion.prevalence_means(k1),  illusion.prevalence_std_errs(k1), 'k')
end
hold off
ylabel('Prevalence')
title('Prevalence')
subplot(1,4,4)
hold on
for k1 = 1:3
    bar(X(k1), illusion.pd_means(k1), colors(k1));
    errorbar(X(k1),  illusion.pd_means(k1),  illusion.pd_std_errs(k1), 'k')
end
hold off
ylabel('Proprioceptive Drift')
title('Proprioceptive Drift')
savefig(f, [folder_path '8-illusionResults.fig']);

% Percentage of "forearm" answer per condition
% Percentage of correct answer per condition
%{
colori = 'kbgr';
figure(10)
subplot(1,2,1)
hold on
for k1 = 1:4
    bar(Xbis(k1), mean_avam_v(k1), colori(k1));
    errorbar(Xbis(k1),  mean_avam_v(k1),  stderr_avam_v(k1), 'k')
end
hold off
ylabel('Perc. of "forearm" answers')
title('Perc. of "forearm" answers')
subplot(1,2,2)
hold on
for k1 = 1:4
    bar(Xbis(k1), mean_c(k1), colori(k1));
    errorbar(Xbis(k1),  mean_c(k1),  stderr_c(k1), 'k')
end
hold off
ylabel('Perc. of correct answers')
title('Perc. of correct answers')
%}

                    
end

