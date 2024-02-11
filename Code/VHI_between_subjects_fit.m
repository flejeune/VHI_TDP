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
function i_fig = VHI_between_subjects_fit( subj, folder_path, i_fig)

load([folder_path 'BetweenSubjects_TDPT' num2str(length(subj)) '.mat'])
load([folder_path 'BetweenSubjects_TDPTmean' num2str(length(subj)) '.mat'])

%% Fit sigmoïd on every TDPT results 
% i.e for every condition (4) of every valid subject
x = [-15.0 -10.0 -5.0 0.0 5.0 10.0 15.0];
for i=1:length(subj)
    disp(['--------------Computing fits for Subj ' num2str(i) '--------------'])
    for j=1:2
        % Compute fitting
        disp(['------------------------Cdt° ' num2str(j) '------------------------'])
        [fitresult,gof] = createFit(x,100*numb_avams_m(:,j,i)'/8);
        if i<6
            fitResults{i}{j}{1,1} = ['Subj ' subj{i}(1) ' - ' subj{i}(3:end)];
        else
            fitResults{i}{j}{1,1} = ['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)];
        end
        % Store results
        fitResults{i}{j}{2,1} = fitresult;
        fitResults{i}{j}{3,1} = gof;
        coefficients(j,:,i) = coeffvalues(fitresult);
        clear fitresult gof output
    end
end

% Save results
save([folder_path 'BetweenSubjects_fitResults' num2str(s) '_createFit.mat'], 'x', 'fitResults', 'coefficients')

%% Fit sigmoïd on averaged TDPT results
[fitobjectPre,gofPre] = createFit(x,mean_numb_avams(:,1)');
[fitobjectPost,gofPost] = createFit(x,mean_numb_avams(:,2)');

coeffsPre = coeffvalues(fitobjectPre);
coeffsPost = coeffvalues(fitobjectPost);

% Save results
save([folder_path 'BetweenSubjects_fitOfAverage' num2str(s) '.mat'], 'fitobjectPre', 'fitobjectPost',...
    'gofPre', 'gofPost', 'coeffsPre', 'coeffsPost');

%% Plot results

for i=1:length(subj)
    if i<13
        f1 = figure(i_fig);
        subplot(3,4,i)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)   
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        if i<6
            title(['Subj ' subj{i}(1) ' - ' subj{i}(3:end)])
        else
            title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        end
        if i==12
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif i<25
        f2 = figure(i_fig + 1);
        subplot(3,4,i-12)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        if i==length(subj)
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif i<37
        f3 = figure(i_fig + 2);
        subplot(3,4,i-24)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        if i==length(subj)
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif i<49
        f4 = figure(i_fig + 3);
        subplot(3,4,i-36)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        if i==length(subj)
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    elseif i<61
        f5 = figure(i_fig + 4);
        subplot(3,4,i-48)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        if i==length(subj)
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    else
        f6 = figure(i_fig + 5);
        subplot(3,4,i-60)
        % Cdt Pre
        plot(x, 100*numb_avams_m(:,1,i)/8,'ok'); hold on;
        hf = plot(fitResults{i}{1}{2}, 'k'); legend off;
        sigData{i}{1} = [hf.XData; hf.YData];
        text(-19, 95, num2str(fitResults{i}{1}{3}.rsquare), 'Color', 'k')
        % Cdt Post
        plot(x, 100*numb_avams_m(:,2,i)/8,'*b');
        hf = plot(fitResults{i}{2}{2}, 'b'); legend off;
        sigData{i}{2} = [hf.XData; hf.YData];
        text(-19, 85, num2str(fitResults{i}{2}{3}.rsquare), 'Color', 'b')
        plot(linspace(-20,20,7),ones(1,7)*50)
        xlabel({'Diff. of the stimulation distance','between the forearm and forehead [mm]'})
        ylabel('Perc. of "forearm" responses')
        ylim([0 100])
        title(['Subj ' subj{i}(1:2) ' - ' subj{i}(4:end)])
        if i==length(subj)
            leg_due = legend('1st-Pre', 'SigPre', '2nd-Post', 'SigPost');
            leg_due.Position(1) = 0.87;
            leg_due.Position(2) = 0.46;
        end
    end
end
savefig(f1, [folder_path '40-BetweenSubjects_individualFit1.fig']);
savefig(f2, [folder_path '41-BetweenSubjects_individualFit2.fig']);
savefig(f3, [folder_path '42-BetweenSubjects_individualFit3.fig']);
savefig(f4, [folder_path '43-BetweenSubjects_individualFit4.fig']);
savefig(f5, [folder_path '44-BetweenSubjects_individualFit5.fig']);
savefig(f6, [folder_path '45-BetweenSubjects_individualFit6.fig']);

i_fig = i_fig + 6;

% Average of Fits
f = figure(i_fig); i_fig = i_fig + 1;
for j=1:2 % for each cdt°, store YData of all subjects in a table. It's easier to compute mean.
    for i=1:length(subj)
        sigXData{j}(i,:) = sigData{i}{j}(1,:); % this is just to check if XData is the same for every fit
        sigYData{j}(i,:) = sigData{i}{j}(2,:); 
    end
end
for j=1:2 % compute the mean fit curve for each cdt°
    mean_sigYData(j,:) = mean(sigYData{j},1);
end
plot(sigXData{1}(1,:),mean_sigYData(1,:),'k');
hold on
plot(sigXData{1}(1,:),mean_sigYData(2,:),'b');
plot(sigXData{1}(1,:),ones(1,length(sigXData{1}(1,:)))*50)
legend('1st - Pre','2nd - Post')
xlabel('Diff. of the stimulation distance between the forearm and forehead [mm]')
ylabel('Perc. of "forearm" responses')
title({['Average Of Fits (over ' num2str(length(subj)) ' subjects)']; 'Forearm answers - SPOT OF STIMULUS'})
savefig(f, [folder_path '46-BetweenSubjects_average_of_fits.fig']);

% Fit Of Average
f = figure(i_fig); i_fig = i_fig + 1;
plot(x,mean_numb_avams(:,1),'ok');
hold on;
plot(fitobjectPre, 'k')
text(-19, 95, num2str(gofPre.rsquare), 'Color', 'k')
plot(x,mean_numb_avams(:,2),'*b');
plot(fitobjectPost, 'b')
text(-19, 95, num2str(gofPost.rsquare), 'Color', 'b')
plot(linspace(-15,15,7),ones(1,7)*50)
legend('1st-Pre', 'SigPre', '2nd-Post', 'Sig 2nd');
xlabel('Diff. of the stimulation distance between the forearm and forehead [mm]')
ylabel('Perc. of "forearm" responses')
title({['Fit Of Averaged Data (over ' num2str(length(subj)) ' subjects)']; 'Forearm answers - SPOT OF STIMULUS'})
ylim([0 100]);
savefig(f, [folder_path '47-BetweenSubjects_fit_of_averages.fig']);

% PSE & AE bar plots
f = figure(i_fig); i_fig = i_fig + 1;
X = categorical({'1st - Pre';'2nd'});
meanCoefficients = mean(coefficients,3);
stderrCoefficients = std(coefficients,0,3)/sqrt(s);
colori = 'kb';
subplot(1,2,1)
hold on
for k = 1:2
    bar(X(k), meanCoefficients(k,1), colori(k));
    errorbar(X(k),  meanCoefficients(k,1),  stderrCoefficients(k,1), 'k')
end
hold off
ylim([0 15])
ylabel('PSE')
title(['PSE mean over ' num2str(s) ' subjects'])
subplot(1,2,2)
hold on
for k = 1:2
    bar(X(k), meanCoefficients(k,2), colori(k));
    errorbar(X(k),  meanCoefficients(k,2),  stderrCoefficients(k,2), 'k')
end
hold off
ylim([0 15])
ylabel('PSE')
title(['EA mean over ' num2str(s) ' subjects'])
savefig(f, [folder_path '48-BetweenSubjects_PSE_EA.fig']);

end

