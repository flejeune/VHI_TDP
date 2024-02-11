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
function VHI_cdt_order_ca( subj, folder_path )

% Compute variables of interest
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
        ca = zeros(6,1);

        for i = 1:length(Stimolo1) % 56
            % Count number of times the LARGER stimulus is felt on the FOREARM
            if (strcmp(Area1{i}, 'AVAMBRACCIO'))
                Dif = Stimolo1(i) - Stimolo2(i);
                if (risposta(i) == 1)
                    switch Dif
                        case 5;   ca(4) = ca(4)+1;
                        case 10;  ca(5) = ca(5)+1;
                        case 15;  ca(6) = ca(6)+1;
                    end
                elseif (risposta(i) == 2)
                    switch Dif
                        case -15; ca(1) = ca(1)+1;
                        case -10; ca(2) = ca(2)+1;
                        case -5;  ca(3) = ca(3)+1;
                    end
                end
            else
                Dif = Stimolo2(i) - Stimolo1(i);
                if (risposta(i) == 2) % If answer is two
                    switch Dif
                        case 5;   ca(4) = ca(4)+1;
                        case 10;  ca(5) = ca(5)+1;
                        case 15;  ca(6) = ca(6)+1;
                    end
                elseif (risposta(i) == 1)
                    switch Dif
                        case -15; ca(1) = ca(1)+1;
                        case -10; ca(2) = ca(2)+1;
                        case -5;  ca(3) = ca(3)+1;
                    end
                end
            end              
        end

        % Save current variable (current subject, current condition) in the general matrix
        ca_m(:,t,s) = ca;
        clear ca
    end    
    cd ..
end

% Convert in percentage
ca_perc = 100*ca_m/8;

% Compute mean
mean_ca = mean(ca_perc,3);

%-------------------------------SAVE DATA---------------------------------%
save([folder_path 'CdtOrder_CAcount' num2str(s) '.mat'],'ca_m');

save([folder_path 'CdtOrder_CAmean' num2str(s) '.mat'],'mean_ca');

end