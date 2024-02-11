%% Authors : 
%   - Marco D'Alonzo, PhD. Senior research associate.
%       marco.dalonzo@unicampus.it
%   - François Le Jeune, PhD. Post-doctoral fellow.
%       francois.le-jeune@hotmail.fr
%
%
% Affiliation of both authors at time of editing : 
%   - NeXT Lab, Università Campus Bio-Medico di Roma (UCBM), Roma, Italy.

%% Main script

subj = {}; % Instert name of folder of participants (ex: {'S1','S2','S3'})
    

% Excel file columns. length(column_name) has to be equal to (length(subj) + 1)
column_name = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' ...
    'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z' 'AA' 'AB' 'AC' 'AD' ...
    'AE' 'AF' 'AG' 'AH' 'AI' 'AJ' 'AK' 'AL' 'AM' 'AN' 'AO' 'AP' 'AQ' ...
    'AR' 'AS' 'AT' 'AU' 'AV' 'AW' 'AX' 'AY' 'AZ' 'BA' 'BB' 'BC' 'BD' ...
    'BE' 'BF' 'BG' 'BH' 'BI' 'BJ' 'BK' 'BL' 'BM' 'BN' 'BO' 'BP' 'BQ' 'BR'};

% Path to the folder containing the folder of the participants
folder_path = '';
if(exist(folder_path, 'dir')==0)
    mkdir(folder_path);
end

% Figures number initialisation
i_fig = 1;
%% 

[i_fig] = VHI_analysis(subj, column_name, folder_path, i_fig);
[i_fig] = VHI_first_plots(length(subj), folder_path, i_fig);
[i_fig] = VHI_fitting(length(subj), folder_path, i_fig);
[i_fig] = VHI_cdt_order_tdpt(subj, folder_path, i_fig);
VHI_cdt_order_ca(subj, folder_path);
[i_fig] = VHI_cdt_order_illusion(subj, folder_path, i_fig);
[i_fig] = VHI_cdt_order_fit( subj, folder_path, i_fig);
[i_fig] = VHI_between_subjects(subj, folder_path, i_fig);
[i_fig] = VHI_between_subjects_fit(subj, folder_path, i_fig);
[i_fig] = VHI_between_subjects_illusion(subj, folder_path, i_fig);
VHI_to_stats(subj, folder_path);
VHI_to_stats_cdt_order(subj, folder_path);
VHI_to_stats_between_subjects(subj, folder_path);
