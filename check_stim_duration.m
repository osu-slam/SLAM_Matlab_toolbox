%% check_stim_duration
% Report the minimum, mean, median, and maximum stimuli duration. Simply
% feed it directory information (i.e. the output of the dir function) and
% watch it go!
%
% Author - MJH
% 
% MM/DD/YY -- CHANGELOG
% 12/17/20 -- File initialized in R2017a. MJH

function [dmin, dmean, dmed, dmax] = check_stim_duration(all_files)
%% Check input
if ~isstruct(all_files)
    flag = 1; 
end

temp = {'name', 'folder', 'date', 'bytes', 'isdir', 'datenum'}; 

if ~isempty(setdiff(fields(all_files), temp))
    flag = 1; 
end
clear temp

if flag
    error('You did not feed me directory information! See instructions.')
end

%% Get duration for each
folders = {all_files.folder}'; names = {all_files.name}'; 
fnames = fullfile(folders, names); 

durations = nan(size(fnames)); 

for ff = 1:length(fnames)
    tempinfo = audioinfo(fnames{ff});    
    durations(ff) = tempinfo.Duration; 
end

clear tempinfo

if any(isnan(durations))
    error('Something broke, likely I could not read duration info!')
end

%% Calculate minimum, mean, median, and max stim duration
dmin = min(durations);
dmean = mean(durations);
dmedian = median(durations);
dmax = max(durations);

end


