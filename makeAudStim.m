function [sound_files,dur_info,trial_output] = makeAudStim(stim_dir, type, trial_input)
%
% sentences are randomly assigned to each condition defined based on syntax
% (OO/OS/SO/SS), target gender (M/F), and/or filler gender (M/F) factors
% and generate auditory outputs according to either pre-defined or randomized order
%
% input arguments:
% stim_dir -- directory where auditory stimuli are located
% type -- stim set type -'P'(for practice), 'A', 'B', or 'C' (for main expt)
%         trial order for tutorial is fixed - no need to specify trial_input
% trial_input -- n by 3 string array with columns [ syntax, target_g, filler_g ]
%              each element of syntax column may include OO, OS, SO, or SS
%              each element of target_g column may include M or F
%              each element of filler_g column may include M or F
%              The number of M and F in target_g and filler_g should be coundter-balanced
%              n is the number of sentences and should be 20 for 'P',
%              48 for 'A' or 'B', and 96 for 'C' (both A and B)
% 
% output arguments:
% sound_files -- n by 1 cell of sound outputs 
% dur_info -- n by 1 vector of output durations in sec 
% trial_output -- n by 1 vector of sentence index

if isstring(type)
    type = char(type);
end

if ~(type=='P' || type=='A' || type=='B' || type=='C')
    error('incorrect type specification');
end

tf_gender_AB=0;  tf_gender_P=0;
load('tf_gender_info_v2.mat','tf_gender_AB','tf_gender_P')
switch type
    case 'P',  tf_gender = tf_gender_P;  isB = 0;
    case 'A',  tf_gender = tf_gender_AB(1:48,:);  isB = 0;
    case 'B',  tf_gender = tf_gender_AB(49:96,:);  isB = 1;
    case 'C',  tf_gender = tf_gender_AB;  isB = 0;
end

nsent = size(trial_input,1);
sound_files = cell(nsent,1);
dur_info = zeros(nsent,1);
trial_output = zeros(nsent,1);

% assign sentences randomly to each condition
ref = ["M","F"];
for i=ref
    for j=ref
        outidx = ( trial_input(:,2)==i&trial_input(:,3)==j );
        sentlist = find( tf_gender(:,1)==i&tf_gender(:,2)==j );
        trial_output(outidx) = Shuffle(sentlist) + isB*48;
    end
end

% duration information
dur_info_AB=0;  dur_info_P=0;
load('duration_info_v2.mat','dur_info_AB','dur_info_P')
switch type
    case 'P',  dur_info_tmp = dur_info_P;
    case 'A',  dur_info_tmp = dur_info_AB(1:48,:);
    case 'B',  dur_info_tmp = dur_info_AB(49:96,:);
    case 'C',  dur_info_tmp = dur_info_AB;
end

ref = ["OO","OS","SO","SS"];
for t=1:nsent
    sentidx = trial_output(t) - isB*48;
    synidx = ( ref==trial_input(t,1) );
    dur_info(t) = dur_info_tmp(sentidx,synidx);
end

% generate auditory outputs
for t=1:nsent
    sentidx = num2digit(trial_output(t),2);
    synidx = char(trial_input(t,1));
    if type=='C'
        if trial_output(t)<=48
            stimfile = [stim_dir '/A' sentidx '_' synidx '.wav'];
        else
            stimfile = [stim_dir '/B' sentidx '_' synidx '.wav'];
        end
    else
        stimfile = [stim_dir '/' type sentidx '_' synidx '.wav'];
    end
    
    [sound_files{t}, ~] = audioread(stimfile);
end

end