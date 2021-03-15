function [inputfiles,outputfiles] = jp_addnoise_hwk(soundfiles, cfg)
%JP_ADDNOISE Adds noise to some soundfiles at specific SNRs.
%
% JP_ADDNOISE(SOUNDFILES, CFG) loops through .wav files in a cell array and
% adds noise to each based on settings in CFG:
%
%   CFG.noisefile  full path to file with noise
%   CFG.prestim    how much noise before stimulus (seconds) [default .5]
%   CFG.poststim   how much noise after stimulus (seconds) [default .5]
%   CFG.snrs       SNRs used to add signal and noise (dB)
%   CFG.fs         sampling frequency
%   CFG.noisedur   noise_only duration when a file is empty [default .5]
%
% If SOUNDFILES is a directory, all of the .wav files in that directory are
% treated as the input files.
%
% The noise and signal files must have the same sampling rate, and are
% assumed to be mono. If the noise file is not long enough to match the
% signal, things will break.
%
% The level of the target (signal) files is unchanged; the noise level is
% adjusted to arrive at different SNRs.
%
% From https://github.com/jpeelle/jp_matlab
% 
% 09/16/19 -- the code has been modified to be used within the experimental
%     code. noise sound to be added is randomly selected from the whole
%     range of noise wav file -- HWK
% 01/29/20 -- 'inputfiles' output has been added. This variable includes
%     input sound files, the rms of which is matched


if ~isfield(cfg, 'prestim') || isempty(cfg.prestim)
    cfg.prestim = 0.5;
end

if ~isfield(cfg, 'poststim') || isempty(cfg.poststim)
    cfg.poststim = 0.5;
end

if ~isfield(cfg, 'noisedur') || isempty(cfg.noisedur)
    cfg.noisedur = 2.56;
end

if ~isfield(cfg, 'fs') || isempty(cfg.fs)
    error('Must specify CFG.fs');
end

if ~isfield(cfg, 'snrs')
    error('Must specify CFG.snrs');
end

if ~isfield(cfg, 'noisefile')
    error('Must specify path to noise file in CFG.noisefile');
end

% error checking
if ~exist(cfg.noisefile, 'file')
    error('Noise file %s not found.', cfg.noisefile);
end

% if soundfiles is a directory, get .wav files
if ischar(soundfiles) && exist(soundfiles, 'dir')
   soundDir = soundfiles;
   D = dir(fullfile(soundfiles, '*.wav'));
   soundfiles = {D.name};

   for i=1:length(soundfiles)
       soundfiles{i} = fullfile(soundDir, soundfiles{i});
   end
   
   % read audio files
    temp = cell(size(soundfiles)); 
    for i = 1:length(soundfiles)
        [temp{i}, fs] = audioread(soundfiles{i}); 
    end
    
    soundfiles = temp;
    
else
    % if string, make a cell
    if ischar(soundfiles)
        soundfiles = cellstr(soundfiles);
    end
    
end

if ~iscell(soundfiles)
    tmp=cell(1,1);
    tmp{1} = soundfiles;
    soundfiles = tmp;
end

% Get noise
[yNoise, fsNoise] = audioread(cfg.noisefile);

ynew_rms = nan(length(soundfiles),length(cfg.snrs));
inputfiles = cell(length(soundfiles),1);  nsignals = length(inputfiles);
outputfiles = cell(length(soundfiles),length(cfg.snrs));

% Loop through soundfiles and add noise
for i = 1:nsignals
    % add noise if there exists an audio file
    fs=cfg.fs;
    if ~isempty(soundfiles{i})
        y = soundfiles{i};  
        inputfiles{i} = y;

        rmsSignal = jp_rms(y);
        dbSignal = jp_mag2db(rmsSignal);

        % get the part of noise we need, and it's RMS and dB
        numSampleNoise = length(y) + cfg.prestim*fs + cfg.poststim*fs;
        spNoise = ceil( rand*(length(yNoise)-numSampleNoise) );

        tmpNoise = yNoise(spNoise:spNoise+numSampleNoise-1);
        tmpNoise = makeFadeInOut(fsNoise, tmpNoise, 1, .05);
        rmsNoise = jp_rms(tmpNoise);

        j=0;
        for thisSNR = cfg.snrs  % going through SNRs
            j=j+1;

            targetDb = dbSignal - thisSNR; % target for noise dB
            targetRMS = 10^(targetDb/20);
            scaleFactor = targetRMS/rmsNoise;

            scaledNoise = tmpNoise * scaleFactor;
            yNew = [zeros(cfg.prestim*fs,1); ... 
                y; zeros(cfg.poststim*fs,1)] + scaledNoise;

            outputfiles{i,j} = yNew;
            ynew_rms(i,j) = jp_rms(yNew);
        end
    else  % add babble track only
        % get the part of noise we need
        numSampleNoise = round((cfg.noisedur + cfg.prestim + ... 
            cfg.poststim) * fs);
        spNoise = ceil( rand*(length(yNoise)-numSampleNoise) );

        tmpNoise = yNoise(spNoise:spNoise+numSampleNoise-1);
        yNew = makeFadeInOut(fsNoise, tmpNoise, 1, .05);
        
        inputfiles{i} = yNew;
        for j = 1:length(cfg.snrs)
            outputfiles{i,j} = yNew;
            ynew_rms(i,j) = jp_rms(yNew);
        end
    end
end

% rematch the rms of input and output soundfiles
mean_rms = mean(ynew_rms(:));
for i = 1:nsignals
    y = inputfiles{i};
    rmsSignal = jp_rms(y);
    g = mean_rms/rmsSignal;
    y2 = y*g;
    inputfiles{i} = y2;
    for j = 1:length(cfg.snrs)
        y = outputfiles{i,j};
        g = mean_rms/ynew_rms(i,j);
        y2 = y*g;
        outputfiles{i,j} = y2;
    end
end

% find a maximum absolute value across signals
maxY = 0;
for i = 1:nsignals
    maxY = max( maxY, max(abs(inputfiles{i})) );
    for j = 1:length(cfg.snrs)
        maxY = max( maxY, max(abs(outputfiles{i,j})) );
    end
end

% finally, normalize all the signals by the maximum value
for i = 1:nsignals
    inputfiles{i} = inputfiles{i}/maxY;
    for j = 1:length(cfg.snrs)
        outputfiles{i,j} = outputfiles{i,j}/maxY;
    end
end

% double-check the clipping issue
for i = 1:nsignals
    if max(inputfiles{i}) > 1 || min(inputfiles{i}) < -1
        fprintf('Clear File %d: MIN = %.3f, MAX = %.3f\n', i, ... 
            min(inputfiles{i}), max(inputfiles{i}));
    end
    for j = 1:length(cfg.snrs)
        if max(outputfiles{i,j}) > 1 || min(outputfiles{i,j}) < -1
            fprintf('Noise File %d: MIN = %.3f, MAX = %.3f\n', i, ... 
                min(outputfiles{i,j}), max(outputfiles{i,j}));
        end
    end
end

% add silent duration prior and post to input audio files to make durations
% of clear and noisy signals equal
for i = 1:nsignals
    inputfiles{i} = [zeros(cfg.prestim*fs,1); inputfiles{i}; ...
        zeros(cfg.poststim*fs,1)]; 
end
