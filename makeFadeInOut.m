function output_tone = makeFadeInOut(srate, input_tone, nchnl, fdur)
%
% create tone of specified frequency, sampling rate, and duration
% play and plot
%
% input arguments:
% srate -- sampling rate in Hz
% input_tone  -- a vector of input tone
% nchnl -- number of channels
% fdur -- fade duration in sec
% 
% output arguments:
% output_tone -- a vector of output tone

fadeLength = length(0:1/srate:fdur);
fadeIn = (1/fadeLength:1/fadeLength:1)';
fadeOut = (1:-1/fadeLength:1/fadeLength)';

output_tone = zeros(length(input_tone),nchnl);
for j=1:nchnl
    output_tone(1:fadeLength,j) = fadeIn.*input_tone(1:fadeLength,j);
    output_tone(fadeLength+1:end-fadeLength,j) = input_tone(fadeLength+1:end-fadeLength,j);
    output_tone(end-fadeLength+1:end,j) = fadeOut.*input_tone(end-fadeLength+1:end,j);
end

end