function [t, tone] = make_tone (srate, freq, dur)
%
% create tone of specified frequency, sampling rate, and duration
% play and plot
%


% input arguments:
% srate -- sampling rate in Hz
% freq  -- frequency of tone in Hz
% 
% output arguments:
% t -- vector of sample times in secs
% tone -- vector of sample values a.u.
% dur = 0.01;

t = 0:1/srate:dur;
w = 2*pi*freq;
fadeLength = length(0:1/srate:0.005);
fadeIn = 1/fadeLength:1/fadeLength:1;
fadeOut = 1:-1/fadeLength:1/fadeLength;
tone = fadeIn.*cos(w*t(1:fadeLength));
tone = [tone cos(w*t(fadeLength+1:end-fadeLength))];
tone = [tone fadeOut.*cos(w*t(end-fadeLength+1:end))];

% soundsc(tone, srate);
% plot(t,tone); shg;
end