function [t, tone] = makeBBtone(srate, carrier_freq, diff_freq, dur, fdur)
%
% create tone of specified frequency, sampling rate, and duration
% play and plot
%
% input arguments:
% srate -- sampling rate in Hz
% carrier_freq  -- carrier frequency of tone in Hz
% diff_freq -- frequency difference between binaural tones in Hz
% dur -- tone duration in sec
% 
% output arguments:
% t -- vector of sample times in secs
% tone -- two vectors of sample values a.u.

if nargin<5
    fdur = 0.05;
end

t = 0:1/srate:dur;
w = 2*pi*[carrier_freq carrier_freq+diff_freq];

tone = zeros(srate*dur+1,2);
for i=1:2
    tone(:,i) = cos(w(i)*t);
end

tone = makeFadeInOut(srate, tone, 2, fdur);

% soundsc(tone, srate);
% plot(t,tone); shg;
tone = tone';

end