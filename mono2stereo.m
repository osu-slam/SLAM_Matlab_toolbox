function audio_out = mono2stereo(audio_in)
% change mono channel inputs into stereo channel outputs
% and transpose 

if iscell(audio_in)
    audio_out = cell(size(audio_in));
    for i=1:size(audio_in,1)
        for j=1:size(audio_in,2)
            alen = size(audio_in{i,j});
            if alen(2)==1 && alen(1)>2
                audio_out{i,j} = [audio_in{i,j} audio_in{i,j}]';
            else
                error('Each audio input should have n x 1 lengths');
            end
        end
    end
else
    alen = size(audio_in);
    if alen(2)==1 && alen(1)>2
        audio_out = [audio_in audio_in]';
    else
        error('Audio input should have n x 1 lengths');
    end
end

end