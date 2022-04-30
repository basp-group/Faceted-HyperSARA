function effChans2Image = util_rearrange_channels(idChannels2Image, nChannelsPerImage)
    nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % num of ouput effective channels: number of images in the estimate model cube
    effChans2Image = cell(nEffChans2Image, 1);
    for iEff = 1:nEffChans2Image
        if iEff < nEffChans2Image
            effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
        else
            effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
        end
        fprintf('\nINFO:Effective channel ID %d: physical channels involved: %d - %d', ...
                iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
    end
