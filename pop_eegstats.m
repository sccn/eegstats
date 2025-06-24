% pop_eegstats() - fit multiple component dipoles using DIPFIT
%
% Usage:
%         >> pop_eegstats(EEG); % pop-up graphical interface
%         >> [EEG, eegstats] = pop_eegstats(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG      - input EEGLAB dataset.
%
% Optional inputs:
% 'thetarange'     - [min max] theta range. Default is [4 6].
% 'alpharange'     - [min max] alpha range. Default is [8 12].
% 'otherranges'    - [min max] alpha range. Default is [8 12].
% 'averagepower'   - ['on' or 'off'] average power across channels. Default is 'off'
% 'channels'       - [list of channels] can be a cell array of string or
%                    integer array. Default is empty and to use all
%                    channels.
% 'winsize'        - [float] window size in seconds for power density estimate.
%                    Default is 2 seconds.
% 'overlap'        - [float] overlap in seconds for power density estimate.
%                    Default is 2 seconds.
% 'csvfile'        - ['string'] file name to save results. Default is to
%                    only display results on the command line.
%
% Individual alpha frequency options:
% 'iaf'            - ['on' or 'off'] save/display individual alpha frequencies
%                    Default is 'on'.
% 'iafminchan'     - ['integer'] minimum number of channels to compute 
%                    individual alpha frequency. Sometimes it is not
%                    possible to find an alpha peak for some channels.
%
% Alpha asymmetry options:
% 'alphaasymmetry' - ['on' or 'off'] save/display alpha asymmetry
%                    (calculated by taking the difference of 10*log10(power) 
%                    in the alpha frequency range defined in 'alpharange'
%                    between the two channels defined in 'asymchans').
%                    Default is 'on'.
% 'asymchans'       - [list of two channels] can be a cell array of string or
%                    integer array. Default F7,F8 or F3,F4 or FP1,FP2 if 
%                    present in the dataset.
%
% Additional parameters for restingIAF function (see its help message).
% 'Fw' (default 11), 'k' (default 5), 'mpow', 'mdiff', 'nfft', 'norm'
% 'taper' (all with default the same as the restingIAF function).
%
% Outputs:
%   EEG     - EEG structure with EEG.etc.eegstats structure containing fields
%        eegstats.power, spectral power computed using the pwelch function
%                               of size frequencies x channels
%        eegstats.powerfreqs, frequency ranges for the power above
%        eegstats.restingIAF.freqs, trimmed vector of frequency bins resolved 
%                               by the PWELCH MATLAB function
%        eegstats.restingIAF.pSum, structure containing summary statistics of 
%                               alpha-band parameters
%        eegstats.restingIAF.iafChan, structure containing channel-wise spectral 
%                               and alpha parameter data
%        eegstats.alpha_asymmetry, alpha power asymetry between elelectrodes 
%                               defined in 'asymchans' at frequency 'alpharange'
%  eegstats - structure with the same fields as above
%
% See also: RESTINGIAF
%
% Author: Arnaud Delorme, Oct. 2021

% Copyright (C) 2021 Arnaud Delorme, SCCN/INC/UCSD
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [EEG, eegstats, com] = pop_eegstats(EEG, varargin)

if nargin < 1
    help pop_eegstats;
    return;
end

com = '';
eegstats = [];
if nargin<2
    
    chanLabel = lower({ EEG(1).chanlocs.labels });
    defaultAssymetry = '';
    if ~isempty( strmatch( 'f7', chanLabel, 'exact') ) && ~isempty( strmatch( 'f8', chanLabel, 'exact') )
        defaultAssymetry = 'F7 F8';
    elseif ~isempty( strmatch( 'f3', chanLabel, 'exact') ) && ~isempty( strmatch( 'f4', chanLabel, 'exact') )
        defaultAssymetry = 'F3 F4';
    elseif ~isempty( strmatch( 'fp1', chanLabel, 'exact') ) && ~isempty( strmatch( 'fp2', chanLabel, 'exact') )
        defaultAssymetry = 'FP1 FP2';
    end
    
    cb_chan1 = 'tmpEEG = get(gcbf, ''userdata''); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''channels''), ''string'',tmpval); clear tmp tmpEEG tmpchanlocs tmpval';
    cb_chan2 = 'tmpEEG = get(gcbf, ''userdata''); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''asymchans''), ''string'',tmpval); clear tmp tmpEEG tmpchanlocs tmpval';
    cb_file1 = [ '[filename, filepath] = uiputfile(''*'', ''Enter a file name'');' ...
                    'if filename(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''csvfile''), ''string'', fullfile(filepath, filename) );' ...
                    'end;' ...
                    'clear filename filepath;' ];
    uilist = { ...
        { 'style' 'text' 'string' 'Frequency ranges' 'fontweight' 'bold' } ...
        {} { 'style' 'text' 'string' 'Theta frequency range (min max in Hz)' } ...
        { 'style' 'edit' 'string' '4 6' 'tag' 'thetarange' } {} ...
        {} { 'style' 'text' 'string' 'Alpha frequency range (min max in Hz)' } ...
        { 'style' 'edit' 'string' '8 12' 'tag' 'alpharange'  } {} ...
        {} { 'style' 'text' 'string' 'Other frequency ranges (separate by ";")' } ...
        { 'style' 'edit' 'string' '18 22; 30 45' 'tag' 'otherranges'  } {} ...
        ...
        { 'style' 'text' 'string' 'Power Spectral Density parameters (PSD)' 'fontweight' 'bold' } ...
        {} { 'style' 'checkbox' 'string' 'Average channels'' power (unchecked: individual channel power)' 'tag' 'averagepower' } ...
        {} {} ...
        {} { 'style' 'text' 'string' 'Channel indices or names (default:all)' } ...
        { 'style' 'edit' 'string' '' 'tag' 'channels' } ...
        { 'style' 'pushbutton' 'string'  '...', 'enable' fastif(isempty(EEG(1).chanlocs), 'off', 'on') 'callback' cb_chan1 }, ...
        {} { 'style' 'text' 'string' 'PSD window size (seconds)' } ...
        { 'style' 'edit' 'string' '2' 'tag' 'winsize' } {}  ...
        {} { 'style' 'text' 'string' 'PSD window overlap (seconds)' } ...
        { 'style' 'edit' 'string' '1' 'tag' 'overlap'  } {}  ...
        ...
        { 'style' 'text' 'string' 'Other measures' 'fontweight' 'bold' } ...
        {} { 'style' 'checkbox' 'string' 'Individual alpha freq. with at least' 'tag' 'iaf' 'value' 1 } ...
        { 'style' 'edit' 'string' '1' 'tag' 'iafminchan'} { 'style' 'text' 'string' 'channel(s)' } ...
        {} { 'style' 'checkbox' 'string' 'Alpha asymmetry between' 'tag' 'alphaasymmetry' 'value' fastif(isempty(defaultAssymetry),0,1) } ...
        { 'style' 'edit' 'string' defaultAssymetry 'tag' 'asymchans' } ...
        { 'style' 'pushbutton' 'string'  '...', 'enable' fastif(isempty(EEG(1).chanlocs), 'off', 'on') 'callback' cb_chan2 } ...
        { } ...
        { 'style' 'text' 'string' 'Export to CSV file (optional)' 'fontweight' 'bold' } { } ...
        { 'style' 'edit' 'string' '' 'tag' 'file1' 'tag' 'csvfile' } ...
        { 'style' 'pushbutton' 'string'  '...', 'callback' cb_file1 } ...
        };
    
    row = [0.2 1.9 0.8 0.6] ;
    [~,~,~,results] = inputgui( 'geometry', { [1] row row row 1 [0.2 3.1 0.05 0.05] row row row 1 row row 1 [1.5 0.1 1.1 0.6] }, ...
        'uilist', uilist, 'helpcom', 'pophelp(''pop_eegstats'')', ...
        'title', 'Compute EEG measures -- pop_eegstats()', 'userdata', EEG(1));
    if isempty(results), return; end

    results.thetarange = str2num(results.thetarange);
    results.alpharange = str2num(results.alpharange);
    results.otherranges = str2num( [ '[' results.otherranges ']' ]);
    results.winsize = str2double(results.winsize);
    results.overlap = str2double(results.overlap);
    results.iaf = fastif(results.iaf, 'on', 'off');
    results.iafminchan = str2double(results.iafminchan);
    results.alphaasymmetry = fastif(results.alphaasymmetry, 'on', 'off');
    results.averagepower   = fastif(results.averagepower  , 'on', 'off');
    
    options = fieldnames(results);
    options(:,2) = struct2cell(results);
    options = options';
    options = options(:)';
else
    options = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    % check that the dipfit settings are the same
    [ EEG, com ] = eeg_eval( 'pop_eegstats', EEG, 'params', options );
    return;
end

% checking parameters
% -------------------
g = finputcheck(options, { ...
    'thetarange'     'float'  []      [4 6];
    'alpharange'     'float'  []      [8 12];
    'otherranges'    'float'  []      [];
    'averagepower'   'string' { 'on' 'off' } 'off';
    'channels'       {'integer' 'string' 'cell'} { {} {} {} } '';
    'winsize'        'float'  []      2;
    'overlap'        'float'  []      1;
    'mpow'           'float'  []      1;
    'Fw'             'float'  []      11;
    'k'              'float'  []      5;
    'mdiff'          'float'  []      0.2;
    'iaf'            'string'  { 'on' 'off' } 'on';
    'nfft'           'integer' []     [];
    'norm'           'boolean' []     false;
    'taper'          'string'  {}     'hamming';
    'iafminchan'     'integer' { }    1;
    'alphaasymmetry'          'string'  { 'on' 'off' } 'on';
    'asymchans'       {'integer' 'string' 'cell'} { {} {} {} } '';
    'csvfile'     'string'   { } '' });

if isstr(g), error(g); end
%EEG     = eeg_checkset(EEG, 'chanlocs_homogeneous');
if ~isempty(g.csvfile)
    [~,~,ext] = fileparts(g.csvfile);
    if ~isequal(lower(ext), '.csv')
        g.csvfile = [ g.csvfile '.csv' ];
    end
end

if isempty(g.channels), g.channels = 1:EEG.nbchan; end
[chaninds g.channels] = eeg_decodechan(EEG.chanlocs, g.channels);
g.asymchans = eeg_decodechan(EEG.chanlocs, g.asymchans);
if length(g.asymchans) ~= 2 && ~isempty(g.asymchans)
    error('Channel to calculate asymetry not found or too many/too few channels');
end
if isempty(g.nfft), g.nfft = 2^nextpow2(g.winsize*EEG.srate); end
    
[iafSum,  iafChan,  freqs] = restingIAF(EEG.data(chaninds,:), length(chaninds), g.iafminchan, [0 EEG(1).srate/2], EEG(1).srate, g.alpharange, g.Fw, g.k, g.mpow, g.mdiff, g.taper  , g.winsize*EEG(1).srate, g.overlap*EEG(1).srate, g.nfft, g.norm);

freqRanges = [ g.thetarange; g.alpharange; g.otherranges];
if ~isempty(g.csvfile)
    fid = fopen(g.csvfile, 'w');
    if fid == -1
        error('Cannot open file for writing');
    end
else
    fid = [];
end

% channel labels on header row
myfprintf(fid, '\tAll_channels');
if ~strcmpi(g.averagepower, 'on')
    for iChan = 1:length(g.channels)
        myfprintf(fid, '\t%s', num2str(g.channels{iChan}));
    end
end
myfprintf(fid, '\n');

% write frequencies
power = [];
for iFreq = 1:size(freqRanges,1)
    myfprintf(fid, '%2.1f-%2.1f Hz', freqRanges(iFreq,1), freqRanges(iFreq,2));
    [~,indBeg] = min(abs(freqs-freqRanges(iFreq,1)));
    [~,indEnd] = min(abs(freqs-freqRanges(iFreq,2)));
    for iChan = 1:length(g.channels)
        power(iFreq, iChan) = mean(iafChan(iChan).pxx(indBeg:indEnd));
    end
    myfprintf(fid, '\t%1.5f', 10*log10(mean(power(iFreq,:))));
    if ~strcmpi(g.averagepower, 'on')
        myfprintf(fid, '\t%1.5f', 10*log10(power(iFreq,:)));
    end
    myfprintf(fid, '\n');
end

% write IAF
if strcmpi(g.iaf, 'on')
    myfprintf(fid, 'Peak alpha frequency');
    myfprintf(fid, '\t%1.2f', iafSum.paf);
    if ~strcmpi(g.averagepower, 'on')
        for iChan = 1:length(g.channels)
            myfprintf(fid, '\t%1.2f', iafChan(iChan).peaks);
        end
    end
    myfprintf(fid, '\n');
    myfprintf(fid, 'Alpha center of gravity');
    myfprintf(fid, '\t%1.2f', iafSum.cog);
    if ~strcmpi(g.averagepower, 'on')
        for iChan = 1:length(g.channels)
            myfprintf(fid, '\t%1.2f', iafChan(iChan).gravs);
        end
    end
    myfprintf(fid, '\n');
end

% write alpha asymetry
if strcmpi(g.alphaasymmetry, 'on') && ~isempty(g.asymchans)
    [~,indBeg] = min(abs(freqs-g.alpharange(1)));
    [~,indEnd] = min(abs(freqs-g.alpharange(2)));
    
    % See Smith et al. (2017) for log-transformation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6449497/
    power1 = 10*log10(mean(iafChan(g.asymchans(1)).pxx(indBeg:indEnd)));
    power2 = 10*log10(mean(iafChan(g.asymchans(2)).pxx(indBeg:indEnd)));
    
    alpha_asymmetry = power1-power2;
    myfprintf(fid, 'Alpha Asymmetry\t%1.5f\n', power1-power2);
    myfprintf(fid, '\n');
else
    alpha_asymmetry = NaN;
end

% close file
if ~isempty(fid)
    fclose(fid);
end

eegstats.power        = power;
eegstats.powerfreqs   = freqRanges;
eegstats.restingIAF.freqs   = freqs;
eegstats.restingIAF.pSum    = iafSum;
eegstats.restingIAF.iafChan = iafChan;
eegstats.alpha_asymmetry = alpha_asymmetry;
eegstats.parameters      = options;

EEG.etc.eegstats = eegstats;
disp('EEG statistics saved into EEG.etc.eegstats')

% write history
com = sprintf('EEG = pop_eegstats(EEG, %s);', vararg2str(options));

function myfprintf(fid, str, varargin)

str = sprintf(str, varargin{:});
fprintf('%s', str);
if ~isempty(fid)
    str(str == 9) = ',';
    fprintf(fid, '%s', str);
end
