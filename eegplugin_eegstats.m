% eegplugin_eegstats() - Plugin to return a variety of measures on EEG
%                        or STUDY data
%
% Usage:
%   >> eegplugin_eegstats(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Authors: Arnaud Delorme

function vers = eegplugin_eegstats(fig, trystrs, catchstrs)

vers = 'eegstats1.0'; % write the name of your plugin here and version
if nargin < 3
    error('eegplugin_eegstats requires 3 arguments');
end

if ~exist('restingIAF')
    p = which('eegplugin_eegstats');
    p = fileparts(p);
    addpath(p);
    addpath(fullfile(p, 'restingIAF'));
end

% find DIPFIT menu handle
tools_m = findobj(fig, 'tag', 'tools'); % find by tag

% we create the menu below
cb = [ 'try, [~,~,~,LASTCOM] = pop_eegstats(EEG);' catchstrs.add_to_hist  ];
uimenu( tools_m, 'label', 'EEG stats', 'CallBack', cb, 'separator', 'on');


