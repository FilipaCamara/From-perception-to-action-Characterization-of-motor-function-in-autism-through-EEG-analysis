% Author: Filipa Câmara
% Purpose: Master's thesis in collaboration between the University of Aveiro (UA) and the Institute of Systems and Robotics, University of Coimbra (ISR-UC)
% Date: October, 2025
% Script for reproducible ERSP advanced-visualization.
% Reproducibility notes:
% - Dependencies: MATLAB + EEGLAB + modified tftopo.m
% - Your modified tftopo.m must be first on MATLAB path (see addpath below)
% - Data: per-subject <ID>.dattimef files with fields:
%         times (1 x nT), freqs (1 x nF), trialinfo (struct array with .condition),
%         and channel variables named chan1..chanN (nF x nT x nTrials, complex TF)
% - A template .set is used only to load chanlocs (labels/positions) for topoplots
% - Baseline window is fixed between −500 to 0 ms for dB conversion


%%
close all; clear; clc;

% Folder that contains YOUR custom functions (modified tftopo.m).
addpath('C:\Users\filip\OneDrive\Desktop\Tese\Test_Filipa\Test_Filipa\Códigos_Matlab\', '-begin');

% Folder with <ID>.dattimef files (general data folder)
dataFolder = 'C:\Users\filip\OneDrive\Desktop\Tese\Test_Filipa\Test_Filipa\STUDY\STUDY_SOFT';

% Initialize EEGLAB so topoplot and related functions are available; then close GUI.
eeglab; close;
which tftopo   %modified tftopo.m

%%
% Channels used to compute the TF average (ERSP panel). .
centralNames = {'C3','C4','Cz'};   

% Plot time windows (ms) and color limits (dB), per phase
WIN_instruction = [0 3000];   CLIM_instruction = [-2.5 2.5];
WIN_execution   = [3500 7500]; CLIM_execution   = [-4 4];

% Time–frequency points drawn by your tftopo (time[ms], freq[Hz])
TFPTS_instruction = [950 9; 1850 9; 2500 9]; %edit to frequency points of interest
TFPTS_execution   = [5200 9; 6500 9; 7000 9];

% Subject lists used per phase and group (as specified in dissertation)
SUBJ_instruction_P = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14'};
SUBJ_instruction_S = {'S2','S3','S4','S7','S8','S10','S11','S13','S14','S15','S16','S17','S19','S20'};
SUBJ_execution_P   = {'P1','P2','P4','P5','P7','P10','P11'};
SUBJ_execution_S   = {'S2','S3','S4','S11','S16','S19','S20'};

% ONLY to obtain channel locations (labels/positions for topoplot).
% The chosen file must exist in dataFolder and contain valid EEG.chanlocs.
EEG = pop_loadset('P1_dance_pair_postICAsoft.set', dataFolder);
fullLocs  = EEG.chanlocs;
allLabels = {fullLocs.labels}';
nChans    = numel(fullLocs);

% Finds index for central channels
centralIdx = find(ismember(allLabels, centralNames));


%% .dattimef files
% Get unique trialinfo.condition values for all subjects.
allSubs = [SUBJ_instruction_S, SUBJ_instruction_P, SUBJ_execution_S, SUBJ_execution_P];
conditions = {};
for i = 1:numel(allSubs)
    f = fullfile(dataFolder, allSubs{i} + ".dattimef");
    if ~exist(f,'file'), continue; end
    S = load(f,'-mat');
    if isfield(S,'trialinfo')
        for k = 1:numel(S.trialinfo)
            cnd = S.trialinfo(k).condition;
            if ~any(strcmpi(conditions, cnd))
                conditions{end+1} = cnd; 
            end
        end
    end
end

%% Main loops (CONDITION × PHASE × GROUP)
for ic = 1:numel(conditions)
    condition = conditions{ic};

    % Phase: instruction
    phaseName = 'instruction';
    tWin  = WIN_instruction;  clim = CLIM_instruction;  tfpts = TFPTS_instruction;

    % Group S
    subjects = SUBJ_instruction_S;
    [tf_avg, times, freqs] = compute_group_db(subjects, dataFolder, condition, nChans);
    figure('Position',[100 100 1200 800],'Color','w');
    tftopo(tf_avg, times, freqs, 'chanlocs', fullLocs, 'showchan',0, ...
        'selchans',centralIdx, 'mapchans',1:nChans, 'mode','ave', ...
        'limits',[tWin(1) tWin(2) 4 40 clim(1) clim(2)], 'maplimits',clim, ...
        'timefreqs',tfpts, 'cbar','on', 'logfreq','native', ...
        'title', sprintf('Group S (%s) — %s', condition, phaseName));
    set(findall(gcf,'-property','FontSize'),'FontSize',16);

    % Group P
    subjects = SUBJ_instruction_P;
    [tf_avg, times, freqs] = compute_group_db(subjects, dataFolder, condition, nChans);
    figure('Position',[100 100 1200 800],'Color','w');
    tftopo(tf_avg, times, freqs, 'chanlocs', fullLocs, 'showchan',0, ...
        'selchans',centralIdx, 'mapchans',1:nChans, 'mode','ave', ...
        'limits',[tWin(1) tWin(2) 4 40 clim(1) clim(2)], 'maplimits',clim, ...
        'timefreqs',tfpts, 'cbar','on', 'logfreq','native', ...
        'title', sprintf('Group P (%s) — %s', condition, phaseName));
    set(findall(gcf,'-property','FontSize'),'FontSize',16);

    %Phase: execution 
    phaseName = 'execution';
    tWin  = WIN_execution;  clim = CLIM_execution;  tfpts = TFPTS_execution;

    % Group S
    subjects = SUBJ_execution_S;
    [tf_avg, times, freqs] = compute_group_db(subjects, dataFolder, condition, nChans);
    figure('Position',[100 100 1200 800],'Color','w');
    tftopo(tf_avg, times, freqs, 'chanlocs', fullLocs, 'showchan',0, ...
        'selchans',centralIdx, 'mapchans',1:nChans, 'mode','ave', ...
        'limits',[tWin(1) tWin(2) 4 40 clim(1) clim(2)], 'maplimits',clim, ...
        'timefreqs',tfpts, 'cbar','on', 'logfreq','native', ...
        'title', sprintf('Group S (%s) — %s', condition, phaseName));
    set(findall(gcf,'-property','FontSize'),'FontSize',16);

    % Group P
    subjects = SUBJ_execution_P;
    [tf_avg, times, freqs] = compute_group_db(subjects, dataFolder, condition, nChans);
    figure('Position',[100 100 1200 800],'Color','w');
    tftopo(tf_avg, times, freqs, 'chanlocs', fullLocs, 'showchan',0, ...
        'selchans',centralIdx, 'mapchans',1:nChans, 'mode','ave', ...
        'limits',[tWin(1) tWin(2) 4 40 clim(1) clim(2)], 'maplimits',clim, ...
        'timefreqs',tfpts, 'cbar','on', 'logfreq','native', ...
        'title', sprintf('Group P (%s) — %s', condition, phaseName));
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
end

%% Helper function 
function [tf_avg, times, freqs] = compute_group_db(subjects, dataFolder, condition, nChans)
% Reads <ID>.dattimef, filters trials by 'condition', converts to dB using
% baseline −500 to 0 ms, computes per-channel ERSP, then averages across subjects.
% Assumes presence of fields: times, freqs, trialinfo, chan1..chanN.

    first = true; nSubs = numel(subjects);
    tf_all = []; times = []; freqs = [];

    for s = 1:nSubs
        f = fullfile(dataFolder, subjects{s} + ".dattimef");
        if ~exist(f,'file'), continue; end
        M = load(f,'-mat');

        % Build logical mask for trials of the requested condition
        idxTr = false(1, numel(M.trialinfo));
        for t = 1:numel(M.trialinfo)
            idxTr(t) = strcmpi(M.trialinfo(t).condition, condition);
        end
        if ~any(idxTr), continue; end

        if first
            times = M.times; freqs = M.freqs;
            [~, i0] = min(abs(times - (-500)));   % baseline start 
            [~, i1] = min(abs(times - 0));        % baseline end   
            [nF, nT] = deal(numel(freqs), numel(times));
            tf_all = nan(nF, nT, nChans, nSubs, 'single');
            first = false;
        end

        % Compute ERSP (dB) per channel for this subject
        for c = 1:nChans
            fld = sprintf('chan%d',c);
            if ~isfield(M,fld), continue; end
            X = M.(fld);                               % [nF x nT x nTrials], complex TF
            pow  = mean(abs(X(:,:,idxTr)).^2, 3, 'omitnan');   % trial-average power
            base = mean(pow(:, i0:i1), 2, 'omitnan');          % baseline spectrum
            tf_all(:,:,c,s) = 10*log10(pow ./ base);           % dB normalisation
        end
    end

    % Average across subjects (ignores NaNs if some subjects missing)
    tf_avg = mean(tf_all, 4, 'omitnan');
end
