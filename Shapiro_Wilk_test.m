%% Normality using Shapiro Wilk (8–10 Hz) — Instruction + Execution
clear; clc; close all;


customFolder = "C:\Users\filip\OneDrive\Desktop\Tese\Test_Filipa\Test_Filipa\Códigos_Matlab\";
addpath(customFolder, '-begin');
dataFolder = 'C:\Users\filip\OneDrive\Desktop\Tese\Test_Filipa\Test_Filipa\STUDY\STUDY_SOFT';


eeglab; close;

%ROI & Frequency
roiChans   = {'C3','Cz','C4'};
freqWin_Hz = [8 10];
base_ms    = [-1000 0];

%Time windows
tw_instr1 = [750 1500];
tw_instr2 = [1500 3000];
tw_exec   = [3500 7500];

%Subjects
instr_P = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14'};
instr_S = {'S2','S3','S4','S7','S8','S10','S11','S13','S14','S15','S16','S17','S19','S20'};
exec_P  = {'P1','P2','P4','P5','P7','P10','P11'};
exec_S  = {'S2','S3','S4','S11','S16','S19','S20'};

%Conditions
conds_instr = struct( ...
   'label', {'Walk Confident','Walk Sad','Dance Solo','Dance Pair'}, ...
   'key',   {'walk_confident_soft','walk_sad_soft','dance_solo_soft','dance_pair_soft'} );

conds_exec = struct( ...
   'label', {'Walk Confident','Walk Natural','Walk Sad','Dance Solo','Dance Pair','Dance Shake'}, ...
   'key',   {'walk_confident_soft','walk_natural_soft','walk_sad_soft','dance_solo_soft','dance_pair_soft','dance_shake_soft'} );

%Channel labels → indices
EEG = pop_loadset('P1_dance_pair_postICAsoft.set', dataFolder);
allLabels = string({EEG.chanlocs.labels});
roiIdx = find(ismember(upper(allLabels), upper(string(roiChans))));

%Time/Freq axes
probe = load(fullfile(dataFolder, exec_P{1} + ".dattimef"), '-mat');
times = probe.times; freqs = probe.freqs;
[~, iB0] = min(abs(times - base_ms(1)));
[~, iB1] = min(abs(times - base_ms(2)));
fMask = freqs >= freqWin_Hz(1) & freqs <= freqWin_Hz(2);

%Colors
blue = hex2rgb('#1d5b90'); burgundy = hex2rgb('#ac1e04');


alpha = 0.05;

%% 
% INSTRUCTION — two time windows (2+2)


timeWins = {tw_instr1, tw_instr2};
for t = 1:2
    timeWin_ms = timeWins{t};
    [~, iT0] = min(abs(times - timeWin_ms(1)));
    [~, iT1] = min(abs(times - timeWin_ms(2)));

    %Group P
    valsP = cell(1,numel(conds_instr)); statsP = repmat(struct('W',NaN,'p',NaN),1,numel(conds_instr));
    for c = 1:numel(conds_instr)
        valsP{c} = extract_vals(dataFolder,instr_P,conds_instr(c).key,roiIdx,iB0,iB1,iT0,iT1,fMask);
        [~,statsP(c).p,statsP(c).W] = swtest(valsP{c},alpha);
    end

    %Group S
    valsS = cell(1,numel(conds_instr)); statsS = repmat(struct('W',NaN,'p',NaN),1,numel(conds_instr));
    for c = 1:numel(conds_instr)
        valsS{c} = extract_vals(dataFolder,instr_S,conds_instr(c).key,roiIdx,iB0,iB1,iT0,iT1,fMask);
        [~,statsS(c).p,statsS(c).W] = swtest(valsS{c},alpha);
    end

    %Plot group P
    titleP = sprintf('Group P — %d–%d ms (ROI %s, %g–%g Hz)', ...
        timeWin_ms(1), timeWin_ms(2), strjoin(roiChans,'/'), freqWin_Hz(1), freqWin_Hz(2));
    make_grid_QQ(valsP,conds_instr,blue,burgundy,[titleP ' — QQ'],statsP);
    make_grid_HIST(valsP,conds_instr,blue,burgundy,[titleP ' — Histograms'],statsP);

    %Plot group S
    titleS = sprintf('Group S — %d–%d ms (ROI %s, %g–%g Hz)', ...
        timeWin_ms(1), timeWin_ms(2), strjoin(roiChans,'/'), freqWin_Hz(1), freqWin_Hz(2));
    make_grid_QQ(valsS,conds_instr,blue,burgundy,[titleS ' — QQ'],statsS);
    make_grid_HIST(valsS,conds_instr,blue,burgundy,[titleS ' — Histograms'],statsS);
end

%% 
%EXECUTION — 3500-7500 ms (3+3) 

timeWin_ms = tw_exec;
[~, iT0] = min(abs(times - timeWin_ms(1)));
[~, iT1] = min(abs(times - timeWin_ms(2)));

%Group P
valsP = cell(1,numel(conds_exec)); statsP = repmat(struct('W',NaN,'p',NaN),1,numel(conds_exec));
for c = 1:numel(conds_exec)
    valsP{c} = extract_vals(dataFolder,exec_P,conds_exec(c).key,roiIdx,iB0,iB1,iT0,iT1,fMask);
    [~,statsP(c).p,statsP(c).W] = swtest(valsP{c},alpha);
end
titleP = sprintf('Group P — Execution (%d–%d ms, ROI %s)', timeWin_ms(1), timeWin_ms(2), strjoin(roiChans,'/'));
make_grid_QQ(valsP,conds_exec,blue,burgundy,[titleP ' — QQ'],statsP);
make_grid_HIST(valsP,conds_exec,blue,burgundy,[titleP ' — Histograms'],statsP);

%Group S
valsS = cell(1,numel(conds_exec)); statsS = repmat(struct('W',NaN,'p',NaN),1,numel(conds_exec));
for c = 1:numel(conds_exec)
    valsS{c} = extract_vals(dataFolder,exec_S,conds_exec(c).key,roiIdx,iB0,iB1,iT0,iT1,fMask);
    [~,statsS(c).p,statsS(c).W] = swtest(valsS{c},alpha);
end
titleS = sprintf('Group S — Execution (%d–%d ms, ROI %s)', timeWin_ms(1), timeWin_ms(2), strjoin(roiChans,'/'));
make_grid_QQ(valsS,conds_exec,blue,burgundy,[titleS ' — QQ'],statsS);
make_grid_HIST(valsS,conds_exec,blue,burgundy,[titleS ' — Histograms'],statsS);


%% ================== HELPER FUNCTIONS ==================
function vals = extract_vals(dataFolder, groupIDs, condKey, roiIdx, iB0, iB1, iT0, iT1, fMask)
    vals = nan(numel(groupIDs),1);
    for s = 1:numel(groupIDs)
        file = fullfile(dataFolder, groupIDs{s} + ".dattimef");
        if ~exist(file,'file'), continue; end
        M = load(file,'-mat');
        if ~isfield(M,'trialinfo'), continue; end
        idx = strcmp({M.trialinfo.condition},condKey);
        if ~any(idx), continue; end
        roiVals = nan(numel(roiIdx),1);
        for k = 1:numel(roiIdx)
            fld = sprintf('chan%d',roiIdx(k));
            if ~isfield(M,fld), continue; end
            X = M.(fld);
            P = mean(abs(X(:,:,idx)).^2,3,'omitnan');
            B = mean(P(:,iB0:iB1),2,'omitnan');
            E = 10*log10(P./B);
            ROI = E(fMask,iT0:iT1);
            roiVals(k) = mean(ROI(:),'omitnan');
        end
        vals(s) = mean(roiVals,'omitnan');
    end
end

function make_grid_QQ(vals_all, conds, blue, burgundy, bigTitle, stats)
    n = numel(conds); ncols = ceil(n/2);
    figure('Color','w','Position',[100 100 1200 800]);
    tl = tiledlayout(2,ncols,'TileSpacing','compact','Padding','compact');
    title(tl,bigTitle,'FontWeight','bold','FontSize',13);
    for i = 1:n
        nexttile; x = vals_all{i}; x = x(~isnan(x));
        if isempty(x), axis off; title(conds(i).label); continue; end
        h = qqplot(x); grid on; box off;
        if startsWith(conds(i).key,'walk'), c=blue; else, c=burgundy; end
        set(h(1),'MarkerFaceColor',c,'MarkerEdgeColor',c);
        set(h(2),'LineWidth',1.8,'Color',[0.3 0.3 0.3]);
        xlabel('Theoretical Quantiles'); ylabel('Sample Quantiles');
        title(conds(i).label);
        Wi = stats(i).W; Pi = stats(i).p;
        text(0.98,0.02,sprintf('W=%.3f, p=%.3f',Wi,Pi),...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',9);
    end
end

function make_grid_HIST(vals_all, conds, blue, burgundy, bigTitle, stats)
    n = numel(conds); ncols = ceil(n/2);
    figure('Color','w','Position',[100 100 1200 800]);
    tl = tiledlayout(2,ncols,'TileSpacing','compact','Padding','compact');
    title(tl,bigTitle,'FontWeight','bold','FontSize',13);
    for i = 1:n
        nexttile; hold on; x = vals_all{i}; x = x(~isnan(x));
        if isempty(x), axis off; title(conds(i).label); continue; end
        if startsWith(conds(i).key,'walk'), c=blue; else, c=burgundy; end
        histogram(x,'Normalization','pdf','FaceColor',c,'EdgeColor',c,'FaceAlpha',0.35);
        [fhat,xi] = ksdensity(x); plot(xi,fhat,'LineWidth',1.6,'Color',c);
        mu=mean(x); sd=std(x);
        xx=linspace(min(x)-sd,max(x)+sd,300);
        yy=(1/(sd*sqrt(2*pi)))*exp(-0.5*((xx-mu)/sd).^2);
        plot(xx,yy,'--','Color',[0.2 0.2 0.2],'LineWidth',1.2);
        legend({'Hist','KDE','Normal fit'},'Location','best'); legend boxoff;
        title(conds(i).label);
        Wi = stats(i).W; Pi = stats(i).p;
        text(0.98,0.02,sprintf('W=%.3f, p=%.3f',Wi,Pi),...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',9);
    end
end

function rgb = hex2rgb(hex)
    if hex(1)=='#', hex=hex(2:end); end
    rgb=[hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))]/255;
end
