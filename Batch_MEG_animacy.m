
%% To calculate the ERFs to test the animacy effect

%% add path
addpath /local_mount/space/crouch/2/Software/fieldtrip-20150923/;
addpath /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/scripts/;
addpath /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/scripts/functions;
ft_defaults;

%% Separate the trials based on verbs' animacy
subIDs = [2,3,4,5,6,7,8,9,11,12,14,17,18,20,21,22,23,24,25,26,27,28,29,31,32,35,36,37,38,39,40,41];
for i=1:length(subIDs)
    subID=subIDs(1,i);
    load(strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/epochs/data_sub',sprintf('%02d',subID)));    
    
    %% only include the MEG data
    cfg=[];
    cfg.channel = 'MEG';
    cfg.latency= [-0.1 1.1];
    data = ft_selectdata(cfg,data_all);
    %trial information: Trl_Nr_inExp, list, block_Nr, trail_Nr_inBlock,
    %stimuli_ID, Verb_cond, verb_animacy, noun_length, noun_concreteness,
    %verb_length, verb_concreteness
    %% separate the trials of the Animate vs. Inanimate conditions
    cfg = [];
    cfg.trials = intersect([find(data.trialinfo(:,6) == 211);find(data.trialinfo(:,6) == 212);find(data.trialinfo(:,6) == 213);find(data.trialinfo(:,6) == 217)],find(data.trialinfo(:,7) == 1));
    HC_Ani = ft_timelockanalysis(cfg,data);
    
    cfg = [];
    cfg.trials = intersect([find(data.trialinfo(:,6) == 211);find(data.trialinfo(:,6) == 212);find(data.trialinfo(:,6) == 213);find(data.trialinfo(:,6) == 217)],find(data.trialinfo(:,7) == 2));
    HC_Inani = ft_timelockanalysis(cfg,data);
    
    cfg = [];
    cfg.trials = intersect([find(data.trialinfo(:,6) == 212);find(data.trialinfo(:,6) == 214);find(data.trialinfo(:,6) == 216);find(data.trialinfo(:,6) == 218)],find(data.trialinfo(:,7) == 1));
    LC_Ani = ft_timelockanalysis(cfg,data);    
    
    cfg = [];
    cfg.trials = intersect([find(data.trialinfo(:,6) == 212);find(data.trialinfo(:,6) == 214);find(data.trialinfo(:,6) == 216);find(data.trialinfo(:,6) == 218)],find(data.trialinfo(:,7) == 2));
    LC_Inani = ft_timelockanalysis(cfg,data);
        
    % combine gradient before grand average
    cfg = [];
    cmb_HC_Ani = ft_combineplanar(cfg,HC_Ani);
    cmb_HC_Inani = ft_combineplanar(cfg,HC_Inani);
    cmb_LC_Ani = ft_combineplanar(cfg,LC_Ani);
    cmb_LC_Inani = ft_combineplanar(cfg,LC_Inani);    
    
    sub_HC_Ani(i).ERF = cmb_HC_Ani;
    sub_HC_Inani(i).ERF = cmb_HC_Inani;  
    sub_LC_Ani(i).ERF = cmb_LC_Ani;
    sub_LC_Inani(i).ERF = cmb_LC_Inani;  
end        
cfg = [];
% cfg.keepindividual = 'yes';
grandERF_HC_Ani = ft_timelockgrandaverage(cfg, sub_HC_Ani(:).ERF);
grandERF_HC_Inani = ft_timelockgrandaverage(cfg, sub_HC_Inani(:).ERF);
grandERF_LC_Ani = ft_timelockgrandaverage(cfg, sub_LC_Ani(:).ERF);
grandERF_LC_Inani = ft_timelockgrandaverage(cfg, sub_LC_Inani(:).ERF);

%baseline correction
cfg = [];
cfg.baseline = [-0.1 0];%relative to verbs
grandERF_HC_Ani = ft_timelockbaseline(cfg, grandERF_HC_Ani);
grandERF_HC_Inani = ft_timelockbaseline(cfg, grandERF_HC_Inani);
grandERF_LC_Ani = ft_timelockbaseline(cfg, grandERF_LC_Ani);
grandERF_LC_Inani = ft_timelockbaseline(cfg, grandERF_LC_Inani);

grandERF_HC_Ani.cfg = rmfield(grandERF_HC_Ani.cfg, 'previous');
grandERF_HC_Inani.cfg = rmfield(grandERF_HC_Inani.cfg, 'previous');
grandERF_LC_Ani.cfg = rmfield(grandERF_LC_Ani.cfg, 'previous');
grandERF_LC_Inani.cfg = rmfield(grandERF_LC_Inani.cfg, 'previous');
outfil = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/ERF/Verbs_animacyERF_grandavg_keeptrials');
save(outfil, 'grandERF_HC_Ani','grandERF_HC_Inani', 'grandERF_LC_Ani','grandERF_LC_Inani');

%% Extract data within N400 time window and conduct stats
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/ERF/Verbs_animacyERF_grandavg_keeptrials';
% LT = {'MEG1512+1513'; 'MEG0132+0133'};
LT = {'MEG0222+0223';'MEG0212+0213';'MEG0132+0133';'MEG0112+0113';'MEG0232+0233';...
    'MEG0242+0243';	'MEG1512+1513';  'MEG0142+0143'; 'MEG1622+1623';	'MEG1612+1613';	'MEG1522+1523';	'MEG1542+1543';	'MEG1532+1533'};
[~, IND, ~] = intersect(grandERF_HC_Ani.label,LT);
timeWind = linspace(-0.1,1.1,1201);
[~,startTime] = min(abs(timeWind-0.3));
[~,endTime] = min(abs(timeWind-0.5));
HC_Ani = grandERF_HC_Ani.individual(:,IND,[startTime:endTime]);
HC_Inani = grandERF_HC_Inani.individual(:,IND,[startTime:endTime]);
LC_Ani = grandERF_LC_Ani.individual(:,IND,[startTime:endTime]);
LC_Inani = grandERF_LC_Inani.individual(:,IND,[startTime:endTime]);
export_data_N400 = [squeeze(mean(mean(HC_Ani,2),3)),squeeze(mean(mean(HC_Inani,2),3)),squeeze(mean(mean(LC_Ani,2),3)),squeeze(mean(mean(LC_Inani,2),3))];
[a,b,c] = ttest(export_data_N400(:,1),export_data_N400 (:,2))
[a,b,c] = ttest(export_data_N400(:,3),export_data_N400 (:,4))
%to export to SPSS
test = export_data_N400 *1.0e+13;


[~, IND, ~] = intersect(grandERF_HC_Ani.label,LT);
timeWind = linspace(-0.1,1.1,1201);
[~,startTime] = min(abs(timeWind-0.5));
[~,endTime] = min(abs(timeWind-1.1));
HC_Ani = grandERF_HC_Ani.individual(:,IND,[startTime:endTime]);
HC_Inani = grandERF_HC_Inani.individual(:,IND,[startTime:endTime]);
LC_Ani = grandERF_LC_Ani.individual(:,IND,[startTime:endTime]);
LC_Inani = grandERF_LC_Inani.individual(:,IND,[startTime:endTime]);
export_data_late = [squeeze(mean(mean(HC_Ani,2),3)),squeeze(mean(mean(HC_Inani,2),3)),squeeze(mean(mean(LC_Ani,2),3)),squeeze(mean(mean(LC_Inani,2),3))];
%to export to SPSS
test = export_data_late*1.0e+13;

timeWind = linspace(-0.1,1.1,1201);
[~,startTime] = min(abs(timeWind-0.5));
[~,endTime] = min(abs(timeWind-0.6));
HC_Ani = grandERF_HC_Ani.individual(:,IND,[startTime:endTime]);
HC_Inani = grandERF_HC_Inani.individual(:,IND,[startTime:endTime]);
LC_Ani = grandERF_LC_Ani.individual(:,IND,[startTime:endTime]);
LC_Inani = grandERF_LC_Inani.individual(:,IND,[startTime:endTime]);
export_data = [squeeze(mean(mean(HC_Ani,2),3)),squeeze(mean(mean(HC_Inani,2),3)),squeeze(mean(mean(LC_Ani,2),3)),squeeze(mean(mean(LC_Inani,2),3))];
%to export to SPSS
test = export_data *1.0e+13;

%% barplots
fig=figure;
mycolor = {'r','b','r','b'};
labs = {'plaus-ani','plaus-inani','anom-ani','anom-inani'};
subplot(2,2,1)
mean_val = mean(export_data_N400,1)
std_val = std(export_data_N400,1)./sqrt(length(export_data_N400));
for i =1:4
    h=bar(i,mean_val(i));    
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
% ylim([0 2e-12])
errorbar(mean_val,std_val,'.');
title ('VerbS, ERF, 300-500ms')
subplot(2,2,2)
mean_val = mean(export_data_late,1)
std_val = std(export_data_late,1)./sqrt(length(export_data_late));
for i =1:4
    h=bar(i,mean_val(i));
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
ylim([0 2e-12])
errorbar(mean_val,std_val,'.')
title ('Verbs, ERF, 500-1100ms')
subplot(2,2,3)
mean_val = mean(export_data,1)
std_val = std(export_data,1)./sqrt(length(export_data));
for i =1:4
    h=bar(i,mean_val(i));
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
ylim([0 2e-12])
errorbar(mean_val,std_val,'.')
title ('Verbs, ERF, 500-600ms')
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/Verb_MEG_ERF_barplots.eps';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-depsc')
close all;

%% plot the waveforms
set(groot, 'DefaultFigureColormap',jet);
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/ERF/Verbs_animacyERF_grandavg;
grandERF_Ani = grandERF_HC_Ani;
grandERF_Ani.avg = (grandERF_HC_Ani.avg+grandERF_LC_Ani.avg)./2;
grandERF_Inani = grandERF_HC_Inani;
grandERF_Inani.avg = (grandERF_HC_Inani.avg+grandERF_LC_Inani.avg)./2;
%context effect
grandERF_HC = grandERF_HC_Ani;
grandERF_HC.avg = (grandERF_HC_Ani.avg+grandERF_HC_Inani.avg)./2;
grandERF_LC = grandERF_LC_Inani;
grandERF_LC.avg = (grandERF_LC_Inani.avg+grandERF_LC_Ani.avg)./2;
grandERF_HCvsLC = grandERF_HC;
grandERF_HCvsLC.avg = -(grandERF_HC.avg - grandERF_LC.avg);


figure; 
cfg = []; 
cfg.fontsize = 6;
cfg.xlim = [-0.1 1.1];
cfg.graphcolor = [1 0 0; 0 0 1]; %red; blue
cfg.layout = 'neuromag306cmb.lay'; cfg.ylim = [-1e-12 3e-12]; 
% cfg.layout = 'neuromag306all.lay'; cfg.ylim = [-1e-12 4e-12]; 
% cfg.layout = 'neuromag306mag.lay'; cfg.ylim = [-1.5e-13 1.5e-13]; 
% ft_multiplotER(cfg,grandERF_Ani,grandERF_Inani)

% cfg.channel = {'MEG0212+0213'}; %left front-temp sensor
cfg.channel = {'MEG1512+1513'}; %left front-temp sensor
% cfg.channel = {'MEG1932+1933'};%left occipital sensor

cfg.graphcolor = ['r','b','m','g']
ft_singleplotER(cfg,grandERF_HC_Ani,grandERF_HC_Inani,grandERF_LC_Ani,grandERF_LC_Inani)

subplot(2,2,1)
cfg.graphcolor = [1 0 0; 0 0 1]; %red; blue
ft_singleplotER(cfg,grandERF_Ani,grandERF_Inani)
title('all: ani vs inani')

subplot(2,2,3)
ft_singleplotER(cfg,grandERF_HC_Ani,grandERF_HC_Inani)
title('HC: ani vs inani')
subplot(2,2,4)
ft_singleplotER(cfg,grandERF_LC_Ani,grandERF_LC_Inani)
title('LC: ani vs inani')

outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/MNE/Figures/MEG_Verb_animacy.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-dpng')
close all;

%% plot the topographies
figure;
cfg = [];
cfg.xlim = [0.3 0.5];
% cfg.xlim = [0.6 1.0];
% cfg.layout = 'neuromag306planar.lay'; cfg.zlim = [-2e-12 6e-12]; 
% cfg.layout = 'neuromag306mag.lay'; cfg.zlim = [-1.5e-13 1.5e-13]; 
cfg.layout = 'neuromag306cmb.lay'; cfg.zlim = [-2e-12 2e-12]; 
% cfg.layout = 'neuromag306all.lay'; cfg.zlim = [-2e-12 2e-12]; 
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.marker = 'off';
cfg.markersymbol = '+';
cfg.highlight='on';
cfg.highlightchannel = {'MEG1512+1513','MEG1932+1933'}; %left front-temp sensor

% HC conditions
grandERF_HC_dif = grandERF_HC_Ani;
grandERF_HC_dif.avg = -(grandERF_HC_Ani.avg-grandERF_HC_Inani.avg);
grandERF_LC_dif = grandERF_LC_Ani;
grandERF_LC_dif.avg = -(grandERF_LC_Ani.avg-grandERF_LC_Inani.avg);
grandERF_dif = grandERF_Ani;
grandERF_dif.avg = -(grandERF_Ani.avg-grandERF_Inani.avg);

cfg.zlim = [-0.5e-12 0.5e-12]; 
figure;
subplot(3,4,1);
ft_topoplotER(cfg,grandERF_dif);
title('all: Ani-Inani');
cfg.highlight='off';
subplot(3,4,2);
ft_topoplotER(cfg,grandERF_HC_dif);
title('HC: Ani-Inani');
subplot(3,4,3);
ft_topoplotER(cfg,grandERF_LC_dif);
title('LC: Ani-Inani');
subplot(3,4,5);
ft_topoplotER(cfg,grandERF_HCvsLC);
title('HC-LC');
cfg.xlim = [0.5 1.1];
subplot(3,4,7);
ft_topoplotER(cfg,grandERF_dif);
title('Ani-Inani: 500-1100ms');
subplot(3,4,8);
ft_topoplotER(cfg,grandERF_HCvsLC);
title('HC-LC: 500-1100ms');

figure;
cfg = [];
cfg.xlim = [0:0.1:1.1];
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [-1.5e-12 1.5e-12];
cfg.colorbar = 'no';
cfg.comment = 'xlim';
cfg.marker  = 'off';
ft_topoplotER(cfg,grandERF_dif);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/MNE/Figures/MEG_Verb_animacy_topo_grad.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-dpng')
close all;

figure;
cfg = [];
cfg.xlim = [0:0.1:1.1];
cfg.layout = 'neuromag306mag.lay';
cfg.ylim = [-1.5e-13 1.5e-13];
cfg.colorbar = 'no';
cfg.comment = 'xlim';
cfg.marker  = 'off';
ft_topoplotER(cfg,grandERF_dif);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/MNE/Figures/MEG_Verb_animacy_topo_mag.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-dpng')
close all;


%% calculate the cross-temporal spatial correlations
grandERF_HC = grandERF_HC_Ani;
grandERF_HC.avg = (grandERF_HC_Ani.avg+grandERF_HC_Inani.avg)./2;
grandERF_LC = grandERF_LC_Ani;
grandERF_LC.avg = (grandERF_LC_Ani.avg+grandERF_LC_Inani.avg)./2;
grandERF_all = grandERF_HC_Ani;
grandERF_all.avg = (grandERF_HC.avg+grandERF_LC.avg)./2;

data = grandERF_all.avg;
timeWind = linspace(-0.1,1.1,1201); %sampling rate: 500Hz
startWind=1;
endWind=length(timeWind);
for timepoint = startWind : endWind
    refData = data(:,timepoint);
    corrData = data;
    R(timepoint-startWind+1,:) = Correlation(refData,corrData);
    clear ref* corr*
end

%% plot: auto-correlations for the difference ERPs
set(groot, 'DefaultFigureColormap',jet);
timeWind = linspace(-100,1100,1201); %sampling rate: 500Hz
figure('units','normalized','outerposition',[0 0 1 1]);
clims = [-1 1];
imagesc(timeWind,fliplr(timeWind),flipud(R),clims)
axis xy;
% colorbar
xlim(ylim)
set(gca, 'XTick', -100:100:1100, 'YTick', -100:100:1100)
grid minor
outfil = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/MNE/Figures/matrix_difERFs_Ani_Inani');
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-dpng')


%% grand average with keeptrials
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/ERF/Verbs_animacyERF_grandavg_keeptrials';
inter_Ani = grandERF_HC_Ani;
inter_Ani.individual = grandERF_HC_Ani.individual - grandERF_LC_Ani.individual;
inter_Inani = grandERF_HC_Inani;
inter_Inani.individual = grandERF_HC_Inani.individual - grandERF_LC_Inani.individual;

grandERPbase_Ani = grandERF_HC_Ani;
grandERPbase_Ani.individual =  (grandERF_HC_Ani.individual + grandERF_LC_Ani.individual)./2;

grandERPbase_Inani = grandERF_HC_Ani;
grandERPbase_Inani.individual =  (grandERF_HC_Inani.individual + grandERF_LC_Inani.individual)./2;

grandERPbase_HC = grandERF_HC_Ani;
grandERPbase_HC.individual =  (grandERF_HC_Ani.individual + grandERF_HC_Inani.individual)./2;

grandERPbase_LC = grandERF_HC_Ani;
grandERPbase_LC.individual =  (grandERF_LC_Ani.individual + grandERF_LC_Inani.individual)./2;

%% cluster based permutation test
% load neuromag306planar_neighb.mat;

%% build neighourhood file
load neuromag306mag_neighb.mat;
new = neighbours;
for i=1:length(neighbours)
    label = neighbours(i).label;
    new(i).label = strcat(label(1:end-1),'2+',label(4:end-1),'3');
    neib = neighbours(i).neighblabel;
    neib2=neib;
    for j=1:length(neib)
        neiblabel = neib{j};
        neib2{j} = strcat(neiblabel(1:end-1),'2+',neiblabel(4:end-1),'3');
    end
    new(i).neighblabel = neib2;
    clear neib neib2 label
end

cfg = [];
cfg.neighbours = new;
cfg.channel = 'MEG*3';
cfg.method = 'montecarlo';
cfg.design           = [1:32 1:32; ones(1,32), ones(1,32) * 2];
cfg.numrandomization = 1000;
cfg.correctm         = 'cluster';
cfg.clusteralpha            = 0.025;
cfg.clustertail             = 0; % one-or two-sided testing
cfg.clusterstatistic = 'maxsum'; % maximum sum of t-values within one cluster is the test statistic
cfg.alpha     = 0.05;
cfg.tail      = 0; % two-sided testing;
cfg.correcttail      = 'prob';
cfg.statistic = 'depsamplesT';
cfg.uvar             = 1; % subject number (unit variable) on line 1 of the design matrix
cfg.ivar             = 2; % condition number (independent variable) on line 2 of the design matrix
cfg.minnbchan        = 2;
cfg.latency = [0 1.1];
cfg.avgovertime = 'no';

data1 = grandERPbase_HC;
data2 = grandERPbase_LC;
stat = ft_timelockstatistics(cfg, data1, data2);
stat.cfg = rmfield(stat.cfg, 'previous');
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/Verbs_animacyERP_MEG_clusterStats_mainContext';
save(outfil, 'stat');


data1 = grandERPbase_Ani;
data2 = grandERPbase_Inani;
stat = ft_timelockstatistics(cfg, data1, data2);
stat.cfg = rmfield(stat.cfg, 'previous');
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/Verbs_animacyERP_MEG_clusterStats_mainAni';
save(outfil, 'stat');

data1 = inter_Ani;
data2 = inter_Inani;
stat = ft_timelockstatistics(cfg, data1, data2);
stat.cfg = rmfield(stat.cfg, 'previous');
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/Verbs_animacyERP_MEG_clusterStats_interaction';
save(outfil, 'stat');

stat.posclusters(1).prob
stat.negclusters(1).prob

% find out the time window and channels that show significant effects
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/Verbs_animacyERP_MEG_clusterStats_mainAni;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/Verbs_animacyERP_MEG_clusterStats_mainContext;
timeWind = stat.time;
sigtime = timeWind(find(mean(stat.mask,1)> 2/length(stat.label)));
sigtime(1)
sigtime(end)
chan = stat.label;
sigchan = chan(find(mean(stat.mask(:,find(mean(stat.mask,1)> 2/length(stat.label))),2)> 0.8))

%when there are multiple clusters
timeWind = stat.time;
pickcluster=find(stat.posclusterslabelmat~=2);
mask = stat.mask;
mask(pickcluster)=0;
sigtime = timeWind(find(mean(mask,1)> 2/length(stat.label)));
sigtime(1)-1.1
sigtime(end)-1.1
chan = stat.label;
sigchan = chan(find(mean(mask(:,find(mean(mask,1)> 2/length(stat.label))),2)> 0.3))

%plot significant clusters
clims=[-1 1];
imagesc(linspace(0,1100,1201),1:13,flipud(stat.negclusterslabelmat),clims)
set(gca, 'YTick', [1:13], 'YTickLabel', flipud(LT))
outfil = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/MNE/Figures/Verbs_animacyERP_MEG_clusterStats_mainContext.png');
print(outfil,'-dpng')

%% plot topos with significant channels
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/ERF/Verbs_animacyERF_grandavg;
grandERF_Ani = grandERF_HC_Ani;
grandERF_Ani.avg = (grandERF_HC_Ani.avg+grandERF_LC_Ani.avg)./2;
grandERF_Inani = grandERF_HC_Inani;
grandERF_Inani.avg = (grandERF_HC_Inani.avg+grandERF_LC_Inani.avg)./2;
%context effect
grandERF_HC = grandERF_HC_Ani;
grandERF_HC.avg = (grandERF_HC_Ani.avg+grandERF_HC_Inani.avg)./2;
grandERF_LC = grandERF_LC_Inani;
grandERF_LC.avg = (grandERF_LC_Inani.avg+grandERF_LC_Ani.avg)./2;
grandERF_HCvsLC = grandERF_HC;
grandERF_HCvsLC.avg = -(grandERF_HC.avg - grandERF_LC.avg);

figure;
cfg = [];
cfg.xlim = [sigtime(1) sigtime(end)];
cfg.layout = 'neuromag306cmb.lay'; cfg.zlim = [-1e-12 1e-12]; 
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.marker = 'off';
cfg.highlight='on';
cfg.highlightchannel = sigchan;

figure;
subplot(3,4,1);
ft_topoplotER(cfg,grandERF_HCvsLC);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/Verb_MEG_ERF_topos.eps';
print(outfil,'-depsc')