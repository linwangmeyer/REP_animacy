%% all the sensors: combine the RSA effects for the two datasets: nonconERP and nonconMM

%% add path
addpath /local_mount/space/crouch/2/Software/fieldtrip-20150923/;
addpath /autofs/cluster/kuperberg/nonconMM/nonconOSC/scripts/;
addpath /autofs/cluster/kuperberg/nonconMM/nonconOSC/scripts/functions;
ft_defaults;

%% For the verbs: all Trials
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani'; %nonconMM: 500HZ, -0.5 - 2.1s; 1301 time points
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_verb_all'; %nonconERP: 512HZ, -0.5 - 2.1s; 1332 time points
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_all';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;

%% For the verbs: only HC trials
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_HC_AniInani';
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_verb_HC';
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_HC';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;

%% For the verbs: only LC trials
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_LC_AniInani';
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_verb_LC';
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_LC';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%% plot barplots
load (strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_all'));
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_HC';
sub_HC_Ani =  Ani;
sub_HC_Inani =  Inani;
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_LC';
sub_LC_Ani =  Ani;
sub_LC_Inani =  Inani;

timeWind = linspace(-0.5,2.1,1301);
[~,StartSample] = min(abs(timeWind-0.45));%0.6
[~,EndSample] = min(abs(timeWind-0.65));%1.0
HC_Ani = mean(sub_HC_Ani(:,[StartSample:EndSample]),2);
HC_Inani = mean(sub_HC_Inani(:,[StartSample:EndSample]),2);
LC_Ani = mean(sub_LC_Ani(:,[StartSample:EndSample]),2);
LC_Inani = mean(sub_LC_Inani(:,[StartSample:EndSample]),2);

export_data_450to650 = [HC_Ani,HC_Inani,LC_Ani,LC_Inani];

fig=figure;
mycolor = {'r','b','r','b'};
subplot(2,2,1)
mean_val = mean(export_data_450to650,1)
std_val = std(export_data_450to650,1)./sqrt(length(export_data_450to650));
for i =1:4
    h=bar(i,mean_val(i));    
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
% ylim([0 2e-12])
errorbar(mean_val,std_val,'.');
title ('Verbs, RSA, 450-650ms')

%Other time windows: 300-500ms
timeWind = linspace(-0.5,2.1,1301);
[~,StartSample] = min(abs(timeWind-0.3));
[~,EndSample] = min(abs(timeWind-0.5));
HC_Ani = mean(sub_HC_Ani(:,[StartSample:EndSample]),2);
HC_Inani = mean(sub_HC_Inani(:,[StartSample:EndSample]),2);
LC_Ani = mean(sub_LC_Ani(:,[StartSample:EndSample]),2);
LC_Inani = mean(sub_LC_Inani(:,[StartSample:EndSample]),2);

export_data = [HC_Ani,HC_Inani,LC_Ani,LC_Inani];

subplot(2,2,3)
mean_val = mean(export_data,1)
std_val = std(export_data,1)./sqrt(length(export_data));
for i =1:4
    h=bar(i,mean_val(i));
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
ylim([0 0.04])
errorbar(mean_val,std_val,'.')
title ('Verbs, RSA, 300-500ms')

%Other time windows: 500-1100ms
timeWind = linspace(-0.5,2.1,1301);
[~,StartSample] = min(abs(timeWind-0.5));
[~,EndSample] = min(abs(timeWind-1.1));
HC_Ani = mean(sub_HC_Ani(:,[StartSample:EndSample]),2);
HC_Inani = mean(sub_HC_Inani(:,[StartSample:EndSample]),2);
LC_Ani = mean(sub_LC_Ani(:,[StartSample:EndSample]),2);
LC_Inani = mean(sub_LC_Inani(:,[StartSample:EndSample]),2);

export_data = [HC_Ani,HC_Inani,LC_Ani,LC_Inani];

subplot(2,2,4)
mean_val = mean(export_data,1)
std_val = std(export_data,1)./sqrt(length(export_data));
for i =1:4
    h=bar(i,mean_val(i));
    setcor = mycolor{i};
    set(h,'FaceColor',setcor);
    hold on
end
ylim([0 0.04])
errorbar(mean_val,std_val,'.')
title ('Verbs, RSA, 500-1100ms')

outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/Verb_EEG_RSA_barplots.eps';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(outfil,'-depsc')
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the nouns: Plausible trials
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_221245;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_221245;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_221245';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the nouns: Plausible but unpredicted trials
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_21245';
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_22245;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_22245';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;



%% For the nouns: HC_pred
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_211;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_211;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_211';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;

%% For the nouns: HC_lexiviol
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_212;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_212;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_212';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;

%% For the nouns: LC_unpred
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_2145;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_2145;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2145';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the nouns: Implausible trials
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_223768;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_223768;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_223768';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%% For the nouns: HC_SRviol
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_2137;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_2137;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2137';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%% For the nouns: LC_SRviol
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/RSA_spatial_AniInani_noun_2168;
Ani = Rsub_Ani;
Inani = Rsub_Inani;
Btw = Rsub_Btw;
clear Rsub*

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/nonconERP_RSA_spatial_AniInani_noun_2168;
for i = 1:size(Rsub_Ani,1)
    Ani(i+32,:) = resample(Rsub_Ani(i,:),size(Ani,2),size(Rsub_Ani,2)); %make sure they have the same time points: 500Hz, 1301 time points
    Inani(i+32,:) = resample(Rsub_Inani(i,:),size(Ani,2),size(Rsub_Ani,2));
    Btw(i+32,:) = resample(Rsub_Btw(i,:),size(Ani,2),size(Rsub_Ani,2));
end
Ani_avg = mean(Ani,1);
Inani_avg = mean(Inani,1);
Btw_avg = mean(Btw,1);
outfil = '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2168';
save(outfil, 'Ani','Inani','Btw','Ani_avg','Inani_avg','Btw_avg');
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the results from here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stats: verbs %%
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_all';
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_HC';
load '/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_verb_LC';
data1_sub =  Ani;
data2_sub = Inani;
data3_sub = Btw;
data1_avg = Ani_avg;
data2_avg = Inani_avg;
data3_avg = Btw_avg;

timeWind = linspace(-0.5,2.1,1301);
[~,StarTime] = min(abs(timeWind-0.45));%0.6
[~,EndTime] = min(abs(timeWind-0.65));%1.0
within1 = mean(data1_sub(:,[StarTime:EndTime]),2);
within2 = mean(data2_sub(:,[StarTime:EndTime]),2);
between = mean(data3_sub(:,[StarTime:EndTime]),2);
[h p ic stat] = ttest(within1,within2)
result(1,1) = stat.tstat;
result(1,2) = p
data=[within1,within2];

%% plot time course
fig = figure;
timeWind = linspace(-0.5,2.1,1301);
plot(timeWind, data1_avg,'r','LineWidth',2);
axis([-0.1 1.1 0 0.1]);
hold on
plot(timeWind, data2_avg,'b','LineWidth',2);
axis([-0.1  1.1 0 0.1]);
ax=gca;
ax.XTick = [-0.1:0.2:1.1];
xt=get(gca,'XTick');
set(gca,'XTick',xt,'XTickLabel',xt*1000);

filename = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/combineDataset_RSA_spatinfo_AniInani_timesource_Verbs_all');
saveas(fig,filename,'epsc')

%% plot scatter
fig = figure;
scatter(within2,within1,'filled','d')
axis([-0.002 0.1 -0.002 0.1]);
title ({'animate vs. inanimate'; strcat('0.45 to 0.65s, p=',num2str(p))})
filename = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/combineDataset_RSA_spatinfo_AniInani_scatter_verb_all.eps');
print(fig,filename,'-depsc')
close all;


%% FDR stats
timeWind = linspace(-0.5,2.1,1301);
[~,StartSample] = min(abs(timeWind-0));
[~,EndSample] = min(abs(timeWind-1.1));
data1 = data1_sub(:,[StartSample:EndSample]);%conv_
data2 = data2_sub(:,[StartSample:EndSample]);%conv_
[h p ci stats] = ttest(data1,data2)

[h,p,ci,stats] = ttest(data1,data2,'alpha',0.05,'dim',1)

[fdrH, crit_p]=fdr_bky(squeeze(p),0.05,'yes');
% [fdrH2, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(p),0.05,'pdep','yes');
plotFig = fdrH;
figure;
plot(timeWind(1,[StartSample:EndSample]),plotFig,'b','LineWidth',2);

%get the significant time points
toi = timeWind(1,[StartSample:EndSample]);
toi(plotFig)


%% Permutation test

%get the original t-values
timeWind = linspace(-0.5,2.1,1301);
[~,StarTime] = min(abs(timeWind-0));
[~,EndTime] = min(abs(timeWind-1.1));

data1 = data1_sub(:,[StartSample:EndSample]);%conv_
data2 = data2_sub(:,[StartSample:EndSample]);%conv_
[h,p,ci,stats] = ttest(data1,data2,'alpha',0.05,'dim',1);
hvalues = squeeze(h);
tvalues = squeeze(stats.tstat);

%find all the clusters and calculate the sum of the t-values
[L,n] = bwlabel(hvalues,4);
for k = 1:n
    if size(find(L==k),1) >= 1& ~isnan(hvalues(find(L==k)))
        cluster(k) = sum(tvalues(find(L==k)));
    end
end
clear h p ci stats k

%permutation test
data_dif=data1-data2;
permN = 10000; %number of permutation
n_subs = 72; %number of subjects
alph = 0.05; %sig threshold
tail = 0;
[Perm,tmx_ptile] = twoDpermute(data_dif,permN,n_subs,alph,tail);


%find and plot the largest cluster
mean(abs(Perm) > max(cluster)) %pos
mean(-abs(Perm) < min(cluster)) %neg
sig_sample = find(L==find(cluster == max(cluster)));
timeWind_roi = timeWind(1,[StarTime:EndTime]);
timeWind_roi(sig_sample)


%plot the permutation distribution
figure;
hist(Perm)

%compare the observed clusters with the distribution, two tails
plotFig = L;
%pos cluster
indClusterPos = find(cluster>=tmx_ptile(2));
for m = 1:length(indClusterPos)
    highlight = find(L==indClusterPos(m));
    plotFig(highlight) = cluster(indClusterPos(m));
    pvaluePos(m) = mean(abs(Perm) > cluster(indClusterPos(m)));
end

%neg cluster
indClusterNeg = find(cluster<=tmx_ptile(1));
for m = 1:length(indClusterNeg)
    highlight = find(L==indClusterNeg(m));
    plotFig(highlight) = cluster(indClusterNeg(m));
    pvalueNeg(m) = mean(-abs(Perm) < cluster(indClusterNeg(m)));
end


%% stats: nouns
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_221245;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_22245;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_223768;

load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_211;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_212;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2145;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2137;
load /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/EEG_RSA/combineDataset_RSA_spatial_AniInani_noun_2168;

data1_sub =  Ani;
data2_sub = Inani;
data3_sub = Btw;
data1_avg = Ani_avg;
data2_avg = Inani_avg;
data3_avg = Btw_avg;

timeWind = linspace(-0.5,2.1,1301);
[~,StarTime] = min(abs(timeWind-1.4));%0.3
[~,EndTime] = min(abs(timeWind-1.6));%0.5
% [~,StarTime] = min(abs(timeWind-1.7));%0.6
% [~,EndTime] = min(abs(timeWind-2.1));%1.0
within1 = mean(data1_sub(:,[StarTime:EndTime]),2);
within2 = mean(data2_sub(:,[StarTime:EndTime]),2);
between = mean(data3_sub(:,[StarTime:EndTime]),2);
[h p ic stat] = ttest(within1,within2)
result(1,1) = stat.tstat;
result(1,2) = p
data=[within1,within2];

%% plot time course
fig = figure;
timeWind = linspace(-0.5,2.1,1301);
plot(timeWind, data1_avg,'r','LineWidth',2);
axis([1.0 2.1 0 0.1]);
hold on
plot(timeWind, data2_avg,'b','LineWidth',2);
axis([1.0  2.1 0 0.1]);
ax=gca;
ax.XTick = [1.1:0.2:2.1];
xt=get(gca,'XTick');
set(gca,'XTick',xt,'XTickLabel',(xt-1.1)*1000);

filename = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/combineDataset_RSA_spatinfo_AniInani_timesource_noun_22245');
saveas(fig,filename,'epsc')


%% FDR stats
timeWind = linspace(-0.5,2.1,1301);
[~,StartSample] = min(abs(timeWind-1.4));
[~,EndSample] = min(abs(timeWind-1.7));
data1 = data1_sub(:,[StartSample:EndSample]);%conv_
data2 = data2_sub(:,[StartSample:EndSample]);%conv_
[h p ci stats] = ttest(data1,data2)

[h,p,ci,stats] = ttest(data1,data2,'alpha',0.05,'dim',1)

[fdrH, crit_p]=fdr_bky(squeeze(p),0.05,'yes');
% [fdrH2, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(squeeze(p),0.05,'pdep','yes');
plotFig = fdrH;
figure;
plot(timeWind(1,[StartSample:EndSample]),plotFig,'b','LineWidth',2);


%get the significant time points
toi = timeWind(1,[StartSample:EndSample]);
toi(plotFig)


%% Permutation test

%get the original t-values
timeWind = linspace(-0.5,2.1,1301);
[~,StarTime] = min(abs(timeWind-1.4));
[~,EndTime] = min(abs(timeWind-1.7));

data1 = data1_sub(:,[StartSample:EndSample]);%conv_
data2 = data2_sub(:,[StartSample:EndSample]);%conv_
[h,p,ci,stats] = ttest(data1,data2,'alpha',0.05,'dim',1);
hvalues = squeeze(h);
tvalues = squeeze(stats.tstat);

%find all the clusters and calculate the sum of the t-values
[L,n] = bwlabel(hvalues,4);
for k = 1:n
    if size(find(L==k),1) >= 1& ~isnan(hvalues(find(L==k)))
        cluster(k) = sum(tvalues(find(L==k)));
    end
end
clear h p ci stats k

%permutation test
data_dif=data1-data2;
permN = 10000; %number of permutation
n_subs = 72; %number of subjects
alph = 0.05; %sig threshold
tail = 0;
[Perm,tmx_ptile] = twoDpermute(data_dif,permN,n_subs,alph,tail);


%find and plot the largest cluster
mean(abs(Perm) > max(cluster)) %pos
mean(-abs(Perm) < min(cluster)) %neg
sig_sample = find(L==find(cluster == max(cluster)));
timeWind_roi = timeWind(1,[StarTime:EndTime]);
timeWind_roi(sig_sample)


%plot the permutation distribution
figure;
hist(Perm)

%compare the observed clusters with the distribution, two tails
plotFig = L;
%pos cluster
indClusterPos = find(cluster>=tmx_ptile(2));
for m = 1:length(indClusterPos)
    highlight = find(L==indClusterPos(m));
    plotFig(highlight) = cluster(indClusterPos(m));
    pvaluePos(m) = mean(abs(Perm) > cluster(indClusterPos(m)));
end

%neg cluster
indClusterNeg = find(cluster<=tmx_ptile(1));
for m = 1:length(indClusterNeg)
    highlight = find(L==indClusterNeg(m));
    plotFig(highlight) = cluster(indClusterNeg(m));
    pvalueNeg(m) = mean(-abs(Perm) < cluster(indClusterNeg(m)));
end