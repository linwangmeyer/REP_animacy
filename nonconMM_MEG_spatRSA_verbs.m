%% nonconMM MEG dataset: spatial RSA for Verbs
%% add path
addpath /local_mount/space/crouch/2/Software/fieldtrip-20150923/;
addpath /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/scripts/;
addpath /autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/scripts/functions;
addpath /run/media/lw604/3T/newMEG/scripts/functions;
ft_defaults;

%% Conduct spatial similarity analysis to MEG data of five conditions relative to verb's onset

%% define trials based on conditions of verbs
trigger = {[211,212,213,217],[212,214,216,218]}; %verb's condition
conds = {'HC','LC'};

for indCond = 1: length(trigger) %loop across conditions
    nsubjects = [1:41];
    nsubjects([1 10 13 15 16 19 30 33 34]) = [];
    for i=1:size(nsubjects,2)
        j = nsubjects(i);
        load(strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/epochs/data_sub',sprintf('%02d',j)));
        
        %% only select data in the selected time window, and select relevant conditions
        for indcond = 1:length(trigger{indCond})
            condTrl{indcond} = find(data_all.trialinfo(:,6) == trigger{indCond}(1,indcond));
        end
        condTrl_all = vertcat(condTrl{:});

        cfg=[];
        cfg.latency = [-0.1 1.1];
        cfg.trials = condTrl_all;
        cfg.channel = 'MEG';
        data = ft_selectdata(cfg,data_all);
        data_sel_all = permute(cat(3,data.trial{:}),[3,1,2]); %get MEG data: trl*sensor*time   
        
        %% calculate the spatial similarity
        R = PearsonCorrSpat_shift_BtwPair(data_sel_all, data_sel_all);
        Rsub(i,:) = R;
        clear data* condTrl_all condTrl index* R
    end
    R_avg = mean(Rsub,1);
    outfil = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/RSA/nonconMM_MEG_RSA_spatial_verb_',conds{indCond});
    save(outfil, 'Rsub','R_avg');
    clear R*
end


%% plot time course
cd /local_mount/space/crouch/1/users/nonconMM/MEG/fieldtrip/MNE/Figures;
conds = {'HC','LC'};

load ('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/RSA/nonconMM_MEG_RSA_spatial_verb_HC');
R_HC = R_avg;
Rsub_HC = Rsub;
clear R_avg;
load ('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/MNE/RSA/nonconMM_MEG_RSA_spatial_verb_LC');
R_LC = R_avg;
Rsub_LC = Rsub;

fig = figure(1);
timeWind = linspace(-0.1,1.1,1201);
plot(timeWind, R_HC,'r','LineWidth',2);
hold on
plot(timeWind, R_LC, 'k','LineWidth',2);
axis([-0.1  1.1 -0.01 0.1]);
ax=gca;
ax.XTick = [-0.1:0.1:1.1];
xt=get(gca,'XTick');
set(gca,'XTick',xt,'XTickLabel',xt*1000);
title('spatial RSA: verbs; n=32')
filename = strcat('/autofs/cluster/kuperberg/nonconMM/MEG/fieldtrip/Figures/nonconMM_MEG_RSA_spatial_verb_timesource.png');
print(fig,filename,'-dpng')


%% stats
timeWind = linspace(-0.1,1.1,1201);
[~,StarTime] = min(abs(timeWind-0.3));
[~,EndTime] = min(abs(timeWind-0.5));
within1 = mean(Rsub_HC(:,[StarTime:EndTime]),2);
within2 = mean(Rsub_LC(:,[StarTime:EndTime]),2);
[h p ic stat] = ttest(within1,within2)
result(1,1) = stat.tstat;
result(1,2) = p
