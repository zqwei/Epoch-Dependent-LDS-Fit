% This analysis is for the temporal correlation of the single neurons at 
% each stage in Li's dataset
% Uni_Latent_dynamic_model

flist        = dir('ANM219033_20131116.mat');
mean_type    = 'Constant_mean';
tol          = 1e-5;
cyc          = 10000;
alpha        = 0.05;
timePoint    = [12 38 64];
xDim         = 3;
binsize      = 50/1000;

fname     = flist.name;
load(fname);
[~,fname, ~] = fileparts(fname);     
Y         = Spikes_correct;
YL        = LSpikes_correct;
valid_trial_type_now = trial_type_correct;
[yDim, T] = size(Y{1});
x_max     = ceil(T*binsize*10)/10;
K         = length(Y);
Y         = cell2mat(Y);
Y         = reshape(Y,yDim,K, T);
Y         = permute(Y, [1 3 2]);
Y         = Y/binsize;
Y         = sqrt(Y);
%     Y(Y==0)   = 1;
%     Y         = log(Y.^2);

m         = ceil(sqrt(yDim));

h         = figure;
for nNeuron = 1: yDim
    subplot(m, m, nNeuron);
    hold on;
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Y(nNeuron,:,valid_trial_type_now)),2),std(squeeze(Y(nNeuron,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Y(nNeuron,:,~valid_trial_type_now)),2),std(squeeze(Y(nNeuron,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
    box off;
    y_S = get(gca,'ylim');
    y_S = y_S(2);
    plot([0.600 0.600],[0 y_S],'--k');
    plot([1.900 1.900],[0 y_S],'--k');
    plot([3.200 3.200],[0 y_S],'--k');
    set(gca,'fontsize',10);
    ylim([0, y_S])
    xlim([0, x_max])
    xlabel('Time (ms)','fontsize',12);
    ylabel('Firing rate (Hz)','fontsize',12);
    hold off;
end

is_suc = false;

while ~is_suc
    try
        Ph      = lds_uni_latent(Y, xDim, 'timePoint', timePoint, 'mean_type',mean_type,'tol',tol,'cyc',cyc,'is_fix_C',false,'is_fix_A',true);
        is_suc  = true;
    catch
        is_suc  = false;
    end
end

[err, y_est, rand_y] = loo (Y, Ph, [0, timePoint, T]);

figure;
for nNeuron = 1: yDim
    subplot(m, m, nNeuron);
    hold on;
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(y_est(nNeuron,:,valid_trial_type_now)),2),std(squeeze(y_est(nNeuron,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(y_est(nNeuron,:,~valid_trial_type_now)),2),std(squeeze(y_est(nNeuron,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
    box off;
    y_S = get(gca,'ylim');
    y_S = y_S(2);
    plot([0.600 0.600],[0 y_S],'--k');
    plot([1.900 1.900],[0 y_S],'--k');
    plot([3.200 3.200],[0 y_S],'--k');
    set(gca,'fontsize',10);
    ylim([-10, y_S])
    xlim([0, x_max])
    xlabel('Time (ms)','fontsize',12);
    ylabel('Firing rate (Hz)','fontsize',12);
    hold off;
end
% % % % setPrintSvg(h, 8*m, 6*m, [fname,'_Est_Dim_3_Firing_rate.svg']);
% % % 
timePoints = [0, timePoint, T];

figure;
for nt_now = 1:4
    timePeriod = timePoints(nt_now) + 1 : timePoints(nt_now+1);
    plot3(squeeze(Ph.Xk_t(1,timePeriod,valid_trial_type_now)),squeeze(Ph.Xk_t(2,timePeriod,valid_trial_type_now)),squeeze(Ph.Xk_t(3,timePeriod,valid_trial_type_now)),'-','color',[0 0 0.25*nt_now])
    hold on; plot3(squeeze(Ph.Xk_t(1,timePeriod,~valid_trial_type_now)),squeeze(Ph.Xk_t(2,timePeriod,~valid_trial_type_now)),squeeze(Ph.Xk_t(3,timePeriod,~valid_trial_type_now)),'-','color',[0.25*nt_now 0 0])
end
hold off;

figure;
for nxDim = 1:xDim
    subplot(1,xDim,nxDim)
    hold on;
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Ph.Xk_t(nxDim,:,valid_trial_type_now)),2),std(squeeze(Ph.Xk_t(nxDim,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
    shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Ph.Xk_t(nxDim,:,~valid_trial_type_now)),2),std(squeeze(Ph.Xk_t(nxDim,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
    box off;
    hold off;
end

