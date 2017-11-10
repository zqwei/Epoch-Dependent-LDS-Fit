% % % 
% % % 
% % % figure;
% % % 
% % % if xDim == 3
% % %     index_valid_trial_type_yes = find(valid_trial_type_now);
% % %     index_valid_trial_type_no  = find(~valid_trial_type_now);
% % %     
% % %     plot3(squeeze(Ph.Xk_t(1,:,index_valid_trial_type_yes(1:10))),...
% % %         squeeze(Ph.Xk_t(2,:,index_valid_trial_type_yes(1:10))),...
% % %         squeeze(Ph.Xk_t(3,:,index_valid_trial_type_yes(1:10))),'-b')
% % %     
% % %     hold on;
% % %     
% % %     plot3(squeeze(Ph.Xk_t(1,:,index_valid_trial_type_no(1:10))),...
% % %         squeeze(Ph.Xk_t(2,:,index_valid_trial_type_no(1:10))),...
% % %         squeeze(Ph.Xk_t(3,:,index_valid_trial_type_no(1:10))),'-r')
% % %     
% % % end
% % % 
% % % 
% % % %% For error trials
% % % % valid_trial_type_now = trial_type_error;
% % % % Y         = Spikes_error;
% % % % [yDim, T] = size(Y{1});
% % % % x_max     = ceil(T*binsize*10)/10;
% % % % K         = length(Y);
% % % % Y         = cell2mat(Y);
% % % % Y         = reshape(Y,yDim,K, T);
% % % % Y         = permute(Y, [1 3 2]);
% % % % Y         = Y/binsize;
% % % % 
% % % % [err, y_est, ~]  = loo(Y, Ph, [0, timePoint, T]);
% % % % figure;
% % % % for nNeuron = 1: yDim
% % % %     subplot(m, m, nNeuron);
% % % %     hold on;
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(y_est(nNeuron,:,valid_trial_type_now))+Ph.d(nNeuron,1),2),std(squeeze(y_est(nNeuron,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(y_est(nNeuron,:,~valid_trial_type_now))+Ph.d(nNeuron,1),2),std(squeeze(y_est(nNeuron,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
% % % %     box off;
% % % %     y_S = get(gca,'ylim');
% % % %     y_S = y_S(2);
% % % %     plot([0.600 0.600],[0 y_S],'--k');
% % % %     plot([1.900 1.900],[0 y_S],'--k');
% % % %     plot([3.200 3.200],[0 y_S],'--k');
% % % %     set(gca,'fontsize',10);
% % % %     ylim([-10, y_S])
% % % %     xlim([0, x_max])
% % % %     xlabel('Time (ms)','fontsize',12);
% % % %     ylabel('Firing rate (Hz)','fontsize',12);
% % % %     hold off;
% % % % end
% % % % 
% % % % h         = figure;
% % % % for nNeuron = 1: yDim
% % % %     subplot(m, m, nNeuron);
% % % %     hold on;
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Y(nNeuron,:,valid_trial_type_now)),2),std(squeeze(Y(nNeuron,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Y(nNeuron,:,~valid_trial_type_now)),2),std(squeeze(Y(nNeuron,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
% % % %     box off;
% % % %     y_S = get(gca,'ylim');
% % % %     y_S = y_S(2);
% % % %     plot([0.600 0.600],[0 y_S],'--k');
% % % %     plot([1.900 1.900],[0 y_S],'--k');
% % % %     plot([3.200 3.200],[0 y_S],'--k');
% % % %     set(gca,'fontsize',10);
% % % %     ylim([0, y_S])
% % % %     xlim([0, x_max])
% % % %     xlabel('Time (ms)','fontsize',12);
% % % %     ylabel('Firing rate (Hz)','fontsize',12);
% % % %     hold off;
% % % % end
% % % % 
% % % % timePoints = [0, timePoint, T];
% % % % 
% % % % for nt_now   = 1:4
% % % %     Y(:,timePoints(nt_now)+1:timePoints(nt_now+1),:) = ...
% % % %         remove_mean(Y(:,timePoints(nt_now)+1:timePoints(nt_now+1),:),Ph.d(:,nt_now));
% % % % end
% % % % 
% % % % Ph_err   = kalman_forward_backward(Y, Ph.A, Ph.Q, Ph.C, Ph.R, Ph.x0, Ph.Q0, timePoints);
% % % % 
% % % % figure;
% % % % for nxDim = 1:3
% % % %     subplot(1,3,nxDim)
% % % %     hold on;
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Ph_err.Xk_t(nxDim,:,valid_trial_type_now)),2),std(squeeze(Ph_err.Xk_t(nxDim,:,valid_trial_type_now)),[],2),{'-b','linewid',2},0.5);
% % % %     shadedErrorBar(((1:T)-0.5)*binsize,mean(squeeze(Ph_err.Xk_t(nxDim,:,~valid_trial_type_now)),2),std(squeeze(Ph_err.Xk_t(nxDim,:,~valid_trial_type_now)),[],2),{'-r','linewid',2},0.5);
% % % %     box off;
% % % %     hold off;
% % % % end
% % % % 
% % % % % first lick time vs latent x
% % % % xDim  = 2;
% % % % 
% % % % figure;
% % % % for nDim = 1:xDim
% % % %     
% % % %     subplot(1,xDim,nDim)
% % % %     
% % % %     first_lick = cellfun(@min_with_bound,lick_times_correct,num2cell(3.2*ones(size(lick_times_correct))),'uni',false);
% % % %     first_lick = [first_lick{:}];
% % % %     hold on;
% % % %     plot(squeeze(Ph.Xk_t(nDim,64,trial_type_correct)),(first_lick(trial_type_correct)-3.2)*1000, 'ob')
% % % %     plot(squeeze(Ph.Xk_t(nDim,64,~trial_type_correct)),(first_lick(~trial_type_correct)-3.2)*1000, 'or')
% % % %     first_lick = cellfun(@min_with_bound,lick_times_error,num2cell(3.2*ones(size(lick_times_error))),'uni',false);
% % % %     first_lick = [first_lick{:}];
% % % %     plot(squeeze(Ph_err.Xk_t(nDim,64,trial_type_error)),(first_lick(trial_type_error)-3.2)*1000, 'sb')
% % % %     plot(squeeze(Ph_err.Xk_t(nDim,64,~trial_type_error)),(first_lick(~trial_type_error)-3.2)*1000, 'sr')
% % % %     hold off
% % % % end
% % % % 
% % % % figure;
% % % % for nDim = 1:xDim
% % % %     
% % % %     subplot(1,xDim,nDim)
% % % %     last_lick = cellfun(@max,lick_times_correct,'uni',false);
% % % %     last_lick = [last_lick{:}];
% % % %     hold on;
% % % %     plot(squeeze(Ph.Xk_t(nDim,64,trial_type_correct)),(last_lick(trial_type_correct)-3.2)*1000, 'ob')
% % % %     plot(squeeze(Ph.Xk_t(nDim,64,~trial_type_correct)),(last_lick(~trial_type_correct)-3.2)*1000, 'or')
% % % %     last_lick = cellfun(@max,lick_times_error,'uni',false);
% % % %     last_lick = [last_lick{:}];
% % % %     plot(squeeze(Ph_err.Xk_t(nDim,64,trial_type_error)),(last_lick(trial_type_error)-3.2)*1000, 'sb')
% % % %     plot(squeeze(Ph_err.Xk_t(nDim,64,~trial_type_error)),(last_lick(~trial_type_error)-3.2)*1000, 'sr')
% % % %     hold off
% % % % end
% % % % 
% % % % 