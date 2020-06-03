
function Weight_output=MultipleRegFuncCalWeight(Finger,Input_dynamics,Output_Dynamics,Figures)


% inputs: Brain_Area; Raw==1; F1_Raw_Hilbert or F1_Hilbert & F1_Filtered ; F1_F2E_Vels; F1_E2F_Vels

Lagged_Kins_S=Output_Dynamics;
Lagged_Projs=Input_dynamics;
Lag_Number=0;
%%
% smoothing the kinematics because of connections points for trials
%Lagged_Kins_S=smooth(Lagged_Kins,50,'lowess');

% Z score of kinematics ? I dont think it is necessary. We have columns of
% ones

% z score of projectedPCs? The ECoG data already z scored and this projection is generated. 

% adding columns of ones

% Lagged_Projs and Lagged_Kins_S are ready for k-fold 
% getting the average of R2 values

fold_num = 10;
idx_test{1} = 1:floor((1/fold_num)*(length(Lagged_Kins_S))); idx_train{1} = setdiff(1:length(Lagged_Kins_S),idx_test{1});
for i=2:fold_num
    idx_test{i} = idx_test{i-1}(end)+1:floor((i)*(1/fold_num)*(length(Lagged_Kins_S))); idx_train{i}= setdiff(1:length(Lagged_Kins_S),idx_test{i});
end


for i=1:fold_num
    
    Lagged_Projs_ones=[ones(size(Lagged_Projs,1),1),Lagged_Projs];
    Projs_lag_test = Lagged_Projs_ones(idx_test{i},:);
    Projs_lag_train = Lagged_Projs_ones(idx_train{i},:);
    ProjsTrain=Projs_lag_train;
    ProjsTest=Projs_lag_test;
    
    velTrain{i}=Lagged_Kins_S(idx_train{i},:);
    velTest{i}=Lagged_Kins_S(idx_test{i},:);
    
    betaVel{i}=pinv(ProjsTrain'*ProjsTrain)*(ProjsTrain'*velTrain{i});
    velPredict_Train{i}=ProjsTrain*betaVel{i};
    velPredict_Test{i}=ProjsTest*betaVel{i};
end


tmp1=[];
for i=1:fold_num
    tmp1 = [tmp1;(betaVel{i}(:,1)')];
end
betaVelcv(1,:) = mean(tmp1);

betaVelcv = betaVelcv';
Weight_output=betaVelcv;

velPredict_cv = Lagged_Projs_ones*betaVelcv;

% 5*2 plots!!
if Figures==1
    
    if Lag_Number~=0
        figure;
        suptitle(['Finger:',num2str(Finger),'; WLag:',num2str(Lag_Number)]);
    elseif Lag_Number==0
        figure;
        suptitle(['Finger:',num2str(Finger),'; W/OLag'])
        
    end
    
    for i=1:fold_num
        
        subplot(fold_num,2,i*2-1)
        plot(velTrain{i}); hold on; plot(velPredict_Train{i},'r-','linewidth',2);
        R_Pre=corrcoef(velTrain{i},velPredict_Train{i});
        Rsq_Pre=R_Pre(1,2).^2;
        title(['Predicted Train Data: R2: ',num2str(Rsq_Pre)])
        
        subplot(fold_num,2,i*2)
        plot(velTest{i}); hold on; plot(velPredict_Test{i},'r-','linewidth',2);
        R_Pre=corrcoef(velTest{i},velPredict_Test{i});
        Rsq_Pre=R_Pre(1,2).^2;
        title(['Predicted Test data: R2: ',num2str(Rsq_Pre)])
    end
    
end
% One final plot for CV
R_Pre_final=corrcoef(Lagged_Kins_S,velPredict_cv);
Rsq_Pre_output=R_Pre_final(1,2).^2;

if Figures==1
    
    figure;
    set(gcf, 'Position', [100, 100, 2000, 500]);
    plot(Lagged_Kins_S,'b','linewidth',1.5)
    %plot((1:1:size(Lagged_Kins_S,1))/5.086263020833333e+02,Lagged_Kins_S,'b','linewidth',1.5)
    hold on
    plot(velPredict_cv,'r-','linewidth',2)
    %plot((1:1:size(velPredict_cv,1))/5.086263020833333e+02,velPredict_cv,'r-','linewidth',2)
    set(gca,'FontSize',16)
%      xlim([0,19])
%      ylim([0,300])
     %xlabel('Time (s)','FontSize',16)
     %ylabel('Velocity (mm/s)','FontSize',16)
    %hold on; vline(1:1:19,'k--')
    legend('Recorded Kinematics','Predicted Kinematics','FontSize',14, 'location','northwest')
    %HighQualityFigs('Finger1_CV_Con')
    
    
    if Lag_Number~=0
        title(['Finger:',num2str(Finger),'; WLag:',num2str(Lag_Number),'; Predicted Output Using CV: R2: ',num2str(Rsq_Pre_output)])
    elseif Lag_Number==0
        title(['Finger:',num2str(Finger),'; Without Lag','; Predicted Output Using CV: R2: ',num2str(Rsq_Pre_output)])
        
    end
    
end
%  output: plots of k-fold results; b values; average of R2 values









