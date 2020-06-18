
clear all
close all
clc
% observing different coupling analysis for these signals:
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')

% using indexing for flexsion/extension
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\FingersKinIndexes.mat')

%% phase-amplitude coupling
% how the phase of LFO is coupled to the power(abs-hilbert) of
% HG-Direct-LFO & HG-Avg-LFO during flexions or extensions

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions
    kk=0;
    for k=1:2:length(Index_Fi)-1
        fprintf(['Flexion:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+1),:);
        Amp_Direct_k=Amp_Direct(Index_Fi(k):Index_Fi(k+1),:);
        Amp_Avg_k=Amp_Avg(Index_Fi(k):Index_Fi(k+1),:);
        kk=kk+1;
        
        for ch=1:256
            %polarplot(Phase_k(1:end,ch),Amp_Direct_k(1:end,ch),'-')
            %polarplot(Phase_k(1:end,ch),Amp_Avg_k(1:end,ch),'-')
            FingerPAC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Real=(abs(sum(Amp_Direct_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            FingerPAC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Real=(abs(sum(Amp_Avg_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            
            AmpCh1=Amp_Direct_k(1:end,ch);
            AmpCh2=Amp_Avg_k(1:end,ch);
            for j=1:2000 % producing shuffled results
                Shuffled_AmpCh1 = AmpCh1(randperm(numel(AmpCh1)));
                FingerPAC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Shuffled_AmpCh1(1:end).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
                Shuffled_AmpCh2 = AmpCh2(randperm(numel(AmpCh2)));
                FingerPAC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Shuffled_AmpCh2(1:end).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
                
            end
            % z-score per channel
            FingerPAC(Fi).Flexion(kk).WholeHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Real, FingerPAC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Shuffled]);
            FingerPAC(Fi).Flexion(kk).AvgHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Real, FingerPAC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Shuffled]);
            
        end
    end
    
    % for extensions
    kk=0;
    for k=2:2:length(Index_Fi)-1
        fprintf(['Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+1),:);
        Amp_Direct_k=Amp_Direct(Index_Fi(k):Index_Fi(k+1),:);
        Amp_Avg_k=Amp_Avg(Index_Fi(k):Index_Fi(k+1),:);
        kk=kk+1;
        
        for ch=1:256
            %polarplot(Phase_k(1:end,ch),Amp_Direct_k(1:end,ch),'-')
            %polarplot(Phase_k(1:end,ch),Amp_Avg_k(1:end,ch),'-')
            FingerPAC(Fi).Extension(kk).WholeHG.ValueCh(ch).Real=(abs(sum(Amp_Direct_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            FingerPAC(Fi).Extension(kk).AvgHG.ValueCh(ch).Real=(abs(sum(Amp_Avg_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            
            AmpCh1=Amp_Direct_k(1:end,ch);
            AmpCh2=Amp_Avg_k(1:end,ch);
            for j=1:2000 % producing shuffled results
                Shuffled_AmpCh1 = AmpCh1(randperm(numel(AmpCh1)));
                FingerPAC(Fi).Extension(kk).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Shuffled_AmpCh1(1:end).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
                Shuffled_AmpCh2 = AmpCh2(randperm(numel(AmpCh2)));
                FingerPAC(Fi).Extension(kk).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Shuffled_AmpCh2(1:end).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
                
            end
            % z-score per channel
            FingerPAC(Fi).Extension(kk).WholeHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).Extension(kk).WholeHG.ValueCh(ch).Real, FingerPAC(Fi).Extension(kk).WholeHG.ValueCh(ch).Shuffled]);
            FingerPAC(Fi).Extension(kk).AvgHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).Extension(kk).AvgHG.ValueCh(ch).Real, FingerPAC(Fi).Extension(kk).AvgHG.ValueCh(ch).Shuffled]);
            
        end
    end
    
end

% plot of observations; example
figure(1)
hist(FingerPAC(1).Flexion(1).WholeHG.ValueCh(1).ZScAll,100)
set(gca,'fontsize',16)

figure(2)
hist(FingerPAC(1).Flexion(1).AvgHG.ValueCh(1).ZScAll,100)
set(gca,'fontsize',16)

%Calculating p-value for all channels within each trial
for Fi=1:5
    % for flexions
    for flexion=1:length(FingerPAC(Fi).Flexion)
        for ch=1:256
            P_value1=length(find(FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).ZScAll(2:end)>...
                FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).ZScAll(1)))...
                /length(FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).ZScAll);
            
            P_value2=length(find(FingerPAC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).ZScAll(2:end)>...
                FingerPAC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).ZScAll(1)))...
                /length(FingerPAC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).ZScAll);
            
            FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).PValue=P_value1;
            FingerPAC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).PValue=P_value2;
        end
    end
    % for extensions
    for extension=1:length(FingerPAC(Fi).Extension)
        for ch=1:256
            P_value1=length(find(FingerPAC(Fi).Extension(extension).WholeHG.ValueCh(ch).ZScAll(2:end)>...
                FingerPAC(Fi).Extension(extension).WholeHG.ValueCh(ch).ZScAll(1)))...
                /length(FingerPAC(Fi).Extension(extension).WholeHG.ValueCh(ch).ZScAll);
            
            P_value2=length(find(FingerPAC(Fi).Extension(extension).AvgHG.ValueCh(ch).ZScAll(2:end)>...
                FingerPAC(Fi).Extension(extension).AvgHG.ValueCh(ch).ZScAll(1)))...
                /length(FingerPAC(Fi).Extension(extension).AvgHG.ValueCh(ch).ZScAll);
            
            FingerPAC(Fi).Extension(extension).WholeHG.ValueCh(ch).PValue=P_value1;
            FingerPAC(Fi).Extension(extension).AvgHG.ValueCh(ch).PValue=P_value2;
        end
          
    end
    
end
%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingerPAC.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger/flexions or extensions with hline in 0.05
Fi=2;
for flexion=1:length(FingerPAC(Fi).Flexion)
    PValuesDirect=[];
    PValuesAvg=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).PValue];
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).PValue];
    end
    FlexionPValuesDirect(:,flexion)=PValuesDirect;
    FlexionPValuesAvg(:,flexion)=PValuesAvg;
%     figure(1);
%     subplot(5,1,flexion)
%     stem(PValuesDirect)
%     xlim([0,256])
%     hline(0.05,'r')
%     set(gca,'fontsize',16)
%     
%     figure(2);
%     subplot(5,1,flexion)
%     stem(PValuesAvg)
%     xlim([0,256])
%     hline(0.05,'r')
%     set(gca,'fontsize',16)
    
end

SigChs=zeros(256,1);
for ch=1:256
    if (FlexionPValuesDirect(ch,:)<0.05)
      SigChs(ch,1)=1;  
    end 
end

% showing significant channels on brain map per finger

load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

for Fi=1:5
    flexion=1;
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).PValue];
    end
    PValuesDirect(find(PValuesDirect<0.05))=-1; %significant
    PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
    
    
end

%% phase-phase coupling
% how the phase of LFO is coupled to the phase of
% HG-Direct-LFO & HG-Avg-LFO during flexions or extensions

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Phase_Direct=angle(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Phase_Avg=angle(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions
    kk=0;
    for k=1:2:length(Index_Fi)-1
        fprintf(['Flexion:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+1),:);
        Phase_Direct_k=Phase_Direct(Index_Fi(k):Index_Fi(k+1),:);
        Phase_Avg_k=Phase_Avg(Index_Fi(k):Index_Fi(k+1),:);
        kk=kk+1;
        
        for ch=1:256
            PPC_Value1=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);
            PPC_Value2=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);
            
            FingerPPC(Fi).Flexion(kk).WholeHG.ValueCh(ch).RealComplex=PPC_Value1;
            FingerPPC(Fi).Flexion(kk).AvgHG.ValueCh(ch).RealComplex=PPC_Value2;
            
            Phase1Ch=Phase_Direct_k(1:end,ch);
            Phase2Ch=Phase_Avg_k(1:end,ch);
            for j=1:2000 % producing shuffled results
                Shuffled_Phase1Ch = Phase1Ch(randperm(numel(Phase1Ch)));
                Shuffled_Phase2Ch = Phase2Ch(randperm(numel(Phase2Ch)));
                
                FingerPPC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Phase_k(1:end,ch)-Shuffled_Phase1Ch(1:end)))))/length(Phase_k);                
                FingerPPC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Phase_k(1:end,ch)-Shuffled_Phase2Ch(1:end)))))/length(Phase_k);                
                
            end
           
        end
    end
    
    % for extensions
    kk=0;
    for k=2:2:length(Index_Fi)-1
        fprintf(['Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+1),:);
        Phase_Direct_k=Phase_Direct(Index_Fi(k):Index_Fi(k+1),:);
        Phase_Avg_k=Phase_Avg(Index_Fi(k):Index_Fi(k+1),:);
        kk=kk+1;
        
        for ch=1:256
            PPC_Value1=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);
            PPC_Value2=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);
            
            FingerPPC(Fi).Extension(kk).WholeHG.ValueCh(ch).RealComplex=PPC_Value1;
            FingerPPC(Fi).Extension(kk).AvgHG.ValueCh(ch).RealComplex=PPC_Value2;
            
            Phase1Ch=Phase_Direct_k(1:end,ch);
            Phase2Ch=Phase_Avg_k(1:end,ch);
            for j=1:2000 % producing shuffled results
                Shuffled_Phase1Ch = Phase1Ch(randperm(numel(Phase1Ch)));
                Shuffled_Phase2Ch = Phase2Ch(randperm(numel(Phase2Ch)));
                
                FingerPPC(Fi).Extension(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Phase_k(1:end,ch)-Shuffled_Phase1Ch(1:end)))))/length(Phase_k);
                FingerPPC(Fi).Extension(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Phase_k(1:end,ch)-Shuffled_Phase2Ch(1:end)))))/length(Phase_k);
                
            end
            
        end
    end
      
end

%Calculating p-value for all channels within each trial
for Fi=1:5
    % for flexions
    for flexion=1:length(FingerPPC(Fi).Flexion)
        for ch=1:256
            P_value1=length(find(abs(FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).RealComplex)>...
                abs(FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).Shuffled)))/...
                (length(FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).Shuffled)+1);
            FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).PValue=P_value1;
            
            P_value2=length(find(abs(FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).RealComplex)>...
                abs(FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).Shuffled)))/...
                (length(FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).Shuffled)+1);
            FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).PValue=P_value2;
            
        end
    end
    % for extensions
    for extension=1:length(FingerPPC(Fi).Extension)
        for ch=1:256
            P_value1=length(find(abs(FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).RealComplex)>...
                abs(FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).Shuffled)))/...
                (length(FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).Shuffled)+1);
            FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).PValue=P_value1;
            
            P_value2=length(find(abs(FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).RealComplex)>...
                abs(FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).Shuffled)))/...
                (length(FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).Shuffled)+1);
            FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).PValue=P_value2;
            
        end      
    end   
end


%Calculating p-value for all channels across trials
for Fi=1:5
    % for flexion
    for ch=1:256
        RealValue1=[];
        ShuffledValue1=[];
        RealValue2=[];
        ShuffledValue2=[];
        for flexion=1:length(FingerPPC(Fi).Flexion)
            RealValue1=[RealValue1; FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).Shuffled];
        end
        Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
            (length(FingerPPC(Fi).Flexion(flexion).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPPC(Fi).PValuesFlexion.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
            (length(FingerPPC(Fi).Flexion(flexion).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPPC(Fi).PValuesFlexion.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
    % for extension
    for ch=1:256
        RealValue1=[];
        ShuffledValue1=[];
        RealValue2=[];
        ShuffledValue2=[];
        for extension=1:length(FingerPPC(Fi).Extension)
            RealValue1=[RealValue1; FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).Shuffled];
        end
        Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
            (length(FingerPPC(Fi).Extension(extension).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPPC(Fi).PValuesExtension.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
            (length(FingerPPC(Fi).Extension(extension).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPPC(Fi).PValuesExtension.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
    
end
%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingerPPC.mat','FingerPPC','-v7.3')


% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPPC(Fi).PValuesFlexion.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPPC(Fi).PValuesFlexion.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion; HG-Avg-LFO')

% for extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPPC(Fi).PValuesExtension.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPPC(Fi).PValuesExtension.AvgHG.ValueCh(ch)];
end
figure(2);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Extension; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Extension; HG-Avg-LFO')

%%
% compass plots of coupled-phase per channel for all finger with 
% significant channel with red lines
Nch_Record=256;

% grid layout
Ch_num_1=1:Nch_Record;
Ch_num_2=reshape(Ch_num_1,[16,16]);
% Ch_num_3 is the grid layout
Ch_num_3=rot90(rot90(Ch_num_2));

for Fi=1:5
    figure;
    set(gcf, 'Position', [100, 100, 1000, 800]);
    suptitle(['coupled-angle per channel; Finger: ',num2str(Fi)]);
    
    for i=1:256
        subplot(16,16,i)
        r=ceil(i/16);
        c=mod(i,16);
        if c==0
            c=8;
        end
        if FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).PValue<0.05
            polarplot([0 FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-r','linewidth',1.5)
            rticks('')
            thetaticks('')
        else
            polarplot([0 FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-b','linewidth',1.5)
            rticks('')
            thetaticks('')
        end
        
        hold on
          
    end
end

% showing significant channels on brain map per finger

load('..\ImagingData\EC171_rh_pial.mat')
load('..\ImagingData\TDT_elecs_all.mat')

for Fi=1:5
    PValuesFi=[];
    for ch=1:256
        PValuesFi=[PValuesFi; FingerPPC(Fi).ValueCh(ch).PValue];
    end
    PValuesFi(find(PValuesFi<0.05))=-1; %significant
    PValuesFi(find(PValuesFi>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesFi),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
    
    
end

save('FingerPPC.mat','FingerPPC')


%% connection of significant channels to Weights of vel/pos prediction
load('HG_AllWeights.mat');



%% rotation dynamics for high gamma in Ch 186, 187
% State space model estimation:
% Xdot=A*X
Samples=amp(:,[49 50]);
% magnify the values
X=1*Samples;
Xdot=diff(X);
X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 
% solving for phase plane
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[4, 0]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);

% rotation dynamics for delta in Ch 186, 187
% State space model estimation:
% Xdot=A*X
Samples=Delta_Fingers(Fi).Hilbert(:,[49 50]);
% magnify the values
X=1*Samples;
Xdot=diff(X);
X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 
% solving for phase plane
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[4, 0]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);


