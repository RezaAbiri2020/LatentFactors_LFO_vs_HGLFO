
% for this commit
% changing the shuffling to phase of lfo
% adding circular shuffling 
% concatanating flexions and extensions for phase-amplitude coupling

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
    Fs_Phase=[];
    Fs_Amp_Direct=[];
    Fs_Amp_Avg=[];
    for k=1:2:length(Index_Fi)-1
        fprintf(['Flexion:',num2str(k),'\n'])
        Fs_Phase=[Fs_Phase; Phase(Index_Fi(k):Index_Fi(k+1),:)];
        Fs_Amp_Direct=[Fs_Amp_Direct; Amp_Direct(Index_Fi(k):Index_Fi(k+1),:)];
        Fs_Amp_Avg=[Fs_Amp_Avg; Amp_Avg(Index_Fi(k):Index_Fi(k+1),:)];
    end
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        %polarplot(Phase(1:end,ch),Amp_Direct(1:end,ch),'-')
        %polarplot(Phase(1:end,ch),Amp_Avg(1:end,ch),'-')
        FingerPAC(Fi).FsWholeHG.ValueCh(ch).Real=(abs(sum(Fs_Amp_Direct(1:end,ch).*exp(1i.*Fs_Phase(1:end,ch)))))/length(Fs_Phase);
        FingerPAC(Fi).FsAvgHG.ValueCh(ch).Real=(abs(sum(Fs_Amp_Avg(1:end,ch).*exp(1i.*Fs_Phase(1:end,ch)))))/length(Fs_Phase);
        
        PhaseCh=Fs_Phase(1:end,ch);
        for j=1:1000 % producing shuffled results
            Shuffled_PhaseCh= PhaseCh(randperm(numel(PhaseCh)));
            FingerPAC(Fi).FsWholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Fs_Amp_Direct(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Fs_Phase);
            FingerPAC(Fi).FsAvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Fs_Amp_Avg(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Fs_Phase);
            
        end
        % z-score per channel
        FingerPAC(Fi).FsWholeHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).FsWholeHG.ValueCh(ch).Real, FingerPAC(Fi).FsWholeHG.ValueCh(ch).Shuffled]);
        FingerPAC(Fi).FsAvgHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).FsAvgHG.ValueCh(ch).Real, FingerPAC(Fi).FsAvgHG.ValueCh(ch).Shuffled]);
        
    end
    
    % for extensions
    Es_Phase=[];
    Es_Amp_Direct=[];
    Es_Amp_Avg=[];
    for k=2:2:length(Index_Fi)-1
        fprintf(['Extension:',num2str(k),'\n'])
        Es_Phase=[Es_Phase; Phase(Index_Fi(k):Index_Fi(k+1),:)];
        Es_Amp_Direct=[Es_Amp_Direct; Amp_Direct(Index_Fi(k):Index_Fi(k+1),:)];
        Es_Amp_Avg=[Es_Amp_Avg; Amp_Avg(Index_Fi(k):Index_Fi(k+1),:)];
    end
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        %polarplot(Phase(1:end,ch),Amp_Direct(1:end,ch),'-')
        %polarplot(Phase(1:end,ch),Amp_Avg(1:end,ch),'-')
        FingerPAC(Fi).EsWholeHG.ValueCh(ch).Real=(abs(sum(Es_Amp_Direct(1:end,ch).*exp(1i.*Es_Phase(1:end,ch)))))/length(Es_Phase);
        FingerPAC(Fi).EsAvgHG.ValueCh(ch).Real=(abs(sum(Es_Amp_Avg(1:end,ch).*exp(1i.*Es_Phase(1:end,ch)))))/length(Es_Phase);
        
        PhaseCh=Es_Phase(1:end,ch);
        for j=1:1000 % producing shuffled results
            Shuffled_PhaseCh= PhaseCh(randperm(numel(PhaseCh)));
            FingerPAC(Fi).EsWholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Es_Amp_Direct(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Es_Phase);
            FingerPAC(Fi).EsAvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Es_Amp_Avg(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Es_Phase);
            
        end
        % z-score per channel
        FingerPAC(Fi).EsWholeHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).EsWholeHG.ValueCh(ch).Real, FingerPAC(Fi).EsWholeHG.ValueCh(ch).Shuffled]);
        FingerPAC(Fi).EsAvgHG.ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).EsAvgHG.ValueCh(ch).Real, FingerPAC(Fi).EsAvgHG.ValueCh(ch).Shuffled]);
        
    end   
end

% plot of observations; example
ch=1; %for Fi=?5
figure;
plot(Fs_Phase(:,ch),Fs_Amp_Direct(:,ch),'o')
figure;
plot(Fs_Phase(:,ch),Fs_Amp_Avg(:,ch),'o')
figure;
histogram2(Fs_Phase(:,ch),Fs_Amp_Direct(:,ch))
figure;
histogram2(Fs_Phase(:,ch),Fs_Amp_Avg(:,ch))

figure(1)
hist(FingerPAC(1).FsWholeHG.ValueCh(1).ZScAll,100)
set(gca,'fontsize',16)

figure(2)
hist(FingerPAC(1).FsAvgHG.ValueCh(1).ZScAll,100)
set(gca,'fontsize',16)

%Calculating p-value for all channels within Fs or Es
for Fi=1:5
    % for flexions
    for ch=1:256
        P_value1=length(find(FingerPAC(Fi).FsWholeHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).FsWholeHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).FsWholeHG.ValueCh(ch).ZScAll);
        
        P_value2=length(find(FingerPAC(Fi).FsAvgHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).FsAvgHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).FsAvgHG.ValueCh(ch).ZScAll);
        
        FingerPAC(Fi).FsWholeHG.ValueCh(ch).PValue=P_value1;
        FingerPAC(Fi).FsAvgHG.ValueCh(ch).PValue=P_value2;
    end
    
    % for extensions
    for ch=1:256
        P_value1=length(find(FingerPAC(Fi).EsWholeHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).EsWholeHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).EsWholeHG.ValueCh(ch).ZScAll);
        
        P_value2=length(find(FingerPAC(Fi).EsAvgHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).EsAvgHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).EsAvgHG.ValueCh(ch).ZScAll);
        
        FingerPAC(Fi).EsWholeHG.ValueCh(ch).PValue=P_value1;
        FingerPAC(Fi).EsAvgHG.ValueCh(ch).PValue=P_value2;
    end
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingerPAC.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger/flexions or extensions with hline in 0.05
Fi=1;
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).FsWholeHG.ValueCh(ch).PValue];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).FsAvgHG.ValueCh(ch).PValue];
end

figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Direct-LFO')

subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Avg-LFO')

% showing significant channels on brain map per finger

load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')
% for flexion
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).FsWholeHG.ValueCh(ch).PValue];
    end
    PValuesDirect(find(PValuesDirect<0.05))=-1; %significant
    PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
      
end
HighQualityFigs('Fs_Direct')

for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).FsAvgHG.ValueCh(ch).PValue];
    end
    PValuesAvg(find(PValuesAvg<0.05))=-1; %significant
    PValuesAvg(find(PValuesAvg>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
      
end
% for extension
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).EsWholeHG.ValueCh(ch).PValue];
    end
    PValuesDirect(find(PValuesDirect<0.05))=-1; %significant
    PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
      
end
HighQualityFigs('Es_Direct')

for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).EsAvgHG.ValueCh(ch).PValue];
    end
    PValuesAvg(find(PValuesAvg<0.05))=-1; %significant
    PValuesAvg(find(PValuesAvg>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesAvg),'rh'); % rho is a 256ch vector
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
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                
                FingerPPC(Fi).Flexion(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);                
                FingerPPC(Fi).Flexion(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);                
                
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
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                
                FingerPPC(Fi).Extension(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);
                FingerPPC(Fi).Extension(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);
                
            end
            
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

%%
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


