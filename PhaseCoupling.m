
% for this commit
% 


clear all
close all
clc
% observing different coupling analysis for these signals:
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')

% using indexing for flexsion/extension
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\FingersKinIndexes.mat')

%% phase-amplitude coupling for (whole movement)

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        %polarplot(Phase(1:end,ch),Amp(1:end,ch),'-')
        PAC_Value1=((sum(Amp_Direct(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        PAC_Value2=((sum(Amp_Avg(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex = PAC_Value1;
        FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex = PAC_Value2;
        
        PhaseCh=Phase(1:end,ch);
        for j=1:2000 % producing shuffled results
            Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
            FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Direct(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase);
            FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Avg(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase);
            
        end
        % z-score per channel
        FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll=...
            zscore([abs(PAC_Value1), FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled]);
        FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll=...
            zscore([abs(PAC_Value2), FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled]);
    end
end

% plot of observations; example
figure
hist(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll,100)
set(gca,'fontsize',16)
figure
hist(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll,100)
set(gca,'fontsize',16)

%Calculating p-value for all channels
for Fi=1:5
    for ch=1:256
        Pvalue1=length(find(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll);
        FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll);
        FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
end

save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC_PA_WholeMove.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
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

HighQualityFigs('WholeMovement_PAC_PA')

% plot of example of significant chanel
Fi=5; ch=27;
Phase=angle(hilbert(LFO_signals(Fi).Filtered));
Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
figure;
histogram2(Phase(:,ch),Amp_Direct(:,ch),150)
xlabel('Phase LFO')
ylabel('Amplitude HG')
zlabel('Number of repetition')
HighQualityFigs('hist2D_PAC_PA')

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    end
    SigChs=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PA_Grid_Direct')

% for avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
    end
    SigChs=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PA_Grid_Avg')

%% phase-amplitude coupling for (whole movement)
% how the phase of LFO is coupled to the phase of
% HG-Direct-LFO & HG-Avg-LFO during whole movement of one finger
for Fi=1%:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Phase_Direct=angle(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Phase_Avg=angle(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        
        PAC_Value1=(sum(exp(1i.*(Phase(1:end,ch)-Phase_Direct(1:end,ch)))))/length(Phase);
        PAC_Value2=(sum(exp(1i.*(Phase_Avg(1:end,ch)-Phase(1:end,ch)))))/length(Phase);
        
        FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex=PAC_Value1;
        FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex=PAC_Value2;
        
        PhaseCh=Phase(1:end,ch);
        PhaseCh2=Phase_Direct(1:end,ch);
        for j=1:1000 % producing shuffled results
            Shuffled_PhaseCh = PhaseCh(randperm(numel(PhaseCh)));
            %Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
            Shuffled_Phase_Direct = PhaseCh2(randperm(numel(PhaseCh2)));
            %Shuffled_Phase_Direct = circshift(Phase_Direct(1:end,ch),randi(length(Phase_Direct(1:end,ch))));
            FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled(j)=...
                (sum(exp(1i.*(Phase(1:end,ch)-Shuffled_Phase_Direct(1:end)))))/length(Phase);
            FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled(j)=...
                (sum(exp(1i.*(Phase_Avg(1:end,ch)-Shuffled_PhaseCh(1:end)))))/length(Phase);
        end
    end
    
end

%Calculating p-value for all channels
for Fi=1%:5
    for ch=1:256
        RealValue1 = FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex;
        ShuffledValue1 = FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled;
        
        Pvalue1=length(find(abs((RealValue1))>abs((ShuffledValue1))))/...
            (length(FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)=Pvalue1;
        
        RealValue2 = FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex;
        ShuffledValue2 = FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled;
        
        Pvalue2=length(find(abs((RealValue2))>abs((ShuffledValue2))))/...
            (length(FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)=Pvalue2;
            
    end
end

save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC_PP_WholeMove.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
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

HighQualityFigs('WholeMovement_PAC_PP')

% plot of example of significant chanel
Fi=5; ch=27;
Phase=angle(hilbert(LFO_signals(Fi).Filtered));
Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
figure;
histogram2(Phase(:,ch),Amp_Direct(:,ch),150)
xlabel('Phase LFO')
ylabel('Amplitude HG')
zlabel('Number of repetition')
HighQualityFigs('hist2D_PAC_PP')

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    end
    SigChs=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PP_Grid_Direct')

% for avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
    end
    SigChs=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PP_Grid_Avg')

%% phase-amplitude coupling for considering a flexion and a extension as one trial
% how the phase of LFO is coupled to the power(abs-hilbert) of HG during flexions or extensions
% it should be the envelope of HG-Direct-LFO or HG-Avg-LFO 

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions and extension
    kk=0;
    for k=1:2:length(Index_Fi)-2
        fprintf(['Flexion & Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+2),:);
        Amp_Direct_k=Amp_Direct(Index_Fi(k):Index_Fi(k+2),:);
        Amp_Avg_k=Amp_Avg(Index_Fi(k):Index_Fi(k+2),:);
        kk=kk+1;
                
        for ch=1:256
            %fprintf(['Ch:',num2str(ch),'\n'])
            %polarplot(Phase(1:end,ch),Amp_Direct(1:end,ch),'-')
            %polarplot(Phase(1:end,ch),Amp_Avg(1:end,ch),'-')
            PAC_Value1=((sum(Amp_Direct_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            PAC_Value2=((sum(Amp_Avg_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).RealComplex = PAC_Value1;
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).RealComplex = PAC_Value2;
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Direct_k(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase_k);
                FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Avg_k(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase_k);
                
            end
            % z-score per channel
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).ZScAll=...
                zscore([abs(PAC_Value1), FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled]);
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).ZScAll=...
                zscore([abs(PAC_Value2), FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled]);     
        end
    end   
end

% plot of observations; example
figure;
histogram2(Phase_k(:,ch),Amp_Direct_k(:,ch))
figure;
histogram2(Phase_k(:,ch),Amp_Avg_k(:,ch))

%Calculating p-value for all channels across trials
for Fi=1:5
    % for flexion and extension
    for ch=1:256
        RealValue1=[];
        ShuffledValue1=[];
        RealValue2=[];
        ShuffledValue2=[];
        for trial=1:length(FingerPAC(Fi).FE)
            RealValue1=[RealValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled];
        end
        Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
            (length(FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
            (length(FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
     
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingerPAC.mat','FingerPAC','-v7.3')
% stem plots of p-values for all channels per finger with hline in 0.05

Fi=5;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Avg-LFO')

HighQualityFigs('FE_PAC')


% showing significant channels on brain map per finger

load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for flexion-extension; for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    end
    SigChs=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('FE_PAC_Grid_Direct')

% for flexion-extension; for Avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
    end
    SigChs=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('FE_PAC_Grid_Avg')


%% phase-amplitude coupling for considering a flexion and a extension as one trial
% how the phase of LFO is coupled to the phase of
% HG-Direct-LFO & HG-Avg-LFO during flexions and extensions

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Phase_Direct=angle(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Phase_Avg=angle(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions and extension
    kk=0;
    for k=1:2:length(Index_Fi)-2
        fprintf(['Flexion & Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+2),:);
        Phase_Direct_k=Phase_Direct(Index_Fi(k):Index_Fi(k+2),:);
        Phase_Avg_k=Phase_Avg(Index_Fi(k):Index_Fi(k+2),:);
        kk=kk+1;
        
        for ch=1:256
            PAC_Value1=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);
            PAC_Value2=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);
            
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).RealComplex=PAC_Value1;
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).RealComplex=PAC_Value2;
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = PhaseCh(randperm(numel(PhaseCh)));
                %Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                
                FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);                
                FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);                
                
            end
           
        end
    end
          
end

%Calculating p-value for all channels across trials and fingers

for ch=1:256
    fprintf(['Ch:',num2str(ch),'\n'])
    RealValue1=[];
    ShuffledValue1=[];
    RealValue2=[];
    ShuffledValue2=[];
    for Fi=1:5
        % for flexion and extension
        for trial=1:length(FingerPAC(Fi).FE)
            RealValue1=[RealValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled];
        end
    end
    Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
        (length(FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled)+1);
    FingerPAC(1).PValuesFE.WholeHG.ValueCh(ch)=Pvalue1;
    
    Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
        (length(FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled)+1);
    FingerPAC(1).PValuesFE.AvgHG.ValueCh(ch)=Pvalue2;
    
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Avg-LFO')

HighQualityFigs('FE_PAC')


if 0
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
        if FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).PValue<0.05
            polarplot([0 FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-r','linewidth',1.5)
            rticks('')
            thetaticks('')
        else
            polarplot([0 FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-b','linewidth',1.5)
            rticks('')
            thetaticks('')
        end
        
        hold on
          
    end
end
end

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for flexion-extension; for direct HG
for Fi=1
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    end
    SigChs=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    figure;
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end

% for flexion-extension; for Avg HG
for Fi=1
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
    end
    SigChs=find(PValuesAvg<0.05); %significant
    %PValuesAvg(find(PValuesAvg>=0.05))=0; % not significant
    figure;
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar  
end


