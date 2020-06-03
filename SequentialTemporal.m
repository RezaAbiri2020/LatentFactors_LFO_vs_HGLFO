
%% load the ERPs
load('DeltaPure_ERPHilbert.mat')
load('HG_ERPHilbert.mat')

%% computation
for Fi=5
    
    %data1=HG_ERPHilbert(Fi).FingerF2E;
    data1=DeltaPure_ERPHilbert(Fi).FingerF2E;
    
    data2=zeros(size(data1));
    for i=1:256
        data2(:,i)=normalize(data1(:,i),'range');
    end
    
    data3=data2';
    figure;
    set(gcf, 'Position', [300, 100, 700, 850]);
    %suptitle(['All Ch; HG-LFO for Finger: ',num2str(Fi)]);
    suptitle(['All Ch; LFO for Finger: ',num2str(Fi)]);
    imagesc(data3)
    colorbar
    
    Rows=[];
    for i=1:10
        [m,n]=max(data3(:,(i-1)*50+1:i*50));
        
        for j=1:length(n)
            if isempty(find(Rows==n(j)))
                Rows=[Rows;n(j)];
            end
        end
        
    end
    
    OrderedData1=data3(Rows,:);
    figure;
    set(gcf, 'Position', [300, 100, 700, 850]);
    %suptitle(['Ordered Ch; HG-LFO for Finger: ',num2str(Fi)]);
    suptitle(['Ordered Ch; LFO for Finger: ',num2str(Fi)]);
    imagesc(OrderedData1)
    colorbar
    yticks('')
    yticks(1:3:length(Rows))
    ytickangle(0)
    yticklabels(num2str(Rows(1:3:end)))
    
    % modification for ordereddata
    OrderedData1=zeros(length(Rows)*3,256);
    for i=1:length(Rows)
        OrderedData1((i-1)*3+1:i*3,Rows(i))=[1;1;1];   
    end
    
    %HG_E2F_Ordered(Fi).FingerData=OrderedData1;
    %HG_E2F_Ordered(Fi).RowsData=Rows;
    
    %Delta_F2E_Ordered(Fi).FingerData=OrderedData1;
    %Delta_F2E_Ordered(Fi).RowsData=Rows;
end

%save('HG_E2F_Ordered.mat','HG_E2F_Ordered')
%save('Delta_F2E_Ordered.mat','Delta_F2E_Ordered')
%close all

%% map delta channel to HG

data1_new=HG_ERPHilbert(Fi).FingerF2E;
data2_new=zeros(size(data1));
    for i=1:256
        data2_new(:,i)=normalize(data1_new(:,i),'range');
    end
    
data2_new=data2_new';

    
OrderedData1_new=data2_new(Rows,:);
figure;
set(gcf, 'Position', [300, 100, 700, 850]);
%suptitle(['Ordered Ch; HG-LFO for Finger: ',num2str(Fi)]);
suptitle(['Ordered Ch of HG-LFO; projection of LFO Ch for Finger: ',num2str(Fi)]);
imagesc(OrderedData1_new)
colorbar
yticks('')
yticks(1:3:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:3:end)))
