% Microglia vs Neuron Interactions/Co-localization Analysis
% Ali Ozgur Argunsah and Lorenzo Gesuita, 2017, Zurich.
% Input: 2 channel RGB z-stack

clear all
close all
clc

addpath(genpath('C:\Users\argunsah\Documents\lorenzoPaperCodes\colocCodes'));

dataFolder  = 'W:\Assembly\Lorenzo\Pics_Confocal\Confocal_Microglia Project Images\Microglia_Time Course for Contacts\P4-P11_SSTCreAi14_Iba1488_Layer5_40x_2048x2048';
D1          = uipickfiles('FilterSpec',dataFolder);

load('ECOC_classifier_Lorenzo_1600x1600_40X.mat');

saveFolder  = 'E:\';

errorvector = zeros(1,size(D1,2))
for i = 1:size(D1,2)
    try
    i
    clear info;
    clear multiCube;
    
    [info, fPath, fName, fExt, multiCube]   = readMicroscopyData(D1{i},3,1,0,0);
    
    close all;
    % Create Info Structure
    clear data;

    regexpInd       = regexp(D1{i},'\');
    %saveFolder      = 'W:\Assembly\Lorenzo\Pics_Confocal\Confocal_Microglia Project Images\Microglia_Time Course for Contacts\Results_050319_L5';
    saveDataName    = fName;

    info.name           = fName;
    info.analysisday    = date;
    info.Analyzer       = 'Automatic';
    info.DataCollected  = 'Martina';
    info.saveDataName   = saveDataName;

    clear I1norm;
    clear I2norm;
    
    I1norm  = squeeze(multiCube(:,:,:,1))/4095;
    I2norm  = squeeze(multiCube(:,:,:,2))/4095;

    cube1 = I1norm;
    cube2 = I2norm;

    info.cube1 = I1norm;
    info.cube2 = I2norm;
    
    % Cell Body and Dendrite Segmentation

    mip2        = max(cube2,[],3);
    mip         = bilateralFilter(imadjust(mip2));
    V           = vesselness2D(mip, .05:0.05:0.25, [.1;.1], .95, true);
    thr         = graythresh(V);
    bwV         = V>thr;
    bwVfill     = regionfill(mip2,bwV);
    bwVfillBWF  = bilateralFilter(bwVfill);

    numIter     = 15;
    gradThresh  = multithresh(bwVfillBWF,numIter);
    C           = imdiffusefilt(bwVfillBWF,'ConductionMethod','quadratic', ...
        'GradientThreshold',gradThresh,'NumberOfIterations',numIter);

    bwCandidates    = C > graythresh(C);
    rpC             = regionprops(bwCandidates,mip2,'All');
    cents           = cat(1, rpC.Centroid);
    Areas2D         = cat(1, rpC.Area);
    meanIntensity2D = cat(1, rpC.MeanIntensity);
    ecc2D           = cat(1, rpC.Eccentricity);
    Extent2D        = cat(1, rpC.Extent);
    EulerNumber2D   = cat(1, rpC.EulerNumber);
    EquivDiameter2D = cat(1, rpC.EquivDiameter);
    
    combinedData = [Areas2D,...
    meanIntensity2D,...
    ecc2D,...
    Extent2D,...
    EulerNumber2D,...
    EquivDiameter2D];

    label             = predict(Mdl,combinedData);
    data.rpC          = rpC;
    data.label        = label;
    
    CC                = bwconncomp(bwCandidates);
    bwCandidatesNew1  = bwCandidates;
    bwCandidatesNew2  = zeros(size(bwCandidates));
    bwCandidatesNew3  = bwCandidatesNew2;
    
    numObjects = size(Areas2D,1);
    
    bwtempAll = zeros(size(bwCandidates));
    [L, num]  = bwlabel(bwCandidates);
    for obj = 1 : num
        if label(obj)>1
            thisBlob = ismember(L, obj);
            if label(obj) == 2
                acVal   = floor(40*(1-meanIntensity2D(obj)))+1;
                try
                    thisBlob = activecontour(mip2 ,thisBlob,acVal,'chan-vese','ContractionBias',-.1,'SmoothFactor',1);
                catch
                    thisBlob = ismember(L, obj);
                end
            end
            if label(obj) == 3
                
                if meanIntensity2D(obj) > 0.75
                    acVal   = floor(20*(1-meanIntensity2D(obj)))+1;
                    try
                        thisBlob = activecontour(mip2 ,thisBlob,acVal,'edge','ContractionBias',.1,'SmoothFactor',1); 
                    catch
                    thisBlob = ismember(L, obj);
                    end
                else
                    acVal   = floor(50*(1-meanIntensity2D(obj)))+1;
                    try
                        thisBlob = activecontour(mip2 ,thisBlob,acVal,'edge','ContractionBias',-.1,'SmoothFactor',1);
                    catch
                    thisBlob = ismember(L, obj);
                    end
                end
            end
            bwtempAll = bwtempAll | thisBlob;
        end
    end

    bwCandidatesNew      = bwtempAll;
    cellBodiesBW         = bwCandidatesNew;   
    cellBodiesBWbig      = bwareaopen(cellBodiesBW,20000);
    cellBodiesBW         = cellBodiesBW - cellBodiesBWbig;
    bwCandidatesNewThick = bwmorph(cellBodiesBW,'thicken',5);

    [bw1,~,~]        = filterWrongNeurons(bwCandidatesNewThick,mip2,Mdl);
    bwCandidatesNewThick = bwCandidatesNewThick - bw1;
    
    [bw1,~,~]        = filterWrongNeurons(cellBodiesBW,mip2,Mdl);
    cellBodiesBW         = cellBodiesBW - bw1;
    
    cellBodiesBW         = bwmorph(cellBodiesBW,'thicken',3);
    dendrites            = (bwV - cellBodiesBW) > 0;
   
    info.cellBodiesBW    = cellBodiesBW;
    info.dendritesBW     = dendrites;   
      
    [lb,center] = adaptcluster_kmeans(max(cube1,[],3));
    clustNum    = max(lb(:));
    
    microglia = lb>1;
    info.microgliaBW    = microglia;
    mip1                = max(cube1,[],3);

    rgb(:,:,1) = double(cellBodiesBW)*255;
    rgb(:,:,2) = double(dendrites)*255;
    rgb(:,:,3) = double(microglia)*255;
    
    rgbRaw(:,:,1) = double(bwperim(cellBodiesBW));
    rgbRaw(:,:,2) = double(mip2);
    rgbRaw(:,:,3) = double(mip1);
    
    figure, imshow(rgbRaw,[]);  
    saveas(gcf,fullfile(saveFolder,sprintf('rgbRawl_%s%s',saveDataName,'.eps')),'epsc');
    saveas(gcf,fullfile(saveFolder,sprintf('rgbRaw_%s%s',saveDataName,'.png')),'png');
     
    figure, imshow(rgb,[]);  
    saveas(gcf,fullfile(saveFolder,sprintf('rgb_%s%s',saveDataName,'.eps')),'epsc');
    saveas(gcf,fullfile(saveFolder,sprintf('rgb_%s%s',saveDataName,'.png')),'png');
     
    cellBodies3D = repmat(double(cellBodiesBW),1,1,size(cube2,3));
    dendrites3D  = repmat(double(dendrites),1,1,size(cube2,3));
    microglia3D  = repmat(double(microglia),1,1,size(cube2,3));

    clear bw_MG;
    clear bw_D;
    clear bw_CB;
    
    bw_MG     = thresholdWithZprofile(cube1,microglia3D);
    bw_D      = thresholdWithZprofile(cube2,dendrites3D);
    bw_CB     = thresholdWithZprofile(cube2,cellBodies3D);
    
    cellBodies3D = cellBodies3D & bw_CB;
    dendrites3D  = dendrites3D  & bw_D;
    microglia3D  = microglia3D  & bw_MG;
    
    info.cellBodies3DBW = cellBodies3D;
    info.dendrites3DBW  = dendrites3D;
    info.microglia3DBW  = microglia3D;
    
    save(fullfile(saveFolder,sprintf('%s_info.mat',saveDataName)), '-struct', 'info');
    catch
        errorvector(i) = 1;
    end
end
%% Analyze Results Here

D1 = uipickfiles('FilterSpec',dataFolder);

load('ECOC_classifier_Lorenzo_1600x1600_40X.mat');

for i = 1:size(D1,2)
    i
    load(D1{i});
    ageInd1     = regexp(D1{i},'\');
    ageInd2     = regexp(D1{i},'_');
    ageList(i)  = str2double(D1{i}(ageInd1(5)+2:ageInd2(2)-1));
    
    mipTemp     = max(cube2,[],3);
    binaryImage = max(cellBodies3DBW,[],3);
    binaryImage = filterFunkyShapes(binaryImage);
    
    cellBodies3DBWNew = repmat(double(binaryImage),1,1,size(cube2,3));
    cellBodies3DBWNew = cellBodies3DBWNew & cellBodies3DBW;

    BW_MG_BODY        = cellBodies3DBWNew & microglia3DBW;
    BW_MG_DEND        = dendrites3DBW     & microglia3DBW;

    clear CC_cellBodies;
    clear L;
    
    CC_cellBodies = bwconncomp(bwmorph(binaryImage,'open'));
    numNeurons(i) = CC_cellBodies.NumObjects;
    
    bwtempAll = zeros(size(binaryImage));
    [L, num]  = bwlabel(bwmorph(binaryImage,'open'));
    for obj = 1 : num
        thisBlob = ismember(L, obj);
        thisBlob3D = repmat(double(thisBlob),1,1,size(cube2,3)) & cellBodies3DBWNew;
    
        thisBlob3D_BW_MG_BODY = thisBlob3D & microglia3DBW;
        
        neuronVolume{i,obj}     = sum(double(thisBlob3D(:)));
        MGonNeuronVolume{i,obj} = sum(double(thisBlob3D_BW_MG_BODY(:)));
    end
    
    dendriteVolume(i)     = sum(double(dendrites3DBW(:)));
    MGonDendriteVolume(i) = sum(double(BW_MG_DEND(:)));
end
%% Make Plots
load('W:\Assembly\Lorenzo\Data_Microglia\L23-Coloc-Results\extractedResults.mat');

colors{1} = [1 0 0];
colors{2} = [0 1 0];
colors{3} = [0 0 1];
colors{4} = [1 0 1];
colors{5} = [1 1 0];

figure,

xvals1 = 0+randi(100,1,length(numNeurons(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(numNeurons(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(numNeurons(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(numNeurons(ageList==15)))/100;
xvals5 = 8+randi(100,1,length(numNeurons(ageList==22)))/100;

plot(xvals1, numNeurons(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, numNeurons(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, numNeurons(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, numNeurons(ageList==15),'.','Color',colors{4},'MarkerSize',30);
plot(xvals5, numNeurons(ageList==22),'.','Color',colors{5},'MarkerSize',30);

plot([0 1], ones(1,2)*mean(numNeurons(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(numNeurons(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(numNeurons(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(numNeurons(ageList==15)),'-','Color','k','LineWidth',6);
plot([8 9], ones(1,2)*mean(numNeurons(ageList==22)),'-','Color','k','LineWidth',6);


axis([-1 8 0 35]);
xticks([0.5 2.5 4.5 6.5 8.5])
xticklabels({'P4-5','P8','P11','P15','P22'});
ylabel('Number of SST Neurons in L2/3');
title('L2/3');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%%
figure,

xvals1 = 0+randi(100,1,length(MGonDendriteVolume(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(MGonDendriteVolume(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(MGonDendriteVolume(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(MGonDendriteVolume(ageList==15)))/100;
xvals5 = 8+randi(100,1,length(MGonDendriteVolume(ageList==22)))/100;

plot(xvals1, MGonDendriteVolume(ageList<=5)./dendriteVolume(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, MGonDendriteVolume(ageList==8)./dendriteVolume(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, MGonDendriteVolume(ageList==15)./dendriteVolume(ageList==15),'.','Color',colors{4},'MarkerSize',30);
plot(xvals5, MGonDendriteVolume(ageList==22)./dendriteVolume(ageList==22),'.','Color',colors{5},'MarkerSize',30);

plot([0 1], ones(1,2)*mean(MGonDendriteVolume(ageList<=5)./dendriteVolume(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(MGonDendriteVolume(ageList==8)./dendriteVolume(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(MGonDendriteVolume(ageList==15)./dendriteVolume(ageList==15)),'-','Color','k','LineWidth',6);
plot([8 9], ones(1,2)*mean(MGonDendriteVolume(ageList==22)./dendriteVolume(ageList==22)),'-','Color','k','LineWidth',6);

ranksum(MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11),...
    MGonDendriteVolume(ageList==15)./dendriteVolume(ageList==15))

axis([-1 8 0 0.06]);
xticks([0.5 2.5 4.5 6.5 8.5])
xticklabels({'P4-5','P8','P11','P15','P22'});
ylabel('MG on Dendrites / Dendrite Volume');
title('L2/3');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%%
neuronVolumeAll     = [];
MGonNeuronVolumeAll = [];
ageListAll          = [];
numNeuronsAll       = [];
ratioOfMGNeurons    = [];

for d = 1:99
        lenObj         = size(cell2mat(neuronVolume(d,:)),2);
        ageListAll     = [ageListAll ones(1,lenObj)*ageList(d)];
        numNeuronsAll  = [numNeuronsAll ones(1,lenObj)*numNeurons(d)];
        temp           = 0;
        for obj = 1:lenObj
            neuronVolumeAll     = [neuronVolumeAll neuronVolume{d,obj}];
            MGonNeuronVolumeAll = [MGonNeuronVolumeAll MGonNeuronVolume{d,obj}];
            if MGonNeuronVolume{d,obj}>0
                temp = temp + 1;
            end
        end
        ratioOfMGNeurons(d) = temp/lenObj;
end

figure,

xvals1 = 0+randi(100,1,length(ratioOfMGNeurons(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(ratioOfMGNeurons(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(ratioOfMGNeurons(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(ratioOfMGNeurons(ageList==15)))/100;
xvals5 = 8+randi(100,1,length(ratioOfMGNeurons(ageList==22)))/100;

plot(xvals1, ratioOfMGNeurons(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, ratioOfMGNeurons(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, ratioOfMGNeurons(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, ratioOfMGNeurons(ageList==15),'.','Color',colors{4},'MarkerSize',30);
plot(xvals5, ratioOfMGNeurons(ageList==22),'.','Color',colors{5},'MarkerSize',30);

ranksum(ratioOfMGNeurons(ageList==11),ratioOfMGNeurons(ageList==15))

plot([0 1], ones(1,2)*mean(ratioOfMGNeurons(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(ratioOfMGNeurons(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(ratioOfMGNeurons(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(ratioOfMGNeurons(ageList==15)),'-','Color','k','LineWidth',6);
plot([8 9], ones(1,2)*mean(ratioOfMGNeurons(ageList==22)),'-','Color','k','LineWidth',6);

axis([-1 8 0 1]);
xticks([0.5 2.5 4.5 6.5 8.5])
xticklabels({'P4-5','P8','P11','P15','P22'});
ylabel('Ratio of Neurons with MG');
title('L2/3');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%%
neuronVolumeAll     = [];
MGonNeuronVolumeAll = [];
ageListAll          = [];
numNeuronsAll       = [];
ratioOfMGNeurons    = [];

for d = 1:99
        lenObj         = size(cell2mat(dendriteVolume(d,:)),2);
        ageListAll     = [ageListAll ones(1,lenObj)*ageList(d)];
        numNeuronsAll  = [numNeuronsAll ones(1,lenObj)*numNeurons(d)];
        temp           = 0;
        for obj = 1:lenObj
            neuronVolumeAll     = [neuronVolumeAll neuronVolume{d,obj}];
            MGonNeuronVolumeAll = [MGonNeuronVolumeAll MGonNeuronVolume{d,obj}];
            if MGonNeuronVolume{d,obj}>0
                temp = temp + 1;
            end
        end
        ratioOfMGNeurons(d) = temp/lenObj;
end

figure,

xvals1 = 0+randi(100,1,length(ratioOfMGNeurons(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(ratioOfMGNeurons(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(ratioOfMGNeurons(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(ratioOfMGNeurons(ageList==15)))/100;
xvals5 = 8+randi(100,1,length(ratioOfMGNeurons(ageList==22)))/100;

plot(xvals1, ratioOfMGNeurons(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, ratioOfMGNeurons(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, ratioOfMGNeurons(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, ratioOfMGNeurons(ageList==15),'.','Color',colors{4},'MarkerSize',30);
plot(xvals5, ratioOfMGNeurons(ageList==22),'.','Color',colors{5},'MarkerSize',30);

ranksum(ratioOfMGNeurons(ageList==11),ratioOfMGNeurons(ageList==15))

plot([0 1], ones(1,2)*mean(ratioOfMGNeurons(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(ratioOfMGNeurons(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(ratioOfMGNeurons(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(ratioOfMGNeurons(ageList==15)),'-','Color','k','LineWidth',6);
plot([8 9], ones(1,2)*mean(ratioOfMGNeurons(ageList==22)),'-','Color','k','LineWidth',6);

axis([-1 8 0 1]);
xticks([0.5 2.5 4.5 6.5 8.5])
xticklabels({'P4-5','P8','P11','P15','P22'});
ylabel('Ratio of Neurons with MG');
title('L2/3');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%% Make Plots Layer 5
load('E:\LorenzoMicrogliaResults\L5\extractedResultsL5.mat');

colors{1} = [1 0 0];
colors{2} = [0 1 0];
colors{3} = [0 0 1];
colors{4} = [1 0 1];
colors{5} = [1 1 0];

figure,

xvals1 = 0+randi(100,1,length(numNeurons(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(numNeurons(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(numNeurons(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(numNeurons(ageList==21)))/100;

plot(xvals1, numNeurons(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, numNeurons(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, numNeurons(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, numNeurons(ageList==21),'.','Color',colors{4},'MarkerSize',30);

plot([0 1], ones(1,2)*mean(numNeurons(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(numNeurons(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(numNeurons(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(numNeurons(ageList==21)),'-','Color','k','LineWidth',6);

axis([-1 6 0 70]);
xticks([0.5 2.5 4.5 6.5])
xticklabels({'P4-5','P8','P11','P21'});
ylabel('Number of SST Neurons in L5');
title('L5');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%%
figure,

xvals1 = 0+randi(100,1,length(MGonDendriteVolume(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(MGonDendriteVolume(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(MGonDendriteVolume(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(MGonDendriteVolume(ageList==21)))/100;

plot(xvals1, MGonDendriteVolume(ageList<=5)./dendriteVolume(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, MGonDendriteVolume(ageList==8)./dendriteVolume(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, MGonDendriteVolume(ageList==21)./dendriteVolume(ageList==21),'.','Color',colors{4},'MarkerSize',30);

plot([0 1], ones(1,2)*mean(MGonDendriteVolume(ageList<=5)./dendriteVolume(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(MGonDendriteVolume(ageList==8)./dendriteVolume(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(MGonDendriteVolume(ageList==21)./dendriteVolume(ageList==21)),'-','Color','k','LineWidth',6);


ranksum(MGonDendriteVolume(ageList==8)./dendriteVolume(ageList==8),...
    MGonDendriteVolume(ageList==11)./dendriteVolume(ageList==11))



axis([-1 6 0 0.06]);
xticks([0.5 2.5 4.5 6.5])
xticklabels({'P4-5','P8','P11','P21'});
ylabel('MG on Dendrites / Dendrite Volume');
title('L5');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
%%
neuronVolumeAll     = [];
MGonNeuronVolumeAll = [];
ageListAll          = [];
numNeuronsAll       = [];
ratioOfMGNeurons    = [];

for d = 1:159
        lenObj         = size(cell2mat(neuronVolume(d,:)),2);
        ageListAll     = [ageListAll ones(1,lenObj)*ageList(d)];
        numNeuronsAll  = [numNeuronsAll ones(1,lenObj)*numNeurons(d)];
        temp           = 0;
        for obj = 1:lenObj
            neuronVolumeAll     = [neuronVolumeAll neuronVolume{d,obj}];
            MGonNeuronVolumeAll = [MGonNeuronVolumeAll MGonNeuronVolume{d,obj}];
            if MGonNeuronVolume{d,obj}>0
                temp = temp + 1;
            end
        end
        ratioOfMGNeurons(d) = temp/lenObj;
end

figure,

xvals1 = 0+randi(100,1,length(ratioOfMGNeurons(ageList<=5)))/100;
xvals2 = 2+randi(100,1,length(ratioOfMGNeurons(ageList==8)))/100;
xvals3 = 4+randi(100,1,length(ratioOfMGNeurons(ageList==11)))/100;
xvals4 = 6+randi(100,1,length(ratioOfMGNeurons(ageList==21)))/100;

plot(xvals1, ratioOfMGNeurons(ageList<=5),'.','Color',colors{1},'MarkerSize',30);
hold on;
plot(xvals2, ratioOfMGNeurons(ageList==8),'.','Color',colors{2},'MarkerSize',30);
plot(xvals3, ratioOfMGNeurons(ageList==11),'.','Color',colors{3},'MarkerSize',30);
plot(xvals4, ratioOfMGNeurons(ageList==21),'.','Color',colors{4},'MarkerSize',30);





plot([0 1], ones(1,2)*mean(ratioOfMGNeurons(ageList<=5)),'-','Color','k','LineWidth',6);
plot([2 3], ones(1,2)*mean(ratioOfMGNeurons(ageList==8)),'-','Color','k','LineWidth',6);
plot([4 5], ones(1,2)*mean(ratioOfMGNeurons(ageList==11)),'-','Color','k','LineWidth',6);
plot([6 7], ones(1,2)*mean(ratioOfMGNeurons(ageList==21)),'-','Color','k','LineWidth',6);

axis([-1 6 0 1]);
xticks([0.5 2.5 4.5 6.5])
xticklabels({'P4-5','P8','P11','P21'});
ylabel('Ratio of Neurons with MG');
title('L5');
set(gca,'TickDir','out');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);

%%
mipTemp    = max(info.cube2,[],3);

figure, imshow(bwCandTemp);
figure, imshow(closeBW);

[bw1,bw2,bw3]        = filterWrongNeurons(closeBW,mipTemp,Mdl);
figure, imshow(bw1);
figure, imshow(bw2);
figure, imshow(bw3);

%%
% Label the blobs.


figure, imshow(binaryImage)
%%










microglia_mip = max(cube1{i},[],3);
thr = multithresh(microglia_mip,15);
temp = zeros(size(microglia_mip));
for i = 1:15
    temp = temp + bwmorph(double(microglia_mip>thr(i)),'skel','Inf');
end
figure, imshow(V1enhanced,[]);
figure, imshow(bilateralFilter(imadjust(mip)),[]);
figure, imshow(bwCandidates1,[]);




%%

bw1_3d = repmat(double(BW1),1,1,size(I1norm,3));
bw2_3d = repmat(double(BW2),1,1,size(I1norm,3));

I1_bw_3d = bw1_3d.*I1norm;
I2_bw_3d = bw2_3d.*I2norm;

IlineBW      = I1_bw_3d(:);
ClineBW      = reshape(Iline ,[size(I1norm,1)*size(I1norm,2) size(I1norm,3)]);

IlineBW2     = I2_bw_3d(:);
ClineBW2     = reshape(Iline2,[size(I2norm,1)*size(I2norm,2) size(I2norm,3)]);

k = 5;

stream  = RandStream('mlfg6331_64');
options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
[idx,c] = kmeans(sort(ClineBW,2),k,'Options',options,'MaxIter',1000,'Distance','cityblock','Display','final','Replicates',5);

allc    = trapz(c');
[s,sx]  = sort(allc);

idxnew = idx;

for i = 1:k
    idxnew(idx==sx(i)) = i;
end
idxnew = (idxnew - 1);%/(k-1);

figure, plot(c');

Cline_bw = reshape(idxnew,[size(I1norm,1) size(I1norm,2) 1]);
figure   , imagesc(Cline_bw);









rgb12(:,:,1) = adapthisteq(max(I2norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');
rgb12(:,:,2) = adapthisteq(max(I1norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');
rgb12(:,:,3) = double(BW1);
figure, imshow(rgb12,[]);

rgb12(:,:,1) = adapthisteq(max(I2norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');
rgb12(:,:,2) = adapthisteq(max(I1norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');
rgb12(:,:,3) = double(BW2);
figure, imshow(rgb12,[]);

bw12_3d = repmat(double(BW1&BW2),1,1,size(I1norm,3));

Iline_bw12_3d       = bw12_3d(:);
Cline_bw12_3d       = reshape(Iline_bw12_3d,[size(I1norm,1)*size(I2norm,2) size(I2norm,3)]);

All_lines_3d_1_11   = Cline_bw12_3d.*Cline;
All_lines_3d_2_12   = Cline_bw12_3d.*Cline2;

[B_1_11,I_1_11]     = sortrows(All_lines_3d_1_11);
[B_2_12,I_2_12]     = sortrows(All_lines_3d_2_12);

startz_1_11         = find(sum(B_1_11,2)>0);
startz_2_12         = find(sum(B_2_12,2)>0);

vals_1_11           = B_1_11(startz_1_11(1):end,:);
vals_2_12           = B_2_12(startz_2_12(1):end,:);

locs_1_11           = I_1_11(startz_1_11(1):end,:);
locs_2_12           = I_2_12(startz_2_12(1):end,:);

% [IDX_1_11 , C_1_11 ] = kmeans(vals_1_11, 50, 'Distance', 'sqeuclidean');
% [IDX_2_12 , C_2_12 ] = kmeans(vals_2_12, 50, 'Distance', 'sqeuclidean');
% [IDX_1_13 , C_1_13 ] = kmeans(vals_1_13, 50, 'Distance', 'sqeuclidean');
% [IDX_3_13 , C_3_13 ] = kmeans(vals_3_13, 50, 'Distance', 'sqeuclidean');

[m_1_11,l_1_11] = max(vals_1_11,[],2);
[m_2_12,l_2_12] = max(vals_2_12,[],2);

d12 = abs(l_1_11 - l_2_12);

[d12_val,d12_locs] = find(d12<=z_param);

candLocs_1_11 = locs_1_11(d12_val);
candLocs_2_12 = locs_2_12(d12_val);

[x_1_11,y_1_11] = convert_1D_to_2D_locs(candLocs_1_11,size(I1norm,1));
BWnew_1_11 = zeros(size(I1norm,1),size(I1norm,1));
for i = 1:length(x_1_11)
    BWnew_1_11(y_1_11(i),x_1_11(i)) = 1;
end

[x_2_12,y_2_12] = convert_1D_to_2D_locs(candLocs_2_12,size(I1norm,1));
BWnew_2_12 = zeros(size(I1norm,1),size(I1norm,1));
for i = 1:length(x_2_12)
    BWnew_2_12(y_2_12(i),x_2_12(i)) = 1;
end


I1_max = adapthisteq(max(I1norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');
I2_max = adapthisteq(max(I2norm,[],3),'NumTiles',[8 8],'ClipLimit',0.005, 'Distribution', 'uniform');

I = imfuse(I1_max,I2_max,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);

MIP_candidates_12 = BWnew_1_11 | BWnew_2_12;
[xloc,yloc]       = find(MIP_candidates_12);

all_x = xloc;
all_y = yloc;

z_traces_red    = zeros(length(xloc),size(I2norm,3));
z_traces_green  = zeros(length(xloc),size(I1norm,3));

for i = 1:length(xloc)
    z_traces_red(i,:)   = I2norm(yloc(i),xloc(i),:);
    z_traces_green(i,:) = I1norm(yloc(i),xloc(i),:);
    %hold on; plot(z_traces_red(i,:),'r'); 
    %hold on; plot(z_traces_green(i,:),'g');
end


% candidates = uint8(255*double(candidates));
%%
I_fused_all_rgb(:,:,2) = I1_max;
I_fused_all_rgb(:,:,1) = I2_max;
I_fused_all_rgb(:,:,3) = MIP_candidates_12;


im_obj = squeeze(I_fused_all_rgb(:,:,3));
cc     = bwconncomp(im_obj);
cc_all = cc;

% numPixels = cellfun(@numel,cc.PixelIdxList);
% [biggest,idx] = max(numPixels);
% BW(CC.PixelIdxList{idx}) = 0;

% figure,
% imagesc(I_fused_all);    

data.info           = rTemp{2};
data.numObjects     = cc.NumObjects;
data.pixelIdxList   = cc.PixelIdxList;
data.version        = 8;


prompt1         = {'Enter the matrix size:'};
name1           = 'Input for Number of Division for the Image Magnification';
numlines1       = 1;
defaultanswer1  = {'4'};
   
answer      = inputdlg(prompt1,name1,numlines1,defaultanswer1);
ind_div     = str2double(answer{1});

block_size  = floor(size(im_obj,1)/ind_div);

for ind = 1:ind_div
    indxy{ind} = (1+(ind-1)*block_size:block_size*ind);
end


status = 'NotDone';

data.INB_MGP    = 0;
data.INB_MGB    = 0;
data.INP_MGP    = 0;
data.INbP_MGP   = 0;

coloc_counter   = 0;

x_pts           = [];
y_pts           = [];
prev_x          = [];
prev_y          = [];
position_labels = [];

colorVectors = {'[1 0 0]','[0 1 0]','[0 0 1]','[1 0 1]'};

col = {'ro','bs','g^','ch'};

while  strcmp(status,'NotDone')
    
        str     = {'INBody-MGProcess','INBody-MGBody','INProcess-MGProcess','INbasalProcess-MGProcess','Done'};
        [s,v]   = listdlg('PromptString','Select:','SelectionMode','single','ListString',str);

        if strcmp(str{s},'Done')
            status = [];
            status = 'Done';
            break;
        end

        counter_int = 0;
        
        for i = 1:ind_div
            for j = 1:ind_div

                newSpineSelectFig = figure;
                imagesc([I(indxy{i},indxy{j},:) , uint8(255*(double(I_fused_all_rgb(indxy{i},indxy{j},:))))]);
                
                if coloc_counter > 0
                    
                    hold on;
                    col = {'ro','bs','g^','ch'};
                    
                    for k = 1:4
                        prev_x = x_pts - (j-1)*block_size + block_size;
                        prev_y = y_pts - (i-1)*block_size;
                        scatter(prev_x(position_labels == k),prev_y(position_labels == k),40,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',colorVectors{k},'LineWidth',2);
                    end 
                end
                
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                title(sprintf('%s %s %s','Choose',str{s},'Interactions'));
                label = str{s};

                [x_pts_temp, y_pts_temp] = getpts_SpineS(newSpineSelectFig);
                
                x_pts_temp = x_pts_temp + (j-1)*block_size - block_size;
                y_pts_temp = y_pts_temp + (i-1)*block_size;
                
                x_pts = [x_pts; x_pts_temp];
                y_pts = [y_pts; y_pts_temp];
                
                close(newSpineSelectFig);

                counter_int = counter_int + length(x_pts_temp);
                cla reset  
            end
        end
        
        position_labels = [position_labels; s*ones(1,counter_int)'];

        switch s
            case 1
                data.INB_MGP  = counter_int;
            case 2
                data.INB_MGB  = counter_int;
            case 3
                data.INP_MGP  = counter_int;
            case 4
                data.INbP_MGP = counter_int;
        end
        
        coloc_counter = coloc_counter + 1;
end

data.labels             = position_labels;
data.positions.x_pts    = all_x;
data.positions.y_pts    = all_y;
data.connComp           = cc_all;

if ismac
        st = regexp(D1{1}, '/', 'split');
        folderName = fullfile('/',st{1:end-1},sprintf('%s%s','results',datestr(now,30))); % Date Format (ISO 8601)  'yyyymmddTHHMMSS'
    else
        st = regexp(D1{1}, '\', 'split');
        folderName = fullfile(st{1:end-1},sprintf('%s%s','results',datestr(now,30))); % Date Format (ISO 8601)  'yyyymmddTHHMMSS'
end

figure,
numberOfXTicks = 4;
h       = stem([data.INB_MGP,data.INB_MGB,data.INP_MGP,data.INbP_MGP],'.');
xData   = get(h,'XData');
set(gca,'Xtick',linspace(xData(1),xData(end),numberOfXTicks))
set(gca, 'XTickLabel',{'INB-MGP','INB-MGB','INP-MGP','INbP-MGP'});
rotateXLabels( gca(), 45 );
axis([0 5 0 max([data.INB_MGP,data.INB_MGB,data.INP_MGP,data.INbP_MGP])+2])
title(st{7});

saveas(h, fullfile('E:\DataTemp\Lorenzo',sprintf('%s%s',st{end},'_numbers.png')), 'png');

%%
figure,
h2 = imshow(uint8(255*double(I_fused_all_rgb)),[]);
hold on

col = {'ro','bs','g^','ch'};

%%
for i = 1:4
    figure;
    labels_numb = data.labels==i;
    labels_numb = sort(labels_numb,'descend');
    scatter(data.positions.x_pts(labels_numb),data.positions.y_pts(labels_numb),col{i});
end
%%
% %%
% for i = 1:3
%     figure;
%     scatter(data.positions.x_pts(data.labels==i),data.positions.y_pts(data.labels==i),col{i});
% end
%%
legend('INB-MGP','INB-MGB','INP-MGP','INbP-MGP','NorthEastOutside');
title(st{7});

saveas(h2, fullfile('E:\DataTemp\Lorenzo',sprintf('%s%s',st{end},'_image.png')), 'png');

%%
% mkdir(folderName);
% af = '/Volumes/workspace/Assembly/Ali/Confocal';
% [status,msg,msgID] = mkdir(af)  

save(fullfile('E:\DataTemp\Lorenzo',sprintf('%s_%s%s',st{end},datestr(now,30),'_data.mat')),'-struct', 'data');

% 
% 
% %%
% % Create a figure and axes
% f = figure('Visible','off');
% ax = axes('Units','pixels');
% imagesc(I_fused_all);    
% % Create pop-up menu 
% popup = uicontrol('Style', 'popup',...
%        'String', {'INBody-MGProcess','INBody-MGBody','INProcess-MGProcess','INProcess-MGBody'},...
%        'Position', [20 340 100 50]);    
% 
% % Create push button
% btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
%     'Position', [20 20 50 20],...
%     'Callback', 'cla');       
% 
% f.Visible = 'on';
% 
% % 
% %    % Create slider
% %     sld = uicontrol('Style', 'slider',...
% %         'Min',1,'Max',50,'Value',41,...
% %         'Position', [400 20 120 20],...
% %         'Callback', @surfzlim); 
% 
% %     % Add a text uicontrol to label the slider.
% %     txt = uicontrol('Style','text',...
% %         'Position',[400 45 120 20],...
% %         'String','Vertical Exaggeration');
% 
% % Make figure visble after adding all components
% % This code uses dot notation to set properties. 
% % Dot notation runs in R2014b and later.
% % For R2014a and earlier: set(f,'Visible','on');
% % 
% % function setmap(source,event)
% %     val = source.Value;
% %     maps = source.String;
% %     % For R2014a and earlier: 
% %     % val = get(source,'Value');
% %     % maps = get(source,'String'); 
% % 
% %     newmap = maps{val};
% %     colormap(newmap);
% % end
% 
% 
% 
% % 
% % imshowpair(I,I_fused_all);
% % % 
% % % new_img_size = 256;
% % % I1_max_small  = imresize(I1_max ,[new_img_size new_img_size]);
% % % I2_max_small  = imresize(I2_max ,[new_img_size new_img_size]);
% % 
% % figure, imshow(I1_max,[]);
% % figure, imshow(I2_max,[]);
% % 
% % 
% % zoom_ind_x{1} = 1:512;
% % 
% % 
% % figure,
% % imshow(I,[]);
% % 
% % newSelectFig = figure;
% % H = uicontrol('Style', 'PushButton', ...
% %                      'String', 'Break', ...
% %                      'Callback', 'delete(gcbf)');
% % 
% % x_pts = [];
% % y_pts = [];
% % 
% % while (ishandle(H)) 
% % 
% %     [x_pts_temp, y_pts_temp] = getpts(newSelectFig);
% %     
% %     x_pts = [x_pts; x_pts_temp];
% %     y_pts = [y_pts; y_pts_temp];
% %     
% % end
% % 
% % close(newSelectFig);
% % cla reset
% % 
% % 
% % % 
% % %     z_profile_I1 = double(squeeze(I1(floor(y_pts_temp),floor(x_pts_temp),:))); 
% % %     z_profile_I2 = double(squeeze(I2(floor(y_pts_temp),floor(x_pts_temp),:)));    
% % % 
% % %     figure, plot(z_profile_I1)
% % %     hold on,
% % %     plot(z_profile_I2,'r');
% % % 
% % %     [val_1,CZ_1] = max(smooth(z_profile_I1));
% % %     [val_2,CZ_2] = max(smooth(z_profile_I2));
% % %     
% %     hold on;
% %     plot(x_pts,y_pts,'r.');
% %  
% %     x_pts = [x_pts; x_pts_temp];
% %     y_pts = [y_pts; y_pts_temp];
% % 
% % end
% % 
% % 
% % handles.newSpineSelectFig = figure;
% % imshow(handles.firstMIP,[]);
% %  
% % if handles.Flags.selectROI == 1
% %     hold on;
% %     plot(handles.x_pts,handles.y_pts,'r.');
% % end
% %  
% % [x_pts, y_pts] = getpts_SpineS(handles.newSpineSelectFig);
% % handles.x_pts = [handles.x_pts; x_pts];
% % handles.y_pts = [handles.y_pts; y_pts];
% %  
% % close(handles.newSpineSelectFig);
% %  
% % cla reset
% % 
% % bref_change = I;
% % [a,b,c] = size(I);
% % r1A    =    [254 255 0];
% % right_wing = rgb2gray(I);
% % right_wing  = imgaussfilt(right_wing, 4);
% % figure; imshow(right_wing);
% % 
% % for x=1:a    
% %     for y=1:b
% %         if     (bref_change(x,y,1) == r1A(1)) && (bref_change(x,y,2) == r1A(2)) && (bref_change(x,y,3) == r1A(3))
% %                 right_wing(x,y,:) = 255;
% %        % else
% %                 %right_wing(x,y,:) = 0;
% %         end
% %     end
% %     %
% % end
% % se = offsetstrel('ball',3,3);
% % dilatedI = imdilate( right_wing ,se);
% % figure; imshow(dilatedI,[])









