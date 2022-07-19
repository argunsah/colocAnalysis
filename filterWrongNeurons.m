function [bwCandidatesNew1,bwCandidatesNew2,bwCandidatesNew3] = filterWrongNeurons(bwCandidates,mip,Mdl)

rpC             = regionprops(bwCandidates,mip,'All');

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
CC                = bwconncomp(bwCandidates);
bwCandidatesNew1  = zeros(size(bwCandidates));
bwCandidatesNew2  = zeros(size(bwCandidates));
bwCandidatesNew3  = zeros(size(bwCandidates));

numObjects = size(Areas2D,1);

for obj = 1:numObjects
    if label(obj) == 1
        bwCandidatesNew1(CC.PixelIdxList{obj}) = 1;
    end
    if label(obj) == 2
        bwCandidatesNew2(CC.PixelIdxList{obj}) = 1;
    end
    if label(obj) == 3
        bwCandidatesNew3(CC.PixelIdxList{obj}) = 1;
    end
end