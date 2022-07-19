function binaryImage = filterFunkyShapes(binaryImage)

labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage,'Area','Perimeter');
% Do size filtering and roundness filtering.
% Get areas and perimeters of all the regions into single arrays.
allAreas = [measurements.Area];
allPerimeters = [measurements.Perimeter];
% Compute circularities.
circularities = allPerimeters.^2 ./ (4*pi*allAreas);
% Find objects that have "round" values of circularities.
maxAllowableArea = 30000;
keeperBlobs = circularities < 15 & allAreas < maxAllowableArea & allAreas > 50; % Whatever values you want.
% Get actual index numbers instead of a logical vector
% so we can use ismember to extract those blob numbers.
roundObjects = find(keeperBlobs);
% Compute new binary image with only the small, round objects in it.
binaryImage = ismember(labeledImage, roundObjects) > 0;