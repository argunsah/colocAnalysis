function result = thresholdWithZprofile(I1norm,BW)

Iline    = I1norm(:);
Cline    = reshape(Iline ,[size(I1norm,1)*size(I1norm,2) size(I1norm,3)]);

BWline   = double(BW(:));
CBW      = reshape(BWline ,[size(I1norm,1)*size(I1norm,2) size(I1norm,3)]);

Cline    = Cline.*CBW;
sumCline = sum(Cline,2);
[ind,val]= find(sumCline>0);

thr      = ones(size(sumCline));
for i = ind
    thr(i)  = graythresh(Cline(i,:));    
end

thrMatrix = squeeze(repmat(thr,1,1,size(Cline,2)));

bwCline = Cline>thrMatrix;

result = reshape(bwCline,[size(I1norm,1) size(I1norm,2) size(I1norm,3)]);