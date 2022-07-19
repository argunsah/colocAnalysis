function [thr1, thr2, locthr1, locthr2, BWnew1,BWnew2] = zProfilesUsingColocCandids(I1norm,I2norm,temp1,temp2)

% I1norm = cube{1};
% I2norm = cube{2};
% 
% BW1 = BW{1};
% BW2 = BW{2};

Iline       = I1norm(:);
Cline       = reshape(Iline ,[size(I1norm,1)*size(I1norm,2) size(I1norm,3)]);
Iline2      = I2norm(:);
Cline2      = reshape(Iline2,[size(I2norm,1)*size(I2norm,2) size(I2norm,3)]);

BW1 = repmat(temp1,[1,1,size(I1norm,3)]);
BW2 = repmat(temp2,[1,1,size(I1norm,3)]);

bw12_3d           = BW1&BW2;
max_bw12_3d       = max(bw12_3d,[],3);

Iline_bw12_3d     = bw12_3d(:);
Cline_bw12_3d     = reshape(Iline_bw12_3d,[size(I1norm,1)*size(I2norm,2) size(I2norm,3)]);

All_lines_3d_1_11 = double(Cline_bw12_3d).*Cline;
All_lines_3d_2_12 = double(Cline_bw12_3d).*Cline2;

[B_1_11,I_1_11]   = sortrows(All_lines_3d_1_11);
[B_2_12,I_2_12]   = sortrows(All_lines_3d_2_12);

startz_1_11       = find(sum(B_1_11,2)>0);
startz_2_12       = find(sum(B_2_12,2)>0);

vals_1_11         = B_1_11(startz_1_11(1):end,:);
vals_2_12         = B_2_12(startz_2_12(1):end,:);
locs_1_11         = I_1_11(startz_1_11(1):end,:);
locs_2_12         = I_2_12(startz_2_12(1):end,:);

[m_1_11,l_1_11] = max(vals_1_11,[],2);
[m_2_12,l_2_12] = max(vals_2_12,[],2);

zs_vals_1_11 = zscore(vals_1_11);
zs_vals_2_12 = zscore(vals_2_12);
% zs_vals_1_11_movmean = movmean(zs_vals_1_11,3,2);
% zs_vals_2_12_movmean = movmean(zs_vals_2_12,3,2);
zs_vals_1_11_movmean = zs_vals_1_11;
zs_vals_2_12_movmean = zs_vals_2_12;

norm_zs_vals_1_11_movmean = zs_vals_1_11_movmean - min(zs_vals_1_11_movmean(:));
norm_zs_vals_2_12_movmean = zs_vals_2_12_movmean - min(zs_vals_2_12_movmean(:));

norm_zs_vals_1_11_movmean = norm_zs_vals_1_11_movmean/max(norm_zs_vals_1_11_movmean(:));
norm_zs_vals_2_12_movmean = norm_zs_vals_2_12_movmean/max(norm_zs_vals_2_12_movmean(:));

locthr1 = zeros(size(norm_zs_vals_1_11_movmean));
locthr2 = zeros(size(norm_zs_vals_2_12_movmean));

for cc = 1:size(norm_zs_vals_1_11_movmean,1)
    t1(cc) = graythresh(norm_zs_vals_1_11_movmean(cc,:));
    t2(cc) = graythresh(norm_zs_vals_2_12_movmean(cc,:));
    
    locthr1(cc,:) = norm_zs_vals_1_11_movmean(cc,:) > t1(cc);
    locthr2(cc,:) = norm_zs_vals_2_12_movmean(cc,:) > t2(cc);
end

gthr1 = graythresh(norm_zs_vals_1_11_movmean(:));
gthr2 = graythresh(norm_zs_vals_2_12_movmean(:));

thr1  = norm_zs_vals_1_11_movmean > gthr1;
thr2  = norm_zs_vals_2_12_movmean > gthr2;

% %%
% cc = 22301;
% figure, plot(norm_zs_vals_2_12_movmean(cc,:))
% hold on, plot(locthr2(cc,:),'r')
% hold on, plot(thr2(cc,:),'.g')
% hold on, plot(thr2(cc,:)&locthr2(cc,:),'k')
% 
% axis([1 45 -1 1.1])
%%




thr1combined = thr1 & locthr1;
thr2combined = thr2 & locthr2;

[x_1_11,y_1_11] = convert_1D_to_2D_locs(locs_1_11,size(I1norm,1));
[x_2_12,y_2_12] = convert_1D_to_2D_locs(locs_2_12,size(I2norm,1));

BWnew1 = zeros(size(I1norm));
BWnew2 = zeros(size(I2norm));

for c = 1:length(x_1_11)
    BWnew1(y_1_11(c),x_1_11(c),:) = thr1combined(c,:);
    BWnew2(y_2_12(c),x_2_12(c),:) = thr2combined(c,:);
end

% figure, imagesc(thr1&thr2);
% 
% [idx1,c1]           = kmeans(sort(norm_zs_vals_1_11_movmean,2),2);
% [sumvals1,indVal1]  = max(sum(c1,2));
% 
% [idx2,c2]           = kmeans(sort(norm_zs_vals_2_12_movmean,2),2);
% [sumvals2,indVal2]  = max(sum(c2,2));
% 
% indx1 = sum(thr1,2);
% 
% candLocs_1_11       = locs_1_11(idx1==1);
% % candLocs_2_12       = locs_2_12(idx1==indVal2);
% 
% [x_1_11,y_1_11]     = convert_1D_to_2D_locs(candLocs_1_11,size(I1norm,1));
% BWnew_1_11          = zeros(size(I1norm));
% 
% for i = 1:length(x_1_11)
%     BWnew_1_11(y_1_11(i),x_1_11(i)) = 1;
% end
% 
% [x_2_12,y_2_12]     = convert_1D_to_2D_locs(candLocs_2_12,size(I1norm,1));
% BWnew_2_12          = zeros(size(I2norm));
% 
% for i = 1:length(x_2_12)
%     BWnew_2_12(y_2_12(i),x_2_12(i)) = 1;
% end














% 
% 
% 
% 
% d12 = abs(l_1_11-l_2_12);
% 
% [d12_val,d12_locs] = find(combined>=2);
% 
% candLocs_1_11 = locs_1_11(d12_val);
% candLocs_2_12 = locs_2_12(d12_val);
% 
% [x_1_11,y_1_11] = convert_1D_to_2D_locs(candLocs_1_11,1600);
% BWnew_1_11 = zeros(1600,1600);
% for i = 1:length(x_1_11)
%     BWnew_1_11(y_1_11(i),x_1_11(i)) = 1;
% end
% 
% [x_2_12,y_2_12] = convert_1D_to_2D_locs(candLocs_2_12,1600);
% BWnew_2_12 = zeros(1600,1600);
% for i = 1:length(x_2_12)
%     BWnew_2_12(y_2_12(i),x_2_12(i)) = 1;
% end