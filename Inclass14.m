%Inclass 14
clear all
x = 0;
%Work with the image stemcells_dapi.tif in this folder

% (1) Make a binary mask by thresholding as best you can
reader = bfGetReader('stemcells_dapi.tif');
sizeMat = [reader.getSizeX, reader.getSizeY];
times = reader.getSizeT;
chan = reader.getSizeC;
zstacks = reader.getSizeZ;
iplane = reader.getIndex(zstacks-1,chan-1,times-1)+1;
img = bfGetPlane(reader,iplane);
%imshow(img, []); %quick check on the image.

%background subtract
img_sm = imfilter(img,fspecial('gaussian',4,2));
img_bg = imopen(img_sm,strel('disk',100));
img_sm_bgsub = imsubtract(img_sm,img_bg);
%figure(x);imshow(img_sm_bgsub,[]);

mask1 = img_sm_bgsub > 85;
mask1 = imerode(mask1,strel('disk',3));
x = x+1;figure(x);imshow(mask1,[]);



% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask
cc = bwconncomp(mask1);
stats = regionprops(cc, 'area');
area = [stats.Area];
fusedCells = area > mean(area) + std(area);
sublist = cc.PixelIdxList(fusedCells);
sublist = cat(1,sublist{:});
fusedMask = false(size(mask1));
fusedMask(sublist) = 1;
%x = x+1;figure(x);imshow(fusedMask,[]);

%erosion
nucmin = imerode(fusedMask, strel('disk',7));
%x = x+1;figure(x);imshow(nucmin,[]);
%define the boundries 
outside = ~imdilate(fusedMask,strel('disk',1));
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin,nucmin | outside);
%pcolor(basin);shading flat;

L = watershed(basin);
%x = x+1;figure(x);imshow(L);colormap('jet');caxis([0 20]);

mask2 = L > 1 | (mask1 - fusedMask);


%repeat to pick up the missing cell clusters from before. 
cc = bwconncomp(mask2);
stats = regionprops(cc, 'area');
area = [stats.Area];
fusedCells = area > mean(area) + std(area);
sublist = cc.PixelIdxList(fusedCells);
sublist = cat(1,sublist{:});
fusedMask = false(size(mask2));
fusedMask(sublist) = 1;


%erosion. 
nucmin = imerode(fusedMask, strel('disk',7));
%define the boundries 
outside = ~imdilate(fusedMask,strel('disk',1));
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin,nucmin | outside);

L = watershed(basin);

mask2 = L > 1 | (mask2 - fusedMask);
x = x+1;figure(x);imshow(mask2);

%(B) with a  distance transform. Which works better in this case?

for ii = 1:3
cc = bwconncomp(mask1);
stats = regionprops(cc, 'area');
area = [stats.Area];
fusedCells = area > mean(area) + std(area);
sublist = cc.PixelIdxList(fusedCells);
sublist = cat(1,sublist{:});
fusedMask = false(size(mask1));
fusedMask(sublist) = 1;

D = bwdist(~fusedMask);
D = -D;
D(~fusedMask) = -inf;
[]);
%define the boundries 


L = watershed(basin);
x = x+1;figure(x);imshow(L,[]);

mask1 = L > 1 | (mask1 - fusedMask);
end
x = x+1;figure(x);imshow(mask1,[]);

%Adam Howard: In this particular instance, it looks like the erosion method
%did a better job of separating the joint cells. 
