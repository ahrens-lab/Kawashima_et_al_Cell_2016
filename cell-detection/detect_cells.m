% detect_cells.m
%
% This scripts detect circular-shaped cells in large stacks based on local
% contrast information.  Detailes are described in Kawashima et al., 2016.
%
% Requires Image Processing Toolbox.  Also requires other funcitons 
% contained in the "functions" folder.
% 
% This script only runs on 64-bit Windows because it uses
% MEX files build by VS2010. For use in other systems please
% re-compile corresponding cpp files.
%
%
% Parameters for cell detection
%
% outdir         -> Output directory of cell recognition result.
% fname          -> Name of the image file to be analyzed.
% br_threshold   -> Minimum brightness threshold for the cell to be detected
% cont_threshold -> Threshold for local contrast. This sets how the cell
%                   should be brighter than surrounding darkest point.
% cell_rad       -> Radius of cells to be detected.
%
%
% Output of the script
%
% cellmask.tif   -> Visual result of cell recognition.
% cell_info.mat  -> Information of detected cells.
%                   [center]   y and x coordinate of the detected cell
%                   [area]     area of the detected cell 
%                   [slice]    Z-plane of the detected cell 
%                   [inds]     one-dimensional indices of the detectec cell
%                   [y_minmax] min and max end of y coordinates
%                   [x_minmax] min and max end of x coordinates
%
% developed by Takashi Kawashima, HHMI Janelia Research Campus
% Sep 19,2016



outdir='C:\Users\kawashimat\Documents\GitHub\Kawashima_et_al_Cell_2016\cell-detection\sample';  %% specify output directory
fname='C:\Users\kawashimat\Documents\GitHub\Kawashima_et_al_Cell_2016\cell-detection\sample\ave.tif'; %% specify image file to recognize

br_threshold=120;
cont_threshold=7;
cell_rad=4;

%%  set filters and matrices for detecting cells;

info = imfinfo(fname);
num_images = numel(info);
dim=[info(1).Height,info(2).Width, numel(info)];

ave_rad=round(cell_rad/2)+1;
[avedisk,  ave_se, r1, c1, maskinds]=make_recog_disk(ave_rad,dim);
[maxdisk,  max_se, r2, c2, maxinds]=make_recog_disk(round(cell_rad*1.5),dim);
onedisk=makeDisk(ave_rad,ave_rad*2+1);
one_se=strel(onedisk);

r=cell_rad*2;
dimp=[dim(1)+r*2 dim(2)+r*2];
oop=zeros(dimp);
oop(r+1:end-r,r+1:end-r)=1;
one_inds=find(oop);

[mdisk,  ~, ~, ~, rankinds]=make_recog_disk(r,dimp);
rank_ones=double(maskones2D_mex(int32([dim(1) dim(2)]),int32(mdisk),int32(size(mdisk))))';

%%  detect cells


mkdir(outdir);
cell_info=struct;
cell_color=zeros([dim(1) dim(2)*2+1 dim(3) 3],'uint8');
cellnum=0;
allmask=zeros(dim(1),dim(2));
imlen=dim(1)*dim(2);

for z=1:dim(3)

    im=double(imread(fname,z));  
    allmask(:)=0;
    contimage = local_contrast_mex(single(im),int32(32),single(cont_threshold));
    contimage = imdilate(imerode(contimage,one_se),one_se);
    contimage = contimage.*uint8((im>br_threshold));

    candidates=find(contimage);

    if ~isempty(candidates)
        %% recognizing cells in the first round   

        imrank = calc_local_rank(im,rank_ones,cell_rad*2,oop,one_inds,rankinds,candidates);
        aveimg = double(local_average_mex(single(imrank),int32(c1),int32(r1),int32(candidates)));
        maximg = double(local_max_mex(single(aveimg),int32(c2),int32(r2),int32(candidates)));

        inds=find(maximg(candidates)>0 & aveimg(candidates) >0.4);
        mask2=zeros(dim(1),dim(2));

        for i=1:length(inds)
            cinds=candidates(inds(i))+maskinds;
            cinds(cinds > imlen | cinds<1)=[];
            mask2(cinds)=1;
        end

        allmask=mask2;

        %% recognizing cells in the second round            

        mask3=ones(size(im),'uint8')-imdilate(uint8(allmask),max_se);
        mask3 = imdilate(imerode(mask3,one_se),one_se);
        candidates2=candidates(mask3(candidates)>0);

        imrank2 = calc_local_rank(im,rank_ones,cell_rad*2,oop,one_inds,rankinds,candidates2);
        aveimg2 = double(local_average_mex(single(imrank2),int32(c1),int32(r1),int32(candidates2)));        
        maximg2 = double(local_max_mex(single(aveimg2),int32(c1),int32(r1),int32(candidates2)));

        inds=find(maximg2(candidates2)>0 & aveimg2(candidates2) >0.4);
        mask2=zeros(dim(1),dim(2));

        for i=1:length(inds)
            cinds=candidates2(inds(i))+maskinds;
            cinds(cinds > imlen | cinds<1)=[];
            mask2(cinds)=1;
        end

        allmask=allmask+mask2;

        %% create each cell ROIs

        [celllabel, totcell]=bwlabel(allmask,8);
        if totcell>0
                cell_info=create_cell_info(cell_info,celllabel, totcell,z);
        end
    else

        totcell=0;
        celllabel=zeros(size(im));
    end

    cell_color(:,:,z,:)=reshape(make2DMask(im,celllabel,candidates),[dim(1) dim(2)*2+1 1 3]);

    cellnum=cellnum+totcell;    
    disp(num2str(z));

end

%%
write_colortiff_mex(fullfile(outdir,'cellmask.tif'),cell_color, int32(size(cell_color)));
save(fullfile(outdir,'cell_info.mat'),'cell_info');

