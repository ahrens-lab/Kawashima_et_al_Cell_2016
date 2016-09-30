% nmotion_estimation.m
%
% This scripts infers deformation of 3D volume through the timeseries.
% Detailes are described in Kawashima et al., 2016.
%
% Requires other funcitons in the "functions" folder.
% This script only runs on 64-bit Windows because it uses
% MEX files build by VS2010. For use in other systems please
% re-compile corresponding cpp files.
%
%
% Input files 
%
% StackDimensions.bin -> UINT 32 binary array of stack dimensions.
%                        [y dimension, xdimension, zdimension, timepoints]
% PlaneXX.stack       -> UINT16 timeseries for each Z-plane
%                        dimensions=[ydim, xdim, timepoints]
%
%
% Parameters
%
% input_dir         -> Data directory
% gridwidth         -> Specifies intervals of reference points and sizes of
%                      square patch around reference for calculating deformation
% tcycle            -> Cycle for averaging time series.
% xypixeldist       -> Pixel size in xy dimension in micron.
% zpixeldist        -> Interval between z-planes in micron.
% brightness_thre   -> Brightness threshold for setting reference points
%
%
% Output
%
% motion.tif        ->  Visual result of 3-dimensional deformation.
%                       XY and Z motions are plotted for each reference point.
%                       Left image shows xy movement by arrows scaled by 10.
%                       Right image shows z movement by color:
%                       [-2, -1, 0, 1, 2] z-plane movement is represented
%                       by [deep blue,cyan, green, yellow,red] dots.
%
% motion_graph.tif  ->  Quantitatition of 3-dimensional deformation
%                       Absolute XY and Z deformation values are
%                       averaged between reference points for each Z plane.
%                       Left image shows quantitation of xy movement.
%                       Right image shows quantitation of Z movement.
%
%
% Samples
%
% 'samples-xyrotation'  -> sample data with xy-rotational deformation
% 'samples-xzrotation'  -> sample data with xz-rotational deformation
%
%
% developed by Takashi Kawashima, HHMI Janelia Research Campus
% Sep 29,2016

clear all;close all;

input_dir = 'sample folder';
gridwidth=12; 
tcycle=3;
xypixeldist=0.406;
zpixeldist=5;
brightness_thre=130;

%% calculating deformations

motion_param=struct;
rwidth=gridwidth*2+1;
zmove=-2:2;

fname='StackDimensions.bin';
f=fopen(fullfile(input_dir,fname),'r');
dim=double(fread(f,'uint32'));
fclose(f);
refstack=zeros(dim(1),dim(2),dim(3));

for zz=1:dim(3)
    fname=['Plane',num2str(zz,'%.2d'),'.stack'];
    f=fopen(fullfile(input_dir,fname),'r');
    fstack=reshape(fread(f,dim(1)*dim(2)*tcycle,'uint16'),[dim(1) dim(2) tcycle]);
    fclose(f);
    refstack(:,:,zz) = mean(fstack,3);
end

for zz=1:dim(3)
    
    %% load files
    
    fname=['Plane',num2str(zz,'%.2d'),'.stack'];
    disp(['Calculating deformation of plane',num2str(zz)]);

    totcycle=floor(dim(4)/tcycle);
    stack=zeros(dim(1),dim(2),totcycle);
    
    f=fopen(fullfile(input_dir,fname),'r');
    fstack=reshape(fread(f,dim(1)*dim(2)*dim(4),'uint16'),[dim(1) dim(2) dim(4)]);
    fclose(f);
    
    for z=1:totcycle
        stack(:,:,z)=mean(fstack(:,:,((z-1)*tcycle+1):(z*tcycle)),3);
    end
    
    %% setting grid reference points
    
    startimg=double(stack(:,:,1));

    ygrid_num=floor((dim(1)-rwidth*2)/gridwidth)+1;
    xgrid_num=floor((dim(2)-rwidth*2)/gridwidth)+1;
    indslist=zeros(ygrid_num*xgrid_num,2);

    [r, c]=find(ones(rwidth)>0);
    r=r-gridwidth-1;c=c-gridwidth-1;
    inds=(c-1)*dim(1)+r;

    for x=1:xgrid_num
        for y=1:ygrid_num
            rinds=(rwidth+(x-1)*gridwidth)*dim(1)+rwidth+(y-1)*gridwidth+1;
            indslist(ygrid_num*(x-1)+y,1)=rinds;
            indslist(ygrid_num*(x-1)+y,2)=mean(startimg(rinds+inds));
        end
    end

    thre_inds=find(indslist(:,2)>brightness_thre);
    indslist2=indslist(thre_inds,1);
    
   %% calculating XY and Z motion for each reference points

    rs=zeros(length(indslist2),3,totcycle);  
    tilt=zeros(length(indslist2),3);    

    for i=1:length(indslist2)
        rinds=indslist2(i);
        source=fft2(reshape(startimg(rinds+inds),[rwidth rwidth]));

        for j=1:totcycle
            target=reshape(stack(rinds+inds+(j-1)*dim(1)*dim(2)),[rwidth rwidth]);

            buff=fftshift(ifft2(source.*conj(fft2(target))));
            [~,maxinds] = max(abs(buff(:)));

            a=-(mod(maxinds,rwidth)-gridwidth-1);
            if a==-gridwidth-1; a=gridwidth;end

            rs(i,1,j)=a; 
            rs(i,2,j)=-(ceil(maxinds/rwidth)-gridwidth-1);   

            moveinds=-rs(i,2,j)*dim(1)-rs(i,1,j);

            zlist=zmove+zz;
            zlist=zlist(zlist > 0 & zlist<=dim(3));
            tt=zeros(1,length(zlist));
            for k=1:length(zlist);     
                     rimg=refstack(:,:,zlist(k));
                     ztargetlist=rimg(rinds+moveinds+inds);                     
                     tt(k)=corr(target(:),ztargetlist(:));                   
            end
            
            [~,zmaxinds] = max(tt);
            rs(i,3,j)=zlist(zmaxinds)-zz;

        end
        
        py=polyfit((1:totcycle)',squeeze(rs(i,1,:)),1);
        px=polyfit((1:totcycle)',squeeze(rs(i,2,:)),1);
        pz=polyfit((1:totcycle)',squeeze(rs(i,3,:)),1);

        tilt(i,1)=py(1)*(totcycle-1);
        tilt(i,2)=px(1)*(totcycle-1);
        tilt(i,3)=pz(1)*(totcycle-1);

    end

    %%  applying median filter between reference points
    
    [r, c]=find(ones(3,3));
    median_inds=(c-2)*ygrid_num+(r-2);
    tilt_med=zeros(length(indslist2),3); 

    for jj=1:3

        move_matrix=zeros(ygrid_num,xgrid_num);
        ones_matrix=zeros(ygrid_num,xgrid_num);

        for ii=1:length(indslist2)
            move_matrix(thre_inds(ii))=tilt(ii,jj);
            ones_matrix(thre_inds(ii))=1;
        end

        for ii=1:length(indslist2)
            tmp = thre_inds(ii)+median_inds;
            inds1=thre_inds(ii)+median_inds(tmp>0 & tmp<=xgrid_num*ygrid_num);
            tilt_med(ii,jj)=median(move_matrix(inds1(ones_matrix(inds1)>0)));
        end
    end
    
    img=repmat(imNormalize99(startimg),[1 1 3]);
    
    motion_param(zz).tilt_med=tilt_med;
    motion_param(zz).indslist2=indslist2;
    motion_param(zz).img=img;
    motion_param(zz).xymove_av=mean(sqrt(tilt_med(:,1).^2 + tilt_med(:,2).^2))*xypixeldist;
    motion_param(zz).xymove_sd=std(sqrt(tilt_med(:,1).^2 + tilt_med(:,2).^2)*xypixeldist,[],1);
    motion_param(zz).zmove_av=mean(abs(tilt_med(:,3)))*zpixeldist;
    motion_param(zz).zmove_sd=std(abs(tilt_med(:,3))*zpixeldist);

end


%%  plot deformation visually

f1=figure(1);set(f1,'Position',[150 150 400 320]);

colorlist=[0 0 1;
           0 1 1;
           0 1 0;
           1 1 0;
           1 0 0;];

for zz=1:dim(3)
    clf(f1);

    indslist2=motion_param(zz).indslist2;
    tilt=motion_param(zz).tilt_med;
    img=motion_param(zz).img;
    
    ha=axes; set(ha,'Position',[0 0 1 1]);
    axis off;

    image(img);hold on;    
    hq=quiver(ceil(indslist2/dim(1)),mod(indslist2,dim(1)),tilt(:,2),tilt(:,1),0,'linewidth',1,'Color',[1 0 0]);
    hold off;
    hU = get(hq,'UData');hV = get(hq,'VData') ;
    set(hq,'UData',10*hU,'VData',10*hV)

    CC=getframe(f1);    
    motion_param(zz).arrows=CC.cdata;
    cla;

    set(ha,'Position',[0 0 1 1]);
    image(img);hold on;

    for i=1:5
        inds=find(round(tilt(:,3))==zmove(i));
        if ~isempty(inds)
            for j=1:length(inds)
                rectangle('Position',[ceil(indslist2(inds(j))/dim(1))-1,mod(indslist2(inds(j)),dim(1))-1, 3, 3],'FaceColor',colorlist(i,:));
            end
        end
    end

    hold off;

    CC=getframe(f1);    
    motion_param(zz).zmotion=CC.cdata;

end

odim=size(motion_param(1).arrows);
outstack=zeros(odim(1),odim(2)*2+1,dim(3),3);

for i=1:dim(3)
    outstack(:,1:odim(2),i,:)=reshape(motion_param(i).arrows,[odim(1) odim(2) 1 3]);
    outstack(:,(1:odim(2))+odim(2)+1,i,:)=reshape(motion_param(i).zmotion,[odim(1) odim(2) 1 3]);
end

write_colortiff_mex(fullfile(input_dir,'motion.tif'),uint8(outstack), int32(size(outstack)));

%% plot overall quantitation

f2=figure(2);
set(f2,'Position',[100 100 700 300]);
subplot(1,2,1);
for z=1:dim(3)
    indslist=motion_param(z).indslist2;
    tilt=motion_param(z).tilt_med;
    plot(ones(1,length(motion_param(z).indslist2))*z,sqrt(tilt(:,1).^2+tilt(:,2).^2)*xypixeldist,'.');hold on;
end
errorbar([motion_param.xymove_av],[motion_param.xymove_sd],'r','LineWidth',2,'LineStyle','none');
scatter(1:dim(3),[motion_param.xymove_av],'ro','fill');hold off;
xlim([0 dim(3)+1]);ylim([-0.5 2]);title('XY motion');
xlabel('Z-plane');ylabel('Mean deviation (micron)');

subplot(1,2,2);
for z=1:dim(3)
    indslist=motion_param(z).indslist2;
    tilt=motion_param(z).tilt_med;
    plot(ones(1,length(motion_param(z).indslist2))*z,abs(tilt(:,3))*zpixeldist,'.');hold on;
end
errorbar([motion_param.zmove_av],[motion_param.zmove_sd],'r','LineWidth',2,'LineStyle', 'none');
scatter(1:dim(3),[motion_param.zmove_av],'ro','fill');hold off;
xlim([0 dim(3)+1]);ylim([-2 8]);title('Z motion');xlabel('Z-plane')
xlabel('Z-plane');ylabel('Mean deviation (micron)');

set(f2,'PaperPositionMode','auto'); 
saveas(f2,fullfile(input_dir,['motion_graph.tif']),'tif');
