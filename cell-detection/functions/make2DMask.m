function [RGBstack]=make2DMask(im,mask,candidates)

    im=single(im);
    dim=size(im);

   %% immax=max(im(:));
   %% immin=min(im(:));

    RGBstack=zeros(dim(1),dim(2)*2+1,3,'uint8');
    im_norm=imNormalize95_candidates(im,candidates);%,[immin immax ]);
    im_norm255=uint8(255*im_norm);
    im_norm150=uint8(150*im_norm);
    
    RGBstack(:,1:dim(2),1)=im_norm255;
    RGBstack(:,1:dim(2),2)=im_norm255;
    RGBstack(:,1:dim(2),3)=im_norm255;
    
    RGBstack(:,dim(2)+2:dim(2)*2+1,1)=im_norm150;
    RGBstack(:,dim(2)+2:dim(2)*2+1,2)=im_norm150;
    RGBstack(:,dim(2)+2:dim(2)*2+1,3)=im_norm150;
    RGBstack(:,dim(2)+1,:)=255;
    
    tot_inds1=find(mask(:)>0);
    cellnum=mask(tot_inds1);

    c=ceil(tot_inds1/dim(1));
    r=mod(tot_inds1,dim(1));r(r==0)=dim(1);

    tot_inds2=zeros(length(r),3);

    tot_inds2(:,1)=dim(1)*dim(2)  +c*dim(1)+r;
    tot_inds2(:,2)=dim(1)*dim(2)*3+(c+1)*dim(1)+r;
    tot_inds2(:,3)=dim(1)*dim(2)*5+(c+2)*dim(1)+r;

    cindex=[122  0    0  ;
            255  0    0  ;
            209  102  0  ;
            163  204  0  ;
            82   230  51 ;
            0    255  102;
            0    178  178;
            0    102  255;
            0    51   255;
            204  0    255;
        ];

    cell_num= false(length(r),10);

    for i=1:10
        cell_num(:,i)=(mod(cellnum,10)==(i-1));
    end

    for i=1:3
        for j=1:10
           RGBstack(tot_inds2(cell_num(:,j),i))=cindex(j,i);
        end
    end

        


return;
    

