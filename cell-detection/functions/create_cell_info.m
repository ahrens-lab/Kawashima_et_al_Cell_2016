function cell_info=create_cell_info(cell_info,celllabel,totcell,z)

if isempty(cell_info)
    cell_info=struct('center',zeros(2,1),'area',0,'inds',0,'y_minmax',zeros(2,1),'x_minmax',zeros(2,1));
    c=1;
else
    c=0;
end

dim=size(celllabel);
tot_inds=find(celllabel>0); yx=ind2sub2d(dim,tot_inds);    
[cell_yx(:,1) cell_yx(:,2) cell_yx(:,3)]=calc_cell_yx_mex(int32(celllabel(tot_inds)),int32(totcell),double(yx),int32(length(tot_inds)));

cellnum=length(cell_info);
if z==1
    cellnum=0;
end

for i=1:totcell
    cell_info(cellnum+i).center=round(cell_yx(i,1:2));
    cell_info(cellnum+i).area  =round(cell_yx(i,3));
    cell_info(cellnum+i).slice =z;
end


rindex=celllabel(tot_inds);
[dummy ri]=sort(rindex,'ascend');
yx2=yx;
yx2(:,:)=yx(ri,:);
tot_inds2=tot_inds(ri);
range=0;

for i=1:totcell
    p=range+1:range+cell_yx(i,3);
    cell_info(cellnum+i).inds=tot_inds2(p);
    cell_info(cellnum+i).y_minmax(1)=min(yx2(p,1)); cell_info(cellnum+i).y_minmax(2)=max(yx2(p,1));
    cell_info(cellnum+i).x_minmax(1)=min(yx2(p,2)); cell_info(cellnum+i).x_minmax(2)=max(yx2(p,2));
    range=range+cell_yx(i,3);
end

%}