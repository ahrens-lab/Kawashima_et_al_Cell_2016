function out=ind2sub2d(dim,inds)

out=zeros(length(inds),2);

out(:,2)=ceil(inds./dim(1));
out(:,1)=inds-(out(:,2)-1)*dim(1);