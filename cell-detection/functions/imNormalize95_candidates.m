function out=imNormalize95_candidates(im,candidates)

if (~isempty(candidates) && length(candidates)>20)    
    temp=sort(double(im(candidates)),'descend');
    th1=temp(round(length(temp)/20));
else
    th1=double(max(im(:)));
end
temp=sort(double(im(:)),'ascend');
th2=double(temp(round(length(im(:))/1000)));

out=(im-th2)/(th1-th2);
out(out>1)=1;
out(out<0)=0;


