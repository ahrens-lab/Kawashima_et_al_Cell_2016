function [mdisk,  se, r, c, maskinds]=make_recog_disk(radius,dim)

mdisk=makeDisk(radius,radius*2+1);
se=strel(mdisk);
[r, c]=find(mdisk);
r=r-radius-1;
c=c-radius-1;
maskinds=r*dim(1)+c;
