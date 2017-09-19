% Look at them all:
last = broaden(r(1:200,1:250,15))/4;
max(last(:))

h = zeros(399,499,1,20);
for i=1:20
    thisf = broaden(r(1:200,1:250,i));
    thisf(thisf<0) = 0;
    thisf = 0.25*thisf/mean(thisf(:));
    thisf(thisf>1) = 1;
    h(:,:,1,i) = thisf;
end

%% Write a Video
a = VideoWriter('watchprogress.avi');
a.FrameRate = 1;
open(a);
writeVideo(a,h(:,:,:,8:end));
close(a);