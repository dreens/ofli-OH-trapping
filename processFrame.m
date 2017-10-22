% Look at them all:
%last = broaden(r(1:200,1:250,15))/4;
%max(last(:))

h = zeros(size(r).*[2 2 1]-[1 1 0]);
h = permute(h,[1 2 4 3]);
for i=1:size(r,3)
    thisf = broaden(r(:,:,i));
    thisf(thisf<0) = 0;
    thisf = 0.5*thisf/mean(thisf(:));
    thisf(thisf>1) = 1;
    h(:,:,1,i) = thisf;
end

%% Write a Video
a = VideoWriter('watch0p5.avi');
a.FrameRate = 1;
open(a);
writeVideo(a,h(:,:,:,8:end));
close(a);