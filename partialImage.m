%% Partial Completion
% Make an OFLI image based on what has been stashed on the data drive of
% the cluster thus far
function ofli2 = partialImage(N,E,X,trap,P)

datadir = '//data/ye/dare4983/splines/';
thisdir = sprintf('N%d_E%d_X%d_%s_P%d',N,E,X,trap,P);

assert(~~exist([datadir thisdir],'dir'),'No partial data found for those parameters');

ofli2 = zeros(N,N,20);

for i=1:N
    dir2file = [datadir thisdir '/i=' num2str(i) '.mat'];
    if exist(dir2file,'file')
        tmp = load(dir2file);
        ofli2(i,:,:) = tmp.stash;
    end
end
cut = 10;
last = broaden(ofli2(:,:,end));
lost = (last==1);
last(last==2) = 20;
last(lost)=0;
red = last/20;
green = last/cut;
blue = 1-last/cut;
red(last<cut)=0;
green( last>cut | last<=0 ) = 0;
blue(  last>cut | last<=0 ) = 0;

red(lost) = 1;
blue(lost) = 0.5;
green(lost) = 0.5;
image = zeros([size(last) 3]);
image(:,:,1) = red;
image(:,:,2) = green;
image(:,:,3) = blue;
imtool(image)

end