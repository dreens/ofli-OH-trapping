%% Partial Completion
% Make an OFLI image based on what has been stashed on the data drive of
% the cluster thus far
function ofli2 = partialImage(N,E,X,trap,P,varargin)

% Two different cases depending on whether the image was run with a version
% that supported different planes at that time.
datadir = '//data/ye/dare4983/splines/';
if isempty(varargin)
    thisdir = sprintf('N%d_E%d_X%d_%s_P%d',N,E,X,trap,P);
else
    thisdir = sprintf('N%d_E%.1f_X%d_%s_P%d_%s',N,E,X,trap,P,varargin{1});
end
assert(~~exist([datadir thisdir],'dir'),'No partial data found for those parameters');

% Two cases depending on whether point by point parfor had been implemented
% yet or not.
ds = dir([datadir thisdir '/dim*']);
dimfile = [datadir thisdir '/' ds(1).name];
if exist(dimfile,'file')
    load(dimfile,'N1','N2');
    fprintf('Image is %dx%d\n',N1,N2);
    ofli2 = zeros(N1*N2,200);
    fprintf('Getting Data from Files... %2d%%',0)
    count = 0;
    for i=1:N1*N2
        if ~mod(i,1000)
            fprintf('\b\b\b%2d%%',round(100*i/(N1*N2)))
        end
        dir2file = [datadir thisdir '/' num2str(mod(i,1000)) '/i=' num2str(i) '.mat'];
        if exist(dir2file,'file')
            count = count + 1;
            try
                tmp = load(dir2file);
                ofli2(i,:) = tmp.stash;
            catch
            end
        end
    end
    prognum = round(100*count/N1/N2);
    fprintf('\nPixel Progress:%d/%d=%2d%%\n',count,N1*N2,prognum);
    ofli2 = reshape(ofli2,[N1 N2 200]);
else
    ofli2 = zeros(3*N,N,20);
    maxi = 1;
    for i=1:3*N
        dir2file = [datadir thisdir '/i=' num2str(i) '.mat'];
        if exist(dir2file,'file')
            try
                tmp = load(dir2file);
                if N~=size(tmp.stash,1)
                    ofli2 = zeros(3*N,size(tmp.stash,1),20);
                end
                ofli2(i,:,:) = tmp.stash;
                maxi = i;
            catch
            end
        end
    end
    ofli2 = ofli2(1:maxi,:,:);
end

cut = 15;
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
if usejava('awt')
    imtool(image)
end
writeloc = sprintf('%s%s/panel%2d.png',datadir,thisdir,prognum);
imwrite(image,writeloc);

end
