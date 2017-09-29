%% Partial Completion
% Make an OFLI image based on what has been stashed on the data drive of
% the cluster thus far
function partialImage(N,E,X,trap,P)

datadir = '//data/ye/dare4983/splines/';
thisdir = sprintf('N%d_E%d_X%d_%s_P%d',N,E,X,trap,P);

assert(exist([datadir thisdir],'dir'),'No partial data found for those parameters');

ofli2 = zeros(N,N,20);

for i=1:N
    dir2file = [datadir thisdir '/i=' num2str(i) '.mat'];
    if exist(dir2file,'file')
        tmp = load(dir2file);
        ofli2(i,:,:) = tmp.stash;
    end
end

out = broaden(ofli2(:,:,end));
imtool(out)

end