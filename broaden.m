%% Broadens a matrix into a 2x2.
% Useful for taking advantage of symmetries.
function bmat = broaden(mat)

bmat = [mat(end:-1:2,end:-1:2) mat(end:-1:2,:) ; ...
        mat(:,end:-1:2) mat];
    
end
