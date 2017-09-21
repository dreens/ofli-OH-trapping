function [vals, terms, dirs] = escbase(~,y,fff) 
    y2 = y(7:12);
    y3 = y(13:18);
    fy = fff(0,y);
    fy = fy(1:6);
    fl = (y2+0.5*y3);
    pj = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;
    ofl = fl - pj(fl,fy);
    ofl = log10(sqrt(sum(ofl.^2)));
    
    val1 = max(abs(y(1:6)))-10;
    val2 = ofl-20;
    vals = [val1 val2];
    terms = [1 1];
    dirs = [0 0];
end