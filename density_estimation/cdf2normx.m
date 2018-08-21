function [x]=cdf2normx(cdf)
x=norminv(cdf);
x(x>8)=8;
x(x<-8)=-8;
end