%transform a time series 
%
% function xn=red(x,varargin)
% 
% %default options
% opts=struct('boxcox','n',...    make a Box-Cox transfrmation
%             'detrend','n', ...   
%             'pad','y', ...      pad where NaN with mean of the serie
%             'anom',0);          seasonal detrend - compute anomalis for periodic series with period>0 
%         
% 
function [xt,sd,p1,H]=red2(x,varargin)

%default options
opts=struct('boxcox','n',...  
            'detrend','n', ...   
            'pad','n', ...   
            'anom',0,...
            'filter','n',...
            'hs',0);
        
opts=parseArgs(varargin,opts);

sz=size(x);
x1(:,1)=x;
N=length(x1);
if opts.detrend=='y' && opts.detrend >0 
y=smooth(x1,floor(N/24));
x=[1:N]';
 if opts.detrend=='y'
[fresult,gof,output] = fit(x(13:N-12),y(13:N-12),'poly1');
ci = confint(fresult,.99);
p1=fresult.p1;
p2=fresult.p2;
H=0;
if (ci(1,1)>=0 & ci(2,1)>=0) | (ci(1,1)<=0 & ci(2,1)<=0)
    x1=x1-p1*x-p2;
    H=1;
end
 elseif opts.detrend==2
     [fresult,gof,output] = fit(x(13:N-12),y(13:N-12),'poly2');

p1=fresult.p1;
p2=fresult.p2;
p3=fresult.p3;

    x1=x1-p1*x.^2-p2*x-p3;
 end
end

if opts.anom > 0
    if opts.filter=='y'
       x1(:,1) = filter_CWT(x1(:,1),'pp',[0 N/2]);
    end
x1(:,1) = ST(x1,opts.anom,opts.hs);
end

if opts.boxcox == 'y'
    if x1>0
        x1=boxcox(x1);
    else
        x1=boxcox(x1-min(x1)+1);
    end
end
sd=nanstd(x1);
if sum(x1)~=0 && sd>1e-6
x1=(x1-nanmean(x1))/sd;
end

if opts.pad == 'y'
    hh=find(isnan(x1)==1);
    x1(hh)=nanmean(x1);
end

xt(1:sz(1),1:sz(2))=x1;

%standarizzation of time series
function [ys,mY2] = ST(y,period,hs)

n1=length(y);
n=ceil(n1/period)*period;
y1=nan(n,1);
y1(1:n1,1)=y;

Y=reshape(y1,period,n/period);
mY=nanmedian(Y');
Y2=reshape(red(y1),period,n/period);
mY2=nanmedian(Y2');


if hs==1
    sY2=nanstd(Y2');
  ys=(y1'-repmat(mY,1,n/period))./repmat(sY2,1,n/period);
else
    ys=(y1'-repmat(mY,1,n/period));
end

ys=ys(1:n1);
end
end

    
    
