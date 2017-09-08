function [re,im] = slepianw0(npts,nw)
% comparable to slepianwavelt but for zero frequency.  
% This really just returns slepian functions for npts points
% nw is the number of complex pair wavetls returned, which 
% is passed as 2*nw to dpss.   
[tapers,v] = dpss(npts,2*nw);
ii=1;
for i=1:nw
    re(:,i)=tapers(:,ii);
    im(:,i)=tapers(:,ii+1);
    ii=ii+2;
end
%  slepianwavelet.m has a normalization here.  Not needed 
% because tapers are eigenvectors already normalized so
% L2 norm is 1


