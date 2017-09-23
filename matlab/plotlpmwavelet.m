function [ nw ] = plotlpmwavelet( w )
%plot routine to accompany lpmwavelet.m - plot the output as 
% subplots of complex pairs. 
% returns number of panels (wavelet pairs) as a pure conveniencep
[m,n]=size(w);
t=1:m;
nw=n/2;
for i=1:nw
    ii=2*(i-1)+1;
    subplot(nw,1,i),plot(t,w(:,ii),'b',t,w(:,ii+1),'r');
end

