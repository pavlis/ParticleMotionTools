function [re,im] = slepianwavelet(norder,ncycles,nw)
npts = 4*norder*ncycles;
[tapers,v] = dpss(npts,nw);
for i=1:npts
	omega = (i-1)*pi/(2*norder);
	rp(i) = cos(omega);
	ip(i) = sin(omega);
end
for i=1:2*nw
	re(:,i) = rp'.*tapers(:,i);
	im(:,i) = ip'.*tapers(:,i);
end

%  normalize so L2 norm is 1

for i=1:2*nw
	re(:,i) = re(:,i)/norm(re(:,i));
	im(:,i) = im(:,i)/norm(im(:,i));
end
