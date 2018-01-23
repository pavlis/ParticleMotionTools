function save_slepians_pf(re,im,nw,ncycles,fname)
%
%  Save slepian multiwavelets computed by parallel routine
%  to a pf file format usable by multiwavelet programs
%
% re - real part 
% im - imaginary part
% nw - nw parameter passed to dpss
% ncycle - ncycles of tapered harmonic of base wavelet
% fname - file name (character) to save results to
%
fid=fopen(fname,'w');
[n m] = size(re);
fprintf(fid,'nsamples    %d\n',n);
fprintf(fid,'nwavelets   %d\n',m);
%
%  These are computed for mwap and mwpm using a scaling by
%  1/fNyquist 
%
f0=2.0*ncycles/n;
fw = 2.0*nw/n;
fprintf(fid,'f0  %f\n',f0);
fprintf(fid,'fw  %f\n',fw);
fprintf(fid,'wavelets &Tbl{\n');
for j=1:m
	for i=1:n
		fprintf(fid,'%f %f\n',re(i,j),im(i,j));
	end
end
fprintf(fid,'}\n');
fclose(fid);