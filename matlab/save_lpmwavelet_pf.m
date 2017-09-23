function [ output_args ] = save_lpmwavelet_pf( u, fc, p, dt, fname )
%Creates a parameter file from multiwavlets stored in u.   It is
% assumed the wavelets are adjacentpairs in columns of u that can
% be placed in real and imaginary parts.   Always plots u as 
% a quality control.   
% Other inputs:
% fc - center frequency used for generating wavelets
% p - time-bandwidth product of these wavelets
% dt - sample interval used with fc in original calculation of wavelets
% (note - the point is fc, p, and dt should be the same as those sent
% to lpmwavelet.m script.)
% fname - output file name (string)
plotlpmwavelet(u);
title('multiwavelet pairs being saved');
figure;
[m,n]=size(u);
nw=n/2;
for i=1:nw
    ii=2*(i-1)+1;
    cw(:,i)=complex(u(:,ii),u(:,ii+1));
end
for i=1:nw
    subplot(2,1,1),plot(abs(cw(:,i)));
    title('time weighting');
    CW=fft(cw(:,i));
    subplot(2,1,2),plot(abs(CW));
    title('frequency weighting');
    pause
end
fid=fopen(fname,'w');
fprintf(fid,'nsamples    %d\n',m);
fprintf(fid,'nwavelets   %d\n',nw);
fNy=1.0/(2*dt);
f0=fc/fNy;
fw=p/(m*dt);
fwnd=fw/fNy;
fprintf(fid,'f0  %f\n',f0);
fprintf(fid,'fw  %f\n',fwnd);
fprintf(fid,'wavelets &Tbl{\n');
re=real(cw);
im=imag(cw);
for j=1:nw
	for i=1:m
		fprintf(fid,'%f %f\n',re(i,j),im(i,j));
	end
end
fprintf(fid,'}\n');
fclose(fid);
