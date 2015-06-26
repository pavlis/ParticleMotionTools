function makewavelet(norder,ncycles,nw,fname)
%
% driver for slepianwavelet to run and save all at once
%
%  norder - order = degree of oversampling (1 is one point per cycle)
%  ncycles - number of cycles of harmonic component
%  nw - nw parameter for dpss (related to time-bandwidth product)
%  fname - output file name
%
[re im] = slepianwavelet(norder,ncycles,nw);
save_slepians_pf(re,im,nw,ncycles,fname);
