# ParticleMotionTools
Collection of software for seismic particle motion analysis with multi-wavelets.

Multi-wavelets are an idea originally devised by Lilly and Park (1995) (Geophysical J. International, 122, pp 1001-1021.)  The ideas are and extension of Slepian funtions that are the foundation of the so called multitaper spectral estimation methods or multiple window methods.   There is a large literature on multitaper methods, but comparatibly little on the multi-wavelet concept.  One reason for that is likely that there is no software base for people to experiment with the multi-wavelet concepts.  This library seeks to remedy that. 

This library is mostly C++ code. The exception is a directory called matlab that builds on Lilly and Park's original implementation to produce their form of wavelets and a variant they discuss but do not use that are essentially the product of slepian tapers and sin/cosine pairs.   The C++ libary is largely an OOP front end to a messy C library the author wrote a decade ago.   The C++ API is invitely cleaner build on concepts of the authors SEISPP library.   The highest level object in that library is an thing call a PMTimeSeries object.  PMTimeSeries abstract the idea of an equally spaced (in time) set of particle motion estimates specified in terms of an ellipse drawn in three-dimensions.

This library is largely useless in isolation.   It is intimately linked to the antelope package available to US academic institutions at no cost from http://www.brtt.com.   Potential users from institutions outside the US will need to make a special request to BRTT to obtain a license for the antelope software if they wish to use this package.   The primary exception is the matlab code which has no links to antelope whatsoever. 
