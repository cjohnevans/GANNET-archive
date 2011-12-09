function MRSplotspec(MRS_struct)
%function MRSplotspec(MRS_struct)
% Plots MRS data loaded by MRSLoadPfiles

%figure

plot(MRS_struct.freq, real(MRS_struct.gabaspec_scaled));
legendtxt = regexprep(MRS_struct.pfile, '_','-');

legend(legendtxt)
set(gca,'XDir','reverse');
oldaxis = axis;
axis([0 5.5  0.5*oldaxis(3) oldaxis(4)])
