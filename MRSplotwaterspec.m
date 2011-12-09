function MRSplotwaterspec(MRS_struct)
%function MRSplotspec(MRS_struct)
% Plots MRS data loaded by LoadPfileGABA

%figure

plot(MRS_struct.freq, real(MRS_struct.waterspec));
set(gca,'XDir','reverse');
oldaxis = axis;
axis([4 5.5  0.5*oldaxis(3) oldaxis(4)])
