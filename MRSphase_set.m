function [MRS_struct] = MRSphase_set(MRS_struct,spectno, phasezero, phasefirst)
% [MRS_struct] = MRSphase_set(MRS_struct,spectno, phasezero, phasefirst)
% CJE 18 Nov 2010
% Explicitly set zeroth and first order phase
% phasezero in degrees
% phasefirst in degrees per ppm.

%Zeroth order
deltaphase0 = phasezero -  MRS_struct.phase(spectno);
MRS_struct.gabaspec(spectno,:) = MRS_struct.gabaspec(spectno,:) *  exp(1i * deltaphase0 * pi / 180);
MRS_struct.phase(spectno)=phasezero;

% First order.
%centre on 17800, approx centre of gaba peak
% 5deg between 1 and 5ppm is 1degppm
% 5ppm -> pt16134 1ppm ->19473
deltaphase1 = phasefirst - MRS_struct.phase_firstorder(spectno);
deg_ppm = 5/(19473-16134);
deltapt = [1:32768] - 17800;
FirstOrderPhase = exp( 1i * deltaphase1 * deg_ppm * deltapt * pi / 180);

MRS_struct.gabaspec(spectno,:) = MRS_struct.gabaspec(spectno,:) .*  FirstOrderPhase;
MRS_struct.phase_firstorder(spectno) = phasefirst;

%figure(80)
%plot(MRS_struct.freq, real(MRS_struct.gabaspec(spectno,:)), 'k');
%set(gca,'XDir','reverse');
%oldaxis = axis;
%axis([0 5 0.5*oldaxis(3) oldaxis(4)])