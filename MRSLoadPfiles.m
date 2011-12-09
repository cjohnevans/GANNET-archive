function [MRS_struct] = MRSLoadPfiles110426(temppfile, FreqPhaseAlign)
%function [MRS_struct] = MRSLoadPfiles110426(pfiles, FreqPhaseAlign)
% CJE110310.  Loads P files for analysis
%
%inputs:
% pfiles     = cell array of P filenames { 'pfile1.7' 'pfile2.7' ... }
% FreqPhaseAlign  = Perform frame-by-frame frequency and phase alignment of spectra and
% Reject pairs of frames based if phase or freq > 3*stdev
%              0=NO, 1=YES


% 091201: Dec 09.  Initial version, based on Richard's code
% 100301: March 10.  
%         Phasing of water data - use first point for each average
%          independently to reduce phase errors.
%         Return a struct, rather than set of arrays
% 101001: Oct 2010, named MRSdatacheck - based on MRSLoadPfiles 100301
%         Frequency realign each transient (water peak) - store as
%          DiffFIDrealign.  Note waterspec, gabaspec, and sumspec
%          are the unaligned versions
%         Generate plots for frequecy shifts, phase variation
%         Fix timeseries swap bug
% 101101: Phase correct each transient first point in FID
% 6 Nov 2010: ###########problem with water fitting############
% 110208: Renamed as MRSLoadPfiles110208
%         Save realigned data by default to gabaspec etc...
% 110303: Make frame-by-frame frequency and phase correction
%          optional - MRS data AND water scans
%         OutlierReject.
% 110310: ** BETA VERSION ON /CUBRIC/SCRATCH IN APRIL 2011 **
%         LB = 2Hz for this version
%         Fit Creatine resonance for each pair of shots (Wadell,2007)
%           for frequency and phase correction
% 110426: ** RELEASED ON 2 AUG 2011 (MRS WORKSHOP) **
%         Return Cr fit info for internal reference
%         Change naming of 'amplitude' where it means 'area'
%         LB back to 4Hz
%         Plot Water freq in Hz

%global freq;
%global filename;
global Comft0;
global Comft1;
global ComftW;
%global pfile;

%LB=4; %RE default
LB = 4; 

MRS_struct.pfile=temppfile;
MRS_struct.LB=LB;
MRS_struct.versionload = '110426';
disp(['Load Version is ' MRS_struct.versionload ]);


global_rescale = 1/1e17;  % rescale all loaded data to get a sensible range

if iscell(temppfile) == 1 % it's a cell array, so work out the number of elements
  numpfiles=numel(temppfile);
  pfiles=temppfile;
else
  numpfiles=1;  % it's just one pfile
  pfiles{1}=temppfile;
end

for ii=1:numpfiles
  %    MRS_struct.pfile{ii}
  %LB=2;
  split = 1;
  i = sqrt(-1);
  %  MRS_struct.phase(ii,:)=0;
  
  %disp('opening p-file %s', MRS_struct.pfile{ii})
  % Open Pfile to read reference scan data.
  fid = fopen(pfiles{ii},'r', 'ieee-be');

  if fid == -1
    tmp = [ 'Unable to locate Pfile ' pfiles{ii} ];
    disp(tmp);
    return;
  end
  % return error message if unable to read file type.

  % create dir for output
  if(exist('MRSload_output','dir') ~= 7)
    mkdir MRSload_output
  end
  

  % Determine size of Pfile header based on Rev number
  status = fseek(fid, 0, 'bof');
  [f_hdr_value, count] = fread(fid, 1, 'real*4');
  rdbm_rev_num = f_hdr_value(1);
  if( rdbm_rev_num == 7.0 );
    pfile_header_size = 39984;  % LX
  elseif ( rdbm_rev_num == 8.0 );
    pfile_header_size = 60464;  % Cardiac / MGD
  elseif (( rdbm_rev_num > 5.0 ) && (rdbm_rev_num < 6.0));
    pfile_header_size = 39940;  % Signa 5.5
  else
    % In 11.0 and later the header and data are stored as little-endian
    fclose(fid);
    fid = fopen(pfiles{ii},'r', 'ieee-le');
    status = fseek(fid, 0, 'bof');
    [f_hdr_value, count] = fread(fid, 1, 'real*4');
    if (f_hdr_value == 9.0)  % 11.0 product release
      pfile_header_size= 61464;
    elseif (f_hdr_value == 11.0);  % 12.0 product release
      pfile_header_size= 66072;
    elseif (f_hdr_value > 11.0) & (f_hdr_value < 100.0)  % 14.0 and later
      status = fseek(fid, 1468, 'bof');
      pfile_header_size = fread(fid,1,'integer*4');
    else
      err_msg = sprintf('Invalid Pfile header revision: %f', f_hdr_value );
      return;
    end
  end

  % Read header information
  status = fseek(fid, 0, 'bof');
  [hdr_value, count] = fread(fid, 102, 'integer*2');
  npasses = hdr_value(33);
  nslices = hdr_value(35);
  nechoes = hdr_value(36);
  nframes = hdr_value(38);
  point_size = hdr_value(42);
  da_xres = hdr_value(52);
  da_yres = hdr_value(53);
  rc_xres = hdr_value(54);
  rc_yres = hdr_value(55);
  start_recv = hdr_value(101);
  stop_recv = hdr_value(102);
  nreceivers = (stop_recv - start_recv) + 1;
  
  
  % Specto Prescan pfiles
  if (da_xres == 1) & (da_yres == 1);
    da_xres = 2048;
  end

  % Determine number of slices in this Pfile:  this does not work for all cases.
  slices_in_pass = nslices/npasses;

  % Compute size (in bytes) of each frame, echo and slice
  data_elements = da_xres*2;
  frame_size = data_elements*point_size;
  echo_size = frame_size*da_yres;
  slice_size = echo_size*nechoes;
  mslice_size = slice_size*slices_in_pass;
  my_slice = 1;
  my_echo = 1;
  my_frame = 1;

  FullData=zeros(nreceivers, da_xres , da_yres-my_frame+1);
  %ChannelSplit=zeros(nreceivers, da_xres,split);

  %Start to read data into Eightchannel structure.
  totalframes=da_yres-my_frame+1;
  chunk=uint16(totalframes/split);

  data_elements2 = data_elements*totalframes*nreceivers;

  %  % Compute offset in bytes to start of frame.
  file_offset = pfile_header_size + ((my_frame-1)*frame_size);

  status = fseek(fid, file_offset, 'bof');

  % read data: point_size = 2 means 16 bit data, point_size = 4 means EDR )
  if (point_size == 2 )
    [raw_data, count] = fread(fid, data_elements2, 'integer*2');
  else
    [raw_data, count] = fread(fid, data_elements2, 'integer*4');
  end

  fclose(fid);
  
  % get some info from the header
  % pfil_hdr = [pfiles{ii} '.hdr' ]
  
  % fid = fopen(pfil_hdr,'r');
  
  % if fid == -1
  %   % no header, so try to run print_raw_headers
  %   stat=system('which print_raw_headers') % see if print_raw_headers is available 
  %   if(stat)
  %     tmp = [ 'print_raw_headers ' pfiles{ii}
  %     system(
  % end

  
  % 110303 CJE
  % calculate Navg from nframes, 8 water frames, 2 phase cycles
  % Needs to be specific to single experiment - for frame rejection
  MRS_struct.Navg(ii) = (nframes-8)*2; 
  
  ShapeData = reshape(raw_data,[2 da_xres totalframes nreceivers]);
  ZeroData = ShapeData(:,:,1,:);
  WaterData = ShapeData(:,:,2:9,:);
  FullData = ShapeData(:,:,10:end,:);

  totalframes = totalframes-9;

  Frames_for_Water = 8;

  FullData = FullData.*repmat([1;i],[1 da_xres totalframes nreceivers]);
  WaterData = WaterData.*repmat([1;i],[1 da_xres Frames_for_Water nreceivers]);

  FullData = squeeze(sum(FullData,1));
  FullData = permute(FullData,[3 1 2]);
  
  
  WaterData = squeeze(sum(WaterData,1));
  WaterData = permute(WaterData,[3 1 2]);
  % at this point, FullData(rx_channel, point, average)
  
  firstpoint=conj(WaterData(:,1,:));
  firstpoint=repmat(firstpoint, [1 da_xres 1]);
  % here firstpoint(rx_channel,[], average)

  
  % CJE March 10 - correct phase of each Water avg independently
  WaterData=WaterData.*firstpoint*global_rescale;;

  %Multiply the Eightchannel data by the firstpointvector 
  % zeroth order phasing of spectra
  % CJE Nov 09: do global rescaling here too
  % don't really need the phasing step here if performing frame-by-frame phasing
  for receiverloop = 1:nreceivers
    FullData(receiverloop,:) = FullData(receiverloop,:)*firstpoint(receiverloop,1,1)*global_rescale;
    % WaterData(receiverloop,:) = WaterData(receiverloop,:)*firstpoint(receiverloop,1,1)*global_rescale;
  end

  % sum over Rx channels
  FullData = squeeze(sum(FullData,1));
  WaterData = squeeze(sum(WaterData,1));

  %Water = WaterData;

  time=(1:1:size(FullData,1))/2000;
  time_zeropad=(1:1:8*size(FullData,1))/16000;
  DataSize = size(FullData,2);
  
  % Finish processing water data.
  ComWater = sum(WaterData,2);
  ComWater = ComWater.*exp(-(time')*LB);
  MRS_struct.waterspec(ii,:)=fftshift(fft(ComWater,size(ComWater,1)*8,1))';

  % CJE Oct 10
  % FT each transient (Zerofill 8 times)
  AllFramesFT=fftshift(fft(FullData,size(FullData,1)*8,1),1);
  
  % work out frequency scale
  Csize=size(AllFramesFT,1);
  freqrange=5000/127.35;
  MRS_struct.freq=(Csize+1-(1:1:Csize))/Csize*freqrange+4.7-freqrange/2.0;

  %%%%%%%%%  Frame-by-frame Frequency Realignment %%%%%%%%
  % find peak location for frequency realignment
  [FrameMax, FrameMaxPos] = max(AllFramesFT, [], 1);
  % align all peaks to _first_ transient (should be closest value set during Prescan)
  FrameShift = FrameMaxPos - FrameMaxPos(1);
  % if(FreqPhaseAlign)       % CJE 110303
  %   for jj=1:totalframes
  % 	 AllFramesFT(:,jj)=circshift(AllFramesFT(:,jj), -FrameShift(jj));
  %   end
  % end
  
  % plot waterfreq in Hz, with f=0 at first point.
  MRS_struct.waterfreq(ii,:) = MRS_struct.freq(FrameMaxPos)*42.58*3;
  MRS_struct.waterfreq(ii,:) = MRS_struct.waterfreq(ii,:) - ...
      MRS_struct.waterfreq(ii,1);    
  
  % CJE Oct10: Calculate sum and difference spectra - no filtering at this point...
  %     OddFramesFTrealign=AllFramesFT(:,1:2:end);
  %     EvenFramesFTrealign=AllFramesFT(:,2:2:end);
  OddFramesFT=AllFramesFT(:,1:2:end);
  EvenFramesFT=AllFramesFT(:,2:2:end);
  
  

    % Do Creatine sum fit in all cases
    lb = 17651; ub=18000;
    freqrange = MRS_struct.freq(lb:ub);
    % initx = [area hwhm f0 phase baseline0 baseline1 ]
    % n.b. initials are in [ area ppm ppm rads baseline0 baseline1 ], but returned values are
    % [area FWHM_Hz Hz deg baseline0 baseline1] 
    Cr_initx = [ 30 0.1 3.0 0 0 0 ];

    %Come up with better initial conditions by fitting to the Cr sum
    CrSumSpec = sum(AllFramesFT(lb:ub,:),2);
    CrSumSpecFit = FitPeaksByFrames(freqrange, CrSumSpec, Cr_initx);
    Cr_initx = CrSumSpecFit;
    %rescale to account for Navg
    Cr_initx(1) = Cr_initx(1)/size(AllFramesFT,2);
    Cr_initx(5) = 0;
    %       Cr_initx(6) = 0;
    
    MRS_struct.Cr_area(ii) = Cr_initx(1);
    MRS_struct.Cr_freq(ii) = Cr_initx(3);
    MRS_struct.fwhmHz(ii) = Cr_initx(2);

  % 110310 Freq and Phase alignment by Cr fitting
  if(FreqPhaseAlign)
    
    %need to rescale area for 1 frame bit messy -
    % convert back to  ppm, rads for new initialiser 
    conv = [1 1/(2*42.576*3) 1/(42.576*3) (pi/180) 1 1 ];
    Cr_initx = Cr_initx .* conv;

    % Water peak freq shift estimation
    % find peak location for frequency realignment
    [FrameMax, FrameMaxPos] = max(AllFramesFT, [], 1);
    % align all peaks to _first_ transient (should be closest value set during Prescan)
    FrameShift = FrameMaxPos - FrameMaxPos(1);
    
    [CrFitParams, rejectframes] = FitPeaksByFrames(freqrange, AllFramesFT(lb:ub,:), Cr_initx);
    CrFreqShift = CrFitParams(:,3);
    CrFreqShift = CrFreqShift - mean(CrFreqShift);
    CrFreqShift = CrFreqShift ./ (3*42.58*(MRS_struct.freq(2) - MRS_struct.freq(1) ));
    CrFreqShift_points = round(CrFreqShift);
    
    % average over ON and OFF spectra - otherwise there is a net freq shift of ON relative to 
    %  Off, causing a subtraction error 110310 CJE
    CrFreqShift_avg = reshape(CrFreqShift, [2 (numel(CrFreqShift)/2) ]);
    CrFreqShift_avg = mean(CrFreqShift_avg, 1);
    CrFreqShift_avg = repmat(CrFreqShift_avg, [2 1]);
    CrFreqShift_avg = round(reshape(CrFreqShift_avg, [ (numel(CrFreqShift_avg)) 1 ]));
    CrPhaseShift = CrFitParams(:,4);
    
    for jj=1:totalframes
      AllFramesFTrealign(:,jj)=AllFramesFT(:,jj) * exp(1i * -CrPhaseShift(jj) * pi /180);
      AllFramesFTrealign(:,jj)=circshift(AllFramesFT(:,jj), -CrFreqShift_avg(jj)); %Cr peak realignment
%	 AllFramesFTrealign(:,jj)=circshift(AllFramesFT(:,jj), -CrFreqShift_points(jj)); %Cr peak realignment
%	 AllFramesFTrealign(:,jj)=circshift(AllFramesFT(:,jj), -FrameShift(jj)); %Water peak realignment
    end
    
    % Need to recalculate these from the f, phase corrected versions...
    OddFramesFTrealign=AllFramesFTrealign(:,1:2:end);
    EvenFramesFTrealign=AllFramesFTrealign(:,2:2:end);

    % CJE 110303: OutlierReject
    % reject based on water fit (not Cr)
    %	 frequpper = mean(MRS_struct.waterfreq(ii,:)) + 3*std(MRS_struct.waterfreq(ii,:));
    %	 freqlower = mean(MRS_struct.waterfreq(ii,:)) - 3*std(MRS_struct.waterfreq(ii,:));
    %	 size(MRS_struct.waterfreq(ii,:))
    %	 frequpper = repmat(frequpper, [1 totalframes]);
    %	 freqlower = repmat(freqlower, [1 totalframes]);
    
    lastreject = -1;
    numreject=0;
    for jj=1:totalframes
      if rejectframes(jj)
	% work out the ON_OFF pair to reject
	% set to zero - will come out in the wash after SUM, DIFF 
	pairnumber = round(jj/2);
	if(pairnumber ~= lastreject) %trap condition where ON and OFF in a pair are rejected
	  OddFramesFTrealign(1:end, pairnumber) = 0;
	  EvenFramesFTrealign(1:end, pairnumber) = 0;
	  %	       AllFramesFTrealign(1:end, (2*pairnumber-1):(2*pairnumber)) = 0;
	  lastreject = pairnumber;
	  numreject = numreject + 2;
	end
      end	   
    end

    for jj=1:(totalframes/2)
      AllFramesFTrealign(:,(2*jj-1)) = OddFramesFTrealign(1:end, jj);
      AllFramesFTrealign(:,(2*jj)) =  EvenFramesFTrealign(1:end, jj); 
    end
    
    MRS_struct.Navg(ii) = MRS_struct.Navg(ii) - numreject;
  else 
    % no realignment
    AllFramesFTrealign=AllFramesFT;
    OddFramesFTrealign=AllFramesFT(:,1:2:end);
    EvenFramesFTrealign=AllFramesFT(:,2:2:end);
    numreject = -1;
  end
  
  DiffSpecRealign(ii,:)=sum(OddFramesFTrealign,2)-sum(EvenFramesFTrealign,2);
  SumSpecRealign(ii,:)=sum(OddFramesFTrealign,2)+sum(EvenFramesFTrealign,2);
  % Do the line broadening on the realigned data
  DiffFIDrealign=ifft(fftshift(DiffSpecRealign(ii,:)));
  SumFIDrealign=ifft(fftshift(SumSpecRealign(ii,:)));
  [B,A] = butter(4, 0.02, 'high');
  %DiffFIDrealign = filtfilt(B,A,DiffFIDrealign);
  %SumFIDrealign = filtfilt(B,A,SumFIDrealign);
  DiffFIDrealign = DiffFIDrealign.*exp(-(time_zeropad)*8*LB); %has already been zeropadded..
  SumFIDrealign = SumFIDrealign.*exp(-(time_zeropad)*8*LB);
  DiffSpecRealign(ii,:)=fftshift(fft(DiffFIDrealign));
  SumSpecRealign(ii,:)=fftshift(fft(SumFIDrealign));
  
  % ditto on non-aligned data (CJE 30/10/11)
  DiffSpec(ii,:)=sum(OddFramesFT,2)-sum(EvenFramesFT,2);
  SumSpec(ii,:)=sum(OddFramesFT,2)+sum(EvenFramesFT,2);     
  DiffFID=ifft(fftshift(DiffSpec(ii,:)));
  SumFID=ifft(fftshift(SumSpec(ii,:)));
  DiffFID = DiffFID.*exp(-(time_zeropad)*8*LB); %has already been zeropadded..
  SumFID = SumFID.*exp(-(time_zeropad)*8*LB);
  DiffSpec(ii,:)=fftshift(fft(DiffFID));
  SumSpec(ii,:)=fftshift(fft(SumFID));

  %  CJE Feb 11:  Return the realigned spectra
  MRS_struct.gabanoalign(ii,:)=DiffSpec(ii,:);
  MRS_struct.gabaspec(ii,:)=DiffSpecRealign(ii,:);
  MRS_struct.sumspec(ii,:)=SumSpecRealign(ii,:);
  % calculate FWHM, based on NAA peak 
  %MRS_struct.fwhmHz(ii)=-fwhm(MRS_struct.freq(18300:18850), ...
%			      real(MRS_struct.sumspec(ii,18300:18850)))*42.58*3;
  %     MRS_struct.PhaseStdevDeg(ii) = std(framephase);
  MRS_struct.FreqStdevHz(ii) = std(MRS_struct.waterfreq(ii,:));
  MRS_struct.phase(ii) = 0; % initial zeroth order phase
  MRS_struct.phase_firstorder(ii) = 0; % initial 1st order phase
  MRS_struct.FreqPhaseAlign = FreqPhaseAlign; %frame-by-frame f align
  %MRS_struct.EvenFrames = EvenFramesFTrealign;
  MRS_struct.Rejects(ii) = numreject;
  
  if(ishandle(101))
    close(101)
  end
  figure(101)

  subplot(2,2,1)
  MRSplotprepostalign(MRS_struct,ii)
  %figure(53)
  subplot(2,2,2)
  if(FreqPhaseAlign)
    rejectframesplot = (1./rejectframes) .*  MRS_struct.waterfreq(ii,:)';
    plot([1:DataSize], MRS_struct.waterfreq(ii,:)', '-', [1:DataSize], rejectframesplot, 'ro')
  else
    plot([1:DataSize], MRS_struct.waterfreq(ii,:)')
  end
  xlabel('frame'); ylabel('f_0')
  title('Water Frequency, Hz');
  
  subplot(2,2,3)

  if(FreqPhaseAlign)     
    plotrealign=[ real(AllFramesFT((lb+50):(ub-150),:)) ; 
		  real(AllFramesFTrealign((lb+50):(ub-150),:)) ];
    imagesc(plotrealign);
  else
    tmp = 'No realignment';
    text(0,0.9, tmp, 'FontName', 'Courier');
  end

  subplot(2,2,4);
  axis off;
  
  tmp = [ 'pfile       : ' MRS_struct.pfile{ii} ];
  tmp = regexprep(tmp, '_','-');
  text(0,0.9, tmp, 'FontName', 'Courier');
  tmp = [ 'Navg        : ' num2str(MRS_struct.Navg(ii)) ];
  text(0,0.8, tmp, 'FontName', 'Courier');
  tmp = sprintf('FWHM (Hz)   : %.2f', MRS_struct.fwhmHz(ii) );
  text(0,0.7, tmp, 'FontName', 'Courier');
  tmp = sprintf('FreqSTD (Hz): %.2f', MRS_struct.FreqStdevHz(ii));
  text(0,0.6, tmp, 'FontName', 'Courier');
  tmp = [ 'LB (Hz)     : ' num2str(MRS_struct.LB,1) ];
  text(0,0.5, tmp, 'FontName', 'Courier');
  %tmp = [ 'Align/Reject: ' num2str(MRS_struct.FreqPhaseAlign) ];
  %text(0,0.5, tmp, 'FontName', 'Courier');
  tmp = [ 'Rejects     : '  num2str(MRS_struct.Rejects(ii)) ];
  text(0,0.4, tmp, 'FontName', 'Courier');
  tmp = [ 'LoadVer     : ' MRS_struct.versionload ];
  text(0,0.3, tmp, 'FontName', 'Courier');
  
  
  pfil_nopath = MRS_struct.pfile{ii};
%  pfil_nopath = pfil_nopath( (length(pfil_nopath)-15) : (length(pfil_nopath)-9) );
  tmp = strfind(pfil_nopath,'/');
  tmp2 = strfind(pfil_nopath,'\');
  if(tmp)
    lastslash=tmp(end);
  elseif (tmp2) 
    %maybe it's Windows...
      lastslash=tmp2(end);
  else
    % it's in the current dir...  
      lastslash=0;
  end
  tmp = strfind(pfil_nopath, '.7');
  dot7 = tmp(end); % just in case there's another .7 somewhere else...
  pfil_nopath = pfil_nopath( (lastslash+1) : (dot7-1) );

 epsfname=[ 'MRSload_output/' pfil_nopath  '.eps' ];
  %   saveas(100+ii, epsfname, 'psc2')
  saveas(101, epsfname, 'psc2');
  
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%% end of main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 110310 CJE fit sum pairs to get freq and phase information of frames
function [FitParams, rejectframe]  = FitPeaksByFrames(freq, FrameData, initx)

%options = optimset('lsqcurvefit');
%options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5 , ...
%		   'MaxFunEvals', 1e12);
% initx = [area hwhm f0 phase baseline0 baseline1] 

nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts, 'MaxIter', 1e5, 'Display','Off');
nframes = size(FrameData,2);

for jj = 1:nframes
  % [fit_param, resnorm, resid, exitflag ]  = ...
  %     lsqcurvefit(@(xdummy,ydummy) LorentzModel(xdummy, ydummy), initx, ...
  % 		  freq', real(FrameData(:,jj)));
  fit_param = nlinfit(freq', real(FrameData(:,jj)), ...
		      @(xdummy, ydummy) LorentzModel(xdummy, ydummy), ...
		      initx, nlinopts);
  FitParams(jj,:) = fit_param;
  fit_plot = LorentzModel(fit_param, freq);
  %  figure(3); plot(freq', real(FrameData(:,jj)), 'g', freq', fit_plot,'b');
  %set(gca,'XDir','reverse');
  %  input('next')
end

% Need to deal with phase wrap:
% Convert to complex number then recalculate phase within 2*pi range
phase_wrapped = FitParams(:,4);
cmplx = cos(phase_wrapped) + 1i * sin(phase_wrapped);
phase_unwrapped = angle(cmplx);

% then fix to be within -pi..pi
offsetpos =  pi*lt(phase_unwrapped, -pi/2);
offsetneg = -pi*gt(phase_unwrapped,  pi/2);
phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
FitParams(:,4) = phase_unwrapped;

% Fix area and linewidth to be positive
FitParams(:,1) = abs(FitParams(:,1));
FitParams(:,2) = abs(FitParams(:,2));

% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = repmat([1 (2*42.576*3) (42.576*3) (180/pi) 1 1 ], [nframes 1]);

FitParams = FitParams .* conv;

% Reject any point where the fit params - area, fwhm, phase
%  or freq are > 3stdev away from the mean
% set reject criteria for all fit parameters
MeanFitParams = mean(FitParams, 1);
UpperLim = repmat(MeanFitParams + 3*std(FitParams,1), [nframes 1]);
LowerLim = repmat(MeanFitParams - 3*std(FitParams,1), [nframes 1]);
%but don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectframe = gt(FitParams, UpperLim);
rejectframe = rejectframe + lt(FitParams, LowerLim);
rejectframe = max(rejectframe,[],2); 

% figure(102); 
% subplot(2,2,1)
% rejectplot=(1./rejectframe).*FitParams(:,1); % odd syntax, but only
%                                             % plot the rejected points!
% plot([1:nframes], FitParams(:,1),'-', [1:nframes], rejectplot,'or'); title('amp');
% subplot(2,2,2)
% rejectplot=(1./rejectframe).*FitParams(:,2);
% plot([1:nframes], FitParams(:,2),'-', [1:nframes], rejectplot,'or'); title('fwhm')

% subplot(2,2,3)
% rejectplot=(1./rejectframe).*FitParams(:,3);
% plot([1:nframes], FitParams(:,3),'-', [1:nframes], rejectplot,'or'); title('freq')

% subplot(2,2,4)
% rejectplot=(1./rejectframe).*FitParams(:,4);
% plot([1:nframes], FitParams(:,4),'-', [1:nframes], rejectplot,'or'); title('phase')




