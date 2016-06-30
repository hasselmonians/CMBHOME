function [dacqInfo] = dacq_set(fname)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%            COPY IN FILE                  %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Pulls relevant information out of the .set file
  %fname = fullfile(dataPath,[prefix,'.set']);
  fid = fopen(fname);
  if fid<0; error('%s not found',fname); end
  lbls = cell(1502,1);
  vals = cell(1502,1);
  for i = 1:1502
      temp = fgetl(fid);
      [lbls{i} vals{i}] = strtok(temp);
  end
  fclose(fid);
  
  % trash junk entries
  goodLbls = cellfun(@ischar,lbls,'unif',1);
  lbls = lbls(goodLbls);
  vals = vals(goodLbls);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%            BUILD REFERENCE DATA          %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % define system codes
  modes = {'sig','ref','-sig','-ref','sig-ref','ref-sig','gnd'};

  % filter details
  filt = readInFiltDefaults(lbls,vals);

  % collect ref channel info (reference only)
  refs = nan(8,1);
  for i = 1:8
   refs(i) = collectVal(lbls,vals,['ref_',num2str(i-1)],0);
  end

  % collect EEG info
  eegCh = nan(64,1);
  saveEEG = nan(64,1);
  for i = 1:64
    eegCh(i) = collectVal(lbls,vals,['EEG_ch_',num2str(i)],0);
    saveEEG(i) = collectVal(lbls,vals,['saveEEG_ch_',num2str(i)],0);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%            BUILD DACQINFO STRUCT         %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  dacqInfo.date = collectVal(lbls,vals,'trial_date',1);
  dacqInfo.time = collectVal(lbls,vals,'trial_time',1);

  % collect channel info
  for i = 1:64
    ch(i).channel = i;
    ch(i).gain = collectVal(lbls,vals,['gain_ch_',num2str(i-1)],0);
    if strcmp(collectVal(lbls,vals,'sw_version',1),' 1.2.2.14')
        ch(i).ref = collectVal(lbls,vals,['b_in_ch_',num2str(i-1)],0); % newer version indicates channel directly
    else
        ch(i).ref = refs(collectVal(lbls,vals,['b_in_ch_',num2str(i-1)],0)+1); % older version refered to ref list
    end
    ch(i).mode = modes{collectVal(lbls,vals,['mode_ch_',num2str(i-1)],0)+1};
    ch(i).filt = filtDetails(i,lbls,vals,filt);
  end
  dacqInfo.chans = ch;
  
	% condense EEG info
  chInds = eegCh(logical(saveEEG)); 
	for chInd = 1:length(chInds)
		if chInds(chInd)==0, continue; end
		eeg(chInd) = ch(chInds(chInd));
	end
	[eeg.saved_4800Hz] = deal(collectVal(lbls,vals,'saveEGF',0));
	dacqInfo.eeg = eeg;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


function val = collectVal(lbls,vals,fld,asStr)
val = NaN;  
tf = ismember(lbls,fld);
  if ~sum(tf)
    warning('dacq_set:bad_field','Field not found: %s',fld);return;
  elseif sum(tf)>1
    warning('dacq_set:bad_field','More than one match found: %s',fld);return;
  end
  val = vals{tf};
  if ~asStr, val = str2num(val); end
end


function filt = readInFiltDefaults(lbls,vals)
  % define types
  filt.types = {'direct','direct+notch','highpass','lowpass','lowpass+notch','custom'};
  filt.respTypes = {'lowpass','highpass','bandpass','bandstop'};
  filt.kindTypes = {'Butterworth','Chebyshev','Bessel'};
  filt.notchFreq = [collectVal(lbls,vals,'notch_frequency',1),' Hz'];
  
  % lowpass default  
  filt.lowpassDef.resp = filt.respTypes{collectVal(lbls,vals,'default_filtresp_lp',0)+1};
  filt.lowpassDef.kind = filt.kindTypes{collectVal(lbls,vals,'default_filtkind_lp',0)+1};
  filt.lowpassDef.freq1 = collectVal(lbls,vals,'default_filtfreq1_lp',0);
  if collectVal(lbls,vals,'default_filtresp_lp',0)>1 % only bandpass and bandstop use both freq1 and freq2
      filt.lowpassDef.freq2 = collectVal(lbls,vals,'default_filtfreq2_lp',0);
  else
      filt.lowpassDef.freq2 = nan;
  end
  filt.lowpassDef.ripple = collectVal(lbls,vals,'default_filtripple_lp',0);
  filt.lowpassDef.dcblock = collectVal(lbls,vals,'default_filtdcblock_lp',0);
  
  % highpass default  
  filt.highpassDef.resp = filt.respTypes{collectVal(lbls,vals,'default_filtresp_hp',0)+1};
  filt.highpassDef.kind = filt.kindTypes{collectVal(lbls,vals,'default_filtkind_hp',0)+1};
  filt.highpassDef.freq1 = collectVal(lbls,vals,'default_filtfreq1_hp',0);
  if collectVal(lbls,vals,'default_filtresp_hp',0)>1 % only bandpass and bandstop use both freq1 and freq2
      filt.highpassDef.freq2 = collectVal(lbls,vals,'default_filtfreq2_hp',0);
  else
      filt.highpassDef.freq2 = nan;
  end
  filt.highpassDef.ripple = collectVal(lbls,vals,'default_filtripple_hp',0);
  filt.highpassDef.dcblock = collectVal(lbls,vals,'default_filtdcblock_hp',0);
end

function filt = filtDetails(i,lbls,vals,filtDef)
    iStr = num2str(i-1);
    type = collectVal(lbls,vals,['filter_ch_',iStr],0);
    filt.type = filtDef.types{type+1};
    switch(type)
        case(0) % 'direct'
            filt.resp = 'n/a';
            filt.kind = 'n/a';
            filt.freq1 = nan;
            filt.freq2 = nan;
            filt.ripple = nan;
            filt.dcblock = nan;
            filt.notch =  'n/a';
        case(1) % direct + notch
            filt.resp = 'n/a';
            filt.kind = 'n/a';
            filt.freq1 = nan;
            filt.freq2 = nan;
            filt.ripple = nan;
            filt.dcblock = nan;
            filt.notch = filtDef.notchFreq;
        case(2) % highpass
            filt.resp = filtDef.highpassDef.resp;
            filt.kind = filtDef.highpassDef.kind;
            filt.freq1 = filtDef.highpassDef.freq1;
            filt.freq2 = filtDef.highpassDef.freq2;
            filt.ripple = filtDef.highpassDef.ripple;
            filt.dcblock = filtDef.highpassDef.dcblock;
            filt.notch = 'n/a';
        case(3) % lowpass
            filt.resp = filtDef.lowpassDef.resp;
            filt.kind = filtDef.lowpassDef.kind;
            filt.freq1 = filtDef.lowpassDef.freq1;
            filt.freq2 = filtDef.lowpassDef.freq2;
            filt.ripple = filtDef.lowpassDef.ripple;
            filt.dcblock = filtDef.lowpassDef.dcblock;
            filt.notch = 'n/a';
        case(4) % lowpass + notch
            filt.resp = filtDef.lowpassDef.resp;
            filt.kind = filtDef.lowpassDef.kind;
            filt.freq1 = filtDef.lowpassDef.freq1;
            filt.freq2 = filtDef.lowpassDef.freq2;
            filt.ripple = filtDef.lowpassDef.ripple;
            filt.dcblock = filtDef.lowpassDef.dcblock;
            filt.notch = filtDef.notchFreq;
        case(5) % custom filter, must read in from file
            resp = collectVal(lbls,vals,['filtresp_ch_',iStr],0);
            filt.resp = filtDef.respTypes{resp+1};
            filt.kind = filtDef.kindTypes{collectVal(lbls,vals,['filtkind_ch_',iStr],0)+1};
            filt.freq1 = collectVal(lbls,vals,['filtfreq1_ch_',num2str(i-1)],0);
            if resp > 1 % only bandpass and bandstop use freq2
                filt.freq2 = collectVal(lbls,vals,['filtfreq2_ch_',num2str(i-1)],0);
            else
                filt.freq2 = nan;
            end
            filt.ripple = collectVal(lbls,vals,['filtripple_ch_',num2str(i-1)],0);
            filt.dcblock = collectVal(lbls,vals,['filtdcblock_ch_',num2str(i-1)],0);
            filt.notch = 'n/a';
    end
end
