function [Spike, spikeWaves, allSpikeWaves, startSpike] = dacq_tetrode(basePath, cutFiles, dotnums, startSpike)
%[root.spike root.wmfs] = dacq_tetrode(dataPath,prefix,varargin)
%
%For a given directory and file prefix, checks for any *.cut files. For all
%if any are found, cross references with .# file to create a CMBHOME.Spike
%object with spike times and spike indices.
%
%Defaults to overwrite existing spike files (specify 0 for third argument
%if undesired).

%Bill, 2012.05.31

Spike = CMBHOME.Spike; %Initialize Spike
packetLength = 216; %Given in DacqUSBFileFormats.PDF

if isempty(strfind(cutFiles{1},'clu'))
    type='Axona';
else
    type='KlustaKwik';
end

%% 
for i = 1:length(cutFiles)
    
    if strcmp(type,'Axona')
        cut = cutFiles{i}(1:end-4);
        inds = strfind(cut,'_');
        dotnum = dotnums{i};
        tetnum = str2num(cut(inds(end)+1:end));
    else
        cut = cutFiles{i};
        inds = strfind(cut,'.');
        tetnum = str2num(cut(inds(end)+1:end));
        dotnum = dotnums{i};
    end

    %open raw event file and error check
    fid = fopen(dotnum,'r','ieee-be'); %file formats PDF specifies big-endian

    if (fid == -1)
      warning('Cannot open file %s, will place hold',file);
    else
      for k = 1:14, fgetl(fid); end %14 full lines of metadata

      headerSize = ftell(fid) + 10; %skip "data_start" then immediately begin binary
      fseek(fid,0,1); % find end of file
      Npackets = ((ftell(fid)-headerSize)/packetLength); % determine total number of packets in file

      fseek(fid,headerSize,-1); %rewind the file to end of header

      % exit smoothly if failed to read in an integer number of packets
      if ~isinteger(Npackets), Npackets = floor(Npackets); end;
    %             warning('Npackets is invalid, skipping waveform read');
    %             Spike=[];
    %             spikeWaves.mean = [];
    %             spikeWaves.std = [];
    %             allSpikeWaves = [];
    %             outFile = {};
    %             return
    %           end

      %For each packet, read the timestamp and 50 8bit samples for all 4 channels                 
      time_stamps = nan(1,Npackets);
      wave1 = int8(nan(50,Npackets));
      wave2 = int8(nan(50,Npackets));
      wave3 = int8(nan(50,Npackets));
      wave4 = int8(nan(50,Npackets));

      for j = 1:Npackets
        time_stamps(j)= fread(fid,1,'int32')/96000; %divide by time base
        wave1(:,j)=fread(fid,50,'int8');

        [~] = fread(fid,1,'int32')/96000;
        wave2(:,j)=fread(fid,50,'int8');

        [~] = fread(fid,1,'int32')/96000;
        wave3(:,j)=fread(fid,50,'int8');

        [~] = fread(fid,1,'int32')/96000;
        wave4(:,j)=fread(fid,50,'int8');
      end
      %disp(num2str(Npackets))
      %% Now cross reference the .cut files to get the spike times
      cfid=fopen([basePath, filesep, cutFiles{i}],'r');
      if cfid == -1
        warning('Cannot open file %s',cutfile);
      else

        % check cut file type by filename:
        %   axona native version has '.cut'
        %   klustakwik has '.clu.'
        if strcmp(type,'Axona')

            % zoom past header
            fseek(cfid,0,-1);
            nu = fgetl(cfid);
            inds = strfind(nu,' ');
            nu=nu(inds(end)+1:end);
            nu=str2num(nu);
            fseek(cfid,0,-1);
            header_lines = 5+nu*3+1;
            for k = 1:header_lines; fgetl(cfid); end

            % pull them out
            isSpike = [];
            while ~feof(cfid)
              a = str2num(fgetl(cfid));
              isSpike = [isSpike a];
            end

            isSpike = isSpike';
        elseif strcmp(type,'KlustaKwik')
            isSpike = [];
            while ~feof(cfid)
              a = str2num(fgetl(cfid));
              isSpike = [isSpike a];
            end
            isSpike = isSpike(:);
        end

        % filter to this session only
        isSpike = isSpike(startSpike(i):startSpike(i)+Npackets-1);
        startSpike(i) = startSpike(i)+Npackets;
        
        %% Run through isSpike and make the .CMBHOME.Spike structure

        % KlustaKwik uses cluster 1 for unsorted events, subtract 1
        % from all avlues
        if strcmp(type,'KlustaKwik'), isSpike = isSpike - 1; end

        clusts=unique(isSpike)';
        clusts(clusts==0) = [];

        if isempty(clusts), 
          Spike(tetnum,1) = CMBHOME.Spike();
          spikeWaves(tetnum,1).mean = nan;
          spikeWaves(tetnum,1).std = nan;
          k=0;
        else
            for k = clusts
                spkInds = (isSpike == k);
                Spike(tetnum,k) = CMBHOME.Spike('i',find(spkInds)','ts',time_stamps(spkInds)','tet_name',cutFiles{i});
                spikeWaves(tetnum,k).mean = [...
                  mean(wave1(:,spkInds),2),...
                  mean(wave2(:,spkInds),2),...
                  mean(wave3(:,spkInds),2),...
                  mean(wave4(:,spkInds),2)];
                spikeWaves(tetnum,k).std = [...
                  std(double(wave1(:,spkInds)),[],2),...
                  std(double(wave2(:,spkInds)),[],2),...
                  std(double(wave3(:,spkInds)),[],2),...
                  std(double(wave4(:,spkInds)),[],2)];
                allSpikeWaves(tetnum,k).raw = cat(3,...
                  double(wave1(:,spkInds)),...
                  double(wave2(:,spkInds)),...
                  double(wave3(:,spkInds)),...
                  double(wave4(:,spkInds)));
                allSpikeWaves(tetnum,k).cell = [tetnum k];
                allSpikeWaves(tetnum,k).ts = time_stamps(spkInds)';
                clear spkInds
            end
        end
        % collect uncut waves as last unit
        %%{
        k=k+1;
        allSpikeWaves(tetnum,k).cell = [tetnum 0];
        spkInds = (isSpike == 0);
        allSpikeWaves(tetnum,k).raw = cat(3,...
          double(wave1(:,spkInds)),...
          double(wave2(:,spkInds)),...
          double(wave3(:,spkInds)),...
          double(wave4(:,spkInds)));
        allSpikeWaves(tetnum,k).ts = time_stamps(spkInds)';
        %}
        fclose(cfid);
      end
      fclose(fid);
    end
end


