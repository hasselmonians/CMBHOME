%SEEscript Processes Tetrode Data and creates a CMBhome file.
%   Detailed explanation goes here
%
%Import files from directory.
%
%NOTE: If a CMBobject already exists in the directory the original will be
%imported into the work space. The second CMBobject will be created in the
%directory, and can be loaded into the work space manually.
files = dir('TT*.txt');
for i=1:length(files);
eval(['load ' files(i).name ' -ascii']);
end


%Create variables Sc* from TT* and saves created variables in the
%current directory.
if exist ('TT1');
Sc1 = TT1;
save ('Sc1', 'Sc1');
end
if exist ('TT2');
Sc2 = TT2;
save ('Sc2', 'Sc2');
end
if exist ('TT3');
Sc3 = TT3;
save ('Sc3', 'Sc3');
end
if exist ('TT4');
Sc4 = TT4;
save ('Sc4', 'Sc4');
end
if exist ('TT5');
Sc5 = TT5;
save ('Sc5', 'Sc5');
end
if exist ('TT6');
Sc6 = TT6;
save ('Sc6', 'Sc6');
end
if exist ('TT7');
Sc7 = TT7;
save ('Sc7', 'Sc7');
end
if exist ('TT8');
Sc8 = TT8;
save ('Sc8', 'Sc8');
end
if exist ('TT9');
Sc9 = TT9;
save ('Sc9', 'Sc9');
end
if exist ('TT10');
Sc10 = TT10;
save ('Sc10', 'Sc10');
end
if exist ('TT11');
Sc11 = TT11;
save ('Sc11', 'Sc11');
end
if exist ('TT12');
Sc12 = TT12;
save ('Sc12', 'Sc12');
end
if exist ('TT13');
Sc13 = TT13;
save ('Sc13', 'Sc13');
end
if exist ('TT14');
Sc14 = TT14;
save ('Sc14', 'Sc14');
end
if exist ('TT15');
Sc15 = TT15;
save ('Sc15', 'Sc15');
end
if exist ('TT16');
Sc16 = TT16;
save ('Sc16', 'Sc16');
end


%Create CMBobject in current directory
CMBHOME.Import.NL('n_tetrodes',16);

%Import CMBobject from current directory
%files = dir('CMBobject*.mat');
%root = load(files.name);
