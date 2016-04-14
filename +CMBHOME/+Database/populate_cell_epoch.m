function populate_cell_epoch(row)
% Part of the batch processing stream to populate the database using
% condor.

dbpe
import CMBHOME.Database.*

%% Uses the output from condor prepared
db = connect('root','secret_password');
db = connect('root','secret_password');

target_file = strcat(dbp,
fid = fopen('to_workon.txt'); % Open the temporary file

row = row+1; % Correct for differences in initialization 

% Now loop through until you get to the correct line
if row == 1
    target = fgetl(fid);
else
    for i = 1:row-1
        fgetl(fid);
        target = fgetl(fid);
    end
end
fclose(fid);

%% Break up the inputs
dbp = '/Volumes/External/Dropbox/';

[fname target] = strtok(target,',');
target = target(2:end);
[label,target] = strtok(target,',');
target = target(2:end);
[epoch,target] = strtok(target,',');
sid = target(2:end);
sid = num2str(sid);

%% Process
fname = [dbp '/UnitRecordingData' fname];
c = load(fname);
root = c.root;
clear c;
root = root.FixDir;
root = root.FixPos;

root.epoch = eval(epoch);
[lv,~] = size((root.cells));

%{
db.prepareStatement('SELECT tetrode_id from cell where (id IN (select cell_id from cell_session where (session_id = "{S}")))',num2str(sid));
[a] = db.query();
tids = a.tetrode_id;
%}

db.prepareStatement('select cell_id from cell_session where (session_id = "{S}")',num2str(sid));
[a] = db.query();
cids = a.cell_id;

label='full session';
for i = 1:lv
    root.cel = root.cells(i,:);
    cid = num2str(cids(i));
    gridness1 = root.Gridness(root.cel);
    gridness2 = root.Gridness2(root.cel);
    gridness3 = root.Gridness(root.cel,'grid3',1);

    si = root.SpatialInformation(root.cel);
    ratemap = root.RateMap(root.cel);
    acorr = root.AutoCorr(root.cel,'supress_plot',1);

    nspk = length(root.cel_x{1});
    f = nspk / diff(root.epoch);
    [fintrins thetaness] = root.IntrinsicFrequency(root.cel,1,[-1 -1],'chronux_spectrum',0);
    theta_skip = root.ThetaSkipping(root.cel);
    if isempty(theta_skip); theta_skip = 'NULL';end

     phasepre_rho = 'NULL';
     phasepre_p = 'NULL';

    lfp_used = root.active_lfp;
    theta = root.cel_thetaphase{1};
    theta = theta(~isnan(theta));
    mean_theta_phase = CMBHOME.Utils.circ_mean(theta);
    theta_phase_mr = CMBHOME.Utils.circ_r(theta);
    
    if isempty(lfp_used)
    	lfp_used = 'NULL';
    	mean_theta_phase = 'NULL';
    	theta_phase_mr = 'NULL';
    end	

    hd = root.headdir;
    hd = deg2rad(hd);
    hd = hd(~isnan(hd));
    hd_dir = CMBHOME.Utils.circ_mean(hd);
    hd_mrl = CMBHOME.Utils.circ_r(hd); 
    u2 = root.HDWatsonsU2(root.cel);

    [f2 v] = root.VelocityRate(root.cel);
    [~,speed_mod_r,~,speed_mod_s,speed_mod_int] = CMBHOME.Utils.LinearRegression(v,f2);


    %% Put in the database
    sv = 'INSERT INTO cell_epoch (cell_id,session_id,epoch_label,gridness1,gridness2,gridness3,si,ratemap,acorr,nspk,f,fintrins,thetaness,theta_skip, phasepre_rho, phasepre_p,lfp_used,mean_theta_phase,theta_phase_mr,hd_dir,hd_mrl,watsonsu2,speed_mod_r,speed_mod_s,speed_mod_int) VALUES(';
    sv2 = sprintf(' ''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'',''%s'')', ...
                    num2str(cellSession.cell_id(i)), ...
                    num2str(cellSession.session_id(i)), ...
                    label, ...
                    num2str(gridness1), ...
                    num2str(gridness2), ...
                    num2str(gridness3), ...
                    num2str(si), ...
                    mat2str(ratemap), ...
                    mat2str(acorr), ...
                    num2str(nspk), ...
                    num2str(f), ...
                    num2str(fintrins), ...
                    num2str(thetaness), ...
                    num2str(theta_skip), ...
                    num2str(phasepre_rho), ...
                    num2str(phasepre_p),...
                    num2str(lfp_used),...
                    num2str(mean_theta_phase),...
                    num2str(theta_phase_mr), ...
                    num2str(hd_dir),...
                    num2str(hd_mrl),...
                    num2str(u2),...
                    num2str(speed_mod_r),...
                    num2str(speed_mod_s),...
                    num2str(speed_mod_int));

    sv = [sv sv2];
    sv = CMBHOME.Database.fns(sv);
    db.prepareStatement(sv);
    [a] = db.query();
    
end