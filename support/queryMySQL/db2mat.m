function [devData err] = db2mat(db,indices,imps)
%function [devData err] = db2mat(db,indices,imps)
%Author: Jonathon Parrish
%Company: Guralp Systems Limited
%Created: May 25th, 2011
%Last Edit: May 27th, 2011
%Purpose: Retrieves a mat file from the database.
%
%Inputs: db - a MYSQLDatabase object
%        indices - double array [te_index,rawdata_id]
%        imps - cell array of packages to load
%
%Outputs: err - empty normally, but contains an error message otherwise
%
%Assumptions: fileName is something sensible supported by the operating
%             system
%
%Dependencies: MYSQLDatabase.m, err2con.m

    err = [];
    devData = [];
    %% Load Packages
    if nargin~=3
        err = ('mat2db: Not enough input arguments');
        return;
    end
    try %#ok<*TRYNC>
        nimps = length(imps);
        for n=1:nimps
            import(char(imps(n)));
        end
    end
    %% Input Validation
    try
        x = db.dbConn;          %#ok<*NASGU> %The database connection
        x = db.hostName;        %The database host
        x = db.schema;          %The schema to use (test7, ims, etc....)
        x = db.userName;        %The userName for the database
        x = db.password;        %The database password
    catch err
        err = err2con(err.stack,err.message);
        return;
    end
    %% Set table data
    tableName = 'rawdata_instance';
    testvar{1} = indices(1);
    testvar{2} = indices(2);
    %% Get from the table
    sqlstatement =['SELECT data ',...
                   'FROM ' tableName ' ',...
                   'WHERE (te_index,rawdata_id) =("{Si}","{Si}")'];
    db.prepareStatement(sqlstatement,testvar);
    result = db.query();
    if isempty(result)
        err = err2con(dbstack,db.errMsg);
    end
    fileIn = int16(cell2mat(result.data)');
    fileOut = repmat(uint8(1),1,length(fileIn));
    for v=1:length(fileIn)
        if fileIn(v)<0
            fileOut(v) = uint8(fileIn(v)+256);
        else
            fileOut(v) = uint8(fileIn(v));
        end
    end
    %% Write to mat file
    fileName = 'tempData.mat';
    fid = fopen(fileName,'w');
    fwrite(fid,fileOut);
    fclose(fid);
    %% Get the variable devData out of the file and return it
    devData = importdata(fileName);
end