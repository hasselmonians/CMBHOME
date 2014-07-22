function err = mat2db(fileName,db,indices,imps)
%function err = mat2db(fileName,db,indices,imps)
%Author: Jonathon Parrish
%Company: Guralp Systems Limited
%Created: May 25th, 2011
%Last Edit: May 27th, 2011
%Purpose: Saves a mat file in a database as a BLOB. Returns an error if
%         unsuccessful
%
%Inputs: fileName - 'testfile.mat'
%        db - a MYSQLDatabase object
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
    %% Load Packages
    if nargin~=4
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
    if ~ischar(fileName)
        err = err2con(dbstack,'fileName must be a character array');
        return;
    end
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
    testvar{3} = fileName;
    %% Send to the table
    sqlstatement = ['INSERT INTO ' db.schema '.' tableName,...
                     '(te_index, rawdata_id, data)',...
                     'VALUES ("{Si}","{Si}","{F}")'];
    db.prepareStatement(sqlstatement,testvar);
    result = db.query();
    if result==-1
        err = err2con(dbstack,db.errMsg);
    end
end