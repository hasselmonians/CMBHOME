function colAdd(db,table,to_add,type,default)
%colAdd(db,'TableName','ColumnName')
%
% Drops the indicated column from the indicated table in the current
% database (db).

import CMBHOME.Database.*

if ~exist('type','var');type = 'FLOAT';end
if ~exist('default','var'),default = 'NULL';end

a = sprintf('ALTER TABLE %s ADD %s %s DEFAULT ''%s'' ',table,to_add,type,default);
a = fns(a);
db.prepareStatement(a);
[~] = db.query();
