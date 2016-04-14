function colDrop(db,table,to_drop)
%colDrop(db,'TableName','ColumnName')
%
% Drops the indicated column from the indicated table in the current
% database (db).

dbpe
import CMBHOME.Database.*

a = sprintf('ALTER TABLE %s DROP COLUMN %s',table,to_drop);
db.prepareStatement(a);
[~] = db.query();
