function out = fns(statement)
% Takes in a command to be passed through queryMySQL and turns 'NULL' to
% NULL in order to make it a command, rather than a string.

import CMBHOME.Database.*

out = strrep(statement,['''','NULL',''''],'NULL');

end