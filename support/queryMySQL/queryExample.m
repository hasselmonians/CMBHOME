%% set path
addpath(fullfile(pwd, 'src'));
javaaddpath('lib/mysql-connector-java-5.1.6/mysql-connector-java-5.1.6-bin.jar');

%% import classes
import edu.stanford.covert.db.MySQLDatabase;

%% create database connection
db = MySQLDatabase('covertlab.stanford.edu', 'test', 'test', 'test');

%% create file
data = char(65:85);
fname = tempname;
fid = fopen(fname, 'wb');
fwrite(fid, data);
fclose(fid);

data2 = char(66:86);
fname2 = tempname;
fid = fopen(fname2, 'wb');
fwrite(fid, data2);
fclose(fid);

data3 = char(85:-1:65);
fname3 = tempname;
fid = fopen(fname3, 'wb');
fwrite(fid, data3);
fclose(fid);

data4 = uint8(85:-1:65);
fname4 = tempname;
fid = fopen(fname4, 'wb');
fwrite(fid, data4);
fclose(fid);

%% send, retrieve blob using sql
db.prepareStatement('INSERT test (DATA, ID) VALUES("abc", "{Si}")', 10001);
db.query();

db.prepareStatement('UPDATE test SET DATA = "{F}" WHERE id = "{Si}"', fname, 10001);
db.query();

db.prepareStatement('SELECT * FROM test WHERE id = "{Si}"', 10001);
result = db.query();

assert(isequal(data, result.data{1}));

%% send, retrieve blob using stored procedures
db.prepareStatement('CALL testBlobIn("{Si}","{F}")', 10001, fname2);
db.query();

db.prepareStatement('CALL testBlobOut("{Si}")', 10001);
result = db.query();

assert(isequal(data2, result.data{1}));

%% Jonathon parish example #1
index = 1;

db.prepareStatement('INSERT INTO test2 (`index`, data) VALUES ("{Si}", "{F}")', index, fname3);
assert(~isequal(-1, db.query));

db.prepareStatement('SELECT * FROM test2 WHERE `index` = "{Si}"', index);
result = db.query();

assert(isequal(data3, result.data{1}));

db.prepareStatement('DELETE FROM test2 WHERE `index` = "{Si}"', index);
assert(~isequal(-1, db.query));

%% Jonathon parish example #2
index = 1;

db.prepareStatement('INSERT INTO test3 (`index`, data) VALUES ("{Si}", "{F}")', index, fname4);
assert(~isequal(-1, db.query));

db.prepareStatement('SELECT * FROM test3 WHERE `index` = "{Si}"', index);
result = db.query();

assert(isequal(data4, result.data{1}));

db.prepareStatement('DELETE FROM test3 WHERE `index` = "{Si}"', index);
assert(~isequal(-1, db.query));

%% Hachi Manzur example
sql = 'select `%% Ans`, `%% Corrects`, `Session No` from `test_4` where `Session No` in (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22) AND `%% Ans`>50';
db.prepareStatement(sql);
[result, colNames] = db.query();
assert(~isequal(-1, result));
assert(~isequal(-1, colNames));

%% Joseph Kuehl example -- date, time, datetime datatypes
%Table with 4 columns:
%
%Name               Datatype   Non-null
%================   ========   ========
%datenull           date       no
%datenonnull        date       yes
%datetimenull       datetime   no
%datetimenonnull    datetime   yes

%clear data
db.prepareStatement('DELETE FROM datetest');
db.query();

%test-1
time1 = '10:54:56';
date1 = '2012-05-30';
db.prepareStatement(sprintf('INSERT INTO `test`.`datetest` (`id`, `timenull`, `timenonnull`, `datenull`, `datenonnull`, `datetimenull`, `datetimenonnull`) VALUES (1, NULL, "%s", NULL, "%s", NULL, "%s %s")', time1, date1, date1, time1));
db.query();

sql = sprintf('SELECT * FROM datetest where `id` = 1');
db.prepareStatement(sql);
result = db.query();

assert(isequal(result.datenull, {''}));
assert(isequal(result.timenull, {''}));
assert(isequal(result.datetimenull, {''}));

assert(isequal(result.datenonnull, {date1}));
assert(isequal(result.timenonnull, {time1}));
assert(isequal(result.datetimenonnull, {[date1 ' ' time1]}));

%test-2
time2 = '00:00:00';
date2 = '0000-00-00';
db.prepareStatement(sprintf('INSERT INTO `test`.`datetest` (`id`, `timenull`, `timenonnull`, `datenull`, `datenonnull`, `datetimenull`, `datetimenonnull`) VALUES (2, NULL, "%s", NULL, "%s", NULL, "%s %s")', time2, date2, date2, time2));
db.query();

sql = sprintf('SELECT * FROM datetest where `id` = 2');
db.prepareStatement(sql);
result = db.query();

assert(isequal(result.datenull, {''}));
assert(isequal(result.timenull, {''}));
assert(isequal(result.datetimenull, {''}));

assert(isequal(result.datenonnull, {date2}));
assert(isequal(result.timenonnull, {time2}));
assert(isequal(result.datetimenonnull, {[date2 ' ' time2]}));

%% Agustin example -- tiny datatype
%clear data
db.prepareStatement('DELETE FROM tinytest');
db.query();

%test
num = 3;
db.prepareStatement('INSERT INTO `test`.`tinytest` (`id`, `datatiny`, `datasmall`, `datamedium`, `dataint`, `databigint`) VALUES (1, "{Si}", "{Si}", "{Si}", "{Si}", "{Si}")', num, num, num, num, num);
db.query();

sql = sprintf('SELECT * FROM tinytest where `id` = 1');
db.prepareStatement(sql);
result = db.query();

assert(isequal(result.datatiny, num));
assert(isequal(result.datasmall, num));
assert(isequal(result.datamedium, num));
assert(isequal(result.dataint, num));
assert(isequal(result.databigint, num));

%% Pierre Martineau example -- bigint datatype
db.prepareStatement('SELECT databigint, dataubigint FROM biginttest');
result = db.query();

assert(isequal(num2str(result.databigint), [
'                   1'
'                 127'
'               32767'
'             8388607'
'          2147483647'
' 9223372036854775807'
' 9223372036854775807'
' 9223372036854775806'
'-9223372036854775806'
'         -2147483646'
'            -8388606'
'              -32766'
'                -126'
    ]));

assert(isequal(num2str(result.dataubigint), [
'                   1'
'                 127'
'               32767'
'             8388607'
'          2147483647'
' 9223372036854775807'
'18446744073709551615'
'18446744073709551614'
'                   0'
'                   0'
'                   0'
'                   0'
'                   0'
    ]));

%% cleanup
delete(fname);
delete(fname2);
delete(fname3);
delete(fname4);
db.close();