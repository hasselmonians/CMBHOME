CREATE TABLE `test` (
  `id` bigint(20) default NULL,
  `data` longtext collate latin1_general_ci
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `test2` (
  `index` bigint(20) NOT NULL default '0',
  `data` longtext collate latin1_general_ci,
  PRIMARY KEY  (`index`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `test3` (
  `index` bigint(20) NOT NULL default '0',
  `data` longblob,
  PRIMARY KEY  (`index`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `test_4` (
  `%% Ans` int(11) NOT NULL,
  `%% Corrects` int(11) NOT NULL,
  `Session No` int(11) NOT NULL,
  PRIMARY KEY  (`Session No`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `datetest` (
  `id` int(11) NOT NULL,
  `timenull` time default NULL,
  `timenonnull` time NOT NULL,
  `datenull` date default NULL,
  `datenonnull` date NOT NULL,
  `datetimenull` datetime default NULL,
  `datetimenonnull` datetime NOT NULL,
  PRIMARY KEY  (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `tinytest` (
  `id` int(11) NOT NULL,
  `datatiny` tinyint(4) default NULL,
  `datasmall` smallint(6) default NULL,
  `datamedium` mediumint(9) default NULL,
  `dataint` int(11) default NULL,
  `databigint` bigint(20) default NULL,
  PRIMARY KEY  (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci

CREATE TABLE `biginttest` (
  `id` int(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `datatiny` tinyint(4) default NULL,
  `datasmall` smallint(6) default NULL,
  `datamedium` mediumint(9) default NULL,
  `dataint` int(11) default NULL,
  `databigint` bigint(20) default NULL,
  `dataubigint` bigint(20) unsigned default NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;
-- this inspired by http://rpbouman.blogspot.com/2010/05/mysql-maximum-value-of-integer.html
INSERT INTO biginttest VALUES 
(NULL, 1, 1, 1, 1, 1, 1),
(NULL, ~0 >> 57, ~0 >> 57, ~0 >> 57, ~0 >> 57, ~0 >> 57, ~0>> 57),
(NULL, ~0 >> 57, ~0 >> 49, ~0 >> 49, ~0 >> 49, ~0 >> 49, ~0>> 49),
(NULL, ~0 >> 57, ~0 >> 49, ~0 >> 41, ~0 >> 41, ~0 >> 41, ~0>> 41),
(NULL, ~0 >> 57, ~0 >> 49, ~0 >> 41, ~0 >> 33, ~0 >> 33, ~0>> 33),
(NULL, ~0 >> 57, ~0 >> 49, ~0 >> 41, ~0 >> 33, ~0 >> 1, ~0>> 1),
(NULL, ~0 >> 57, ~0 >> 49, ~0 >> 41, ~0 >> 33, ~0 >> 1, ~0),
(NULL, (~0 >> 57) -1, (~0 >> 49) -1, (~0 >> 41) -1, (~0 >> 33) -1, (~0 >> 1) -1, ~0 -1),
(NULL, -(~0 >> 57) +1, -(~0 >> 49) +1, -(~0 >> 41) +1, -(~0 >> 33) +1, -(~0 >> 1) +1, 0),
(NULL, -(~0 >> 57) +1, -(~0 >> 49) +1, -(~0 >> 41) +1, -(~0 >> 33) +1, -(~0 >> 33) +1, 0),
(NULL, -(~0 >> 57) +1, -(~0 >> 49) +1, -(~0 >> 41) +1, -(~0 >> 41) +1, -(~0 >> 41) +1, 0),
(NULL, -(~0 >> 57) +1, -(~0 >> 49) +1, -(~0 >> 49) +1, -(~0 >> 49) +1, -(~0 >> 49) +1, 0),
(NULL, -(~0 >> 57) +1, -(~0 >> 57) +1, -(~0 >> 57) +1, -(~0 >> 57) +1, -(~0 >> 57) +1, 0);

CREATE DEFINER=`root`@`%` PROCEDURE `testBlobIn`(IN _ID bigint, IN _Data longtext)
BEGIN
  UPDATE test
  SET DATA = _Data
  WHERE id = _ID;
  SELECT _ID id;
END

CREATE DEFINER=`root`@`%` PROCEDURE `testBlobOut`(IN _ID bigint)
BEGIN
  SELECT *
  FROM test
  WHERE id = _ID;
END
