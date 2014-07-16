% Help for file Nlx2MatTS
% Version 4.1.1
%
% Purpsose
% The purpose of this program is to extract Neuralynx data from a file and load it into the matlab environment with any combination of fields 
% from the file within a time or record range or using a time or record list.
% 
% Composition
% This program is a matlab mex file compilation of matlab and C++ source code.
% C++ source code is used to read in data from a file and then input the
% data into the matlab enviroment in the form of matlab arrays or matrices.  
% 
% Note: All Timestamps are in microseconds.
% 
% NumFields = 1;
% Field List
%     1. Timestamps   
% 
% NumModes = 5;
% Extraction Modes
%     1. Extract All - This will extract every record from the file into the matlab environment;
%     2. Extract Record Index Range = This will extract every Record whos index is within a range specified by Paramter 5.
%     3. Extract Record Index List = This will extract every Record whos index in the file is the same index that is specified by Paramter 5.
%     4. Extract Timestamp Range = This will extract every Record whos timestamp is within a range specified by Paramter 5.
%     5. Extract Timestamp List = This will extract every Record with the same timestamp that is specified by Paramter 5.
%     
% Input Parameters
% 
%     Parameter1 - Input File Name and Path
%                  Requirements
%                  1. The file must be a Neuralynx Timestamp data file or a
%                  MClust Timestamp data file (.t)
%                  2. The filename must be a string and a row vector.
%                  Example:  
%                      Filename = 'C:\TestDirector\TestFile.nts'
%                  
%     Parameter2 - Field Selection Array - This denotes which fields to extract into matlab, any combination is possible.
%                  Requirements
%                  1. The Array must contain only 1's and 0's.  A 1 will extract the field, a 0 will not extract the field.
%                  2. The Array must be of size NumFields.
%                  Example:   We want to extract fields for Timestamps.
%                      FieldSelection(1) = 1;
%                  
%     Parameter3 - Extract Header Value - This value indicates whether to extract the header from the file into the matlab environment (if one exists).
%                  Requirements
%                  1. The value must be a 1 or a 0.  A 1 will extract the header, a 0 will not extract the header.
%                  Example:   We want to extract the header.
%                      ExtractHeader = 1;
%                  
%     Parameter4 - Extraction Mode - This value indicates what records to extract from the file.  See Extract Modes Above for Description.
%                  Requirements
%                  1. The value must be in the range 1..NumModes.
%                  Example:   We want to extract all the records from the file into matlab.
%                      ExtractMode = 1;
% 
%     Parameter5 - Extraction Mode Array - This array indicates what records to extract from the file.  Each Selection requires a different array.  They are as follows:
%                  1. Extract All - This mode does not require an Array, this parameter should be left blank.
%                  2. Extract Record Index Range.
%                     Requirements
%                         1. The Array must have 2 elements only.
%                         2. The 2 elements must be in increasing order.
%                     Example: We wish to extract the first 10 records from the file.
%                         ModeArray(1) = 0;
%                         ModeArray(2) = 9;
%                  3. Extract Record Index List.
%                     Requirements
%                         1. The Array may have any number of elements( its recommended that the size does not exceed the number of records in the file ).
%                         2. The elements must be in increasing order.
%                     Example: We wish to extract the 1st, 5th and 10th record from the file.
%                         ModeArray(1) = 0;
%                         ModeArray(2) = 4;
%                         ModeArray(3) = 9;
%                  4. Extract Timestamp Range.
%                     Requirements
%                         1. The Array must have 2 elements only.
%                         2. The 2 elements must be in increasing order.
%                     Example: We wish to extract records with time range 10000000 to 20000000.
%                         ModeArray(1) = 10000000;
%                         ModeArray(2) = 20000000;
%                  5. Extract Timestamp List.
%                     Requirements
%                         1. The Array may have any number of elements( its recommended that the size does not exceed the number of records in the file ).
%                         2. The elements must be in increasing order.
%                     Example: We wish to extract records with the timestamps 10000000, 15000000 and 20000000
%                         ModeArray(1) = 10000000;
%                         ModeArray(2) = 15000000;
%                         ModeArray(3) = 20000000;
% 
% 
% Output Parameters
%
%   The 1 variable Nlx2MatTS_v3.dll outputs is as follows:
%   TimeStamps and NlxHeader in that order.
%   TimeStamps is a one dimensional matrices,
%   NlxHeader is a one dimensional cell array, while all others are numeric.
%
%   Note: NlxHeaders are Neuralynx Headers generated by Neuralynx Software, however they may not always be present in a Neuralynx file.
%       If a NlxHeader is not found an empty cell matrix is returned.
%       Headers for MClust(.t) files are not processed due to a lack of formatting.
%
%       
%     The number of output parameters depends on the number of fields selected to extract from the file plus the extraction of the header(if selected).
% 
%     Field and Header Variables - The number of field variables must exactly match the number of fields chosen from parameter 2.
%         
%         Example: If we continue with the example from Parameter2 and we know that we want the Timestamp field extracted, the output will be as follows:
% 
%         [TimestampVariable, HeaderVariable] = Function( Filename, FieldSelection, ExtractHeader, ExtractMode );
%         
%         Example: We want to extract all fields, no header and use record list extraction mode.
% 
%         FieldSelection(1) = 1;
% 
%         ExtractHeader = 0;
%         
%         ExtractMode = 3;
% 
%         ModeArray(1) = 0;
%         ModeArray(2) = 4;
%         ModeArray(3) = 9;
% 
%         [Timestamp] = Nlx2MatTS( Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
%         
%
% Matlab Tips
%
% Some helpful matlab usage tips and commands( commands will be written inside single quotes ' ':
%
%  'whos' - This command will list all current matlab variables, their dimensions and byte sizes. For Example:
%
% >> whos
%    Name                     Size         Bytes  Class
%
%  ChannelNumbers           1x101          808  double array
%  NlxHeader               17x1           2366  cell array
%  NumberValidSamples       1x101          808  double array
%  SampleFrequencies        1x101          808  double array
%  Samples                512x101       413696  double array
%  TimeStamps               1x101          808  double array
%  ans                      1x1              8  double array
%
% Grand total is 52535 elements using 419302 bytes
%
%
%  Indexing Matrices - Using the colon(:) operator.
%       The colon operator can be used to access all data in a specified dimension.  For Example:
%       If you wanted to access all 512 Samples from the above display for record 50.  The following
%       syntax would be required.
%
% >> Samples( :,50 );
%
%       This will list all 512 values for record number 50 vertically.
%
%   Formating numeric data display.
%         When timestamps become large matlab defaults to a scientific notation of numeric values.
%         In order to view the complete values use the following command: 'format long g' For Example:
%         
% >> TimeStamps(2)
% 
% ans =
% 
%   1.0000e+013
% 
% >> format long g
% >> TimeStamps(2)
% 
% ans =
% 
%             10000000000000
% 
% >> 