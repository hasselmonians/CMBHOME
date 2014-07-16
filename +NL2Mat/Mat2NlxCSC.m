% Help for file Mat2NlxCSC
% Version 4.1.1
%
% Purpsose
% The purpose of this program is to take extracted Neuralynx data from the matlab environment and create the appropriate Neuralynx file.
% 
% Composition
% This program is a matlab mex file compilation of matlab and C++ source code.
% C++ source code is used to extract data from the matlab environment in
% the form of matlab arrays or matrices and write that data to a Neuralynx file.
% 
% Note: All Timestamps must be in microseconds.
% 
% NumFields = 6;
% Field List
%     1. Timestamps   
%     2. Channel Numbers
%     3. Sample Frequency
%     4. Number of Valid Samples
%     5. Samples
%     6. Header
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
%     Parameter1 - Output File Name and Path
%                  Requirements
%                  1. The file must be a Neuralynx CSC(Continuously Sampled Channel) data file.
%                       If writing to a new file, the correct extension is ".ncs"
%                  2. The filename must be a string and a row vector.
%                  Example:  
%                      Filename = 'C:\TestDirectory\TestFile.ncs'
%
%     Parameter2 - Append File Variable - This denotes whether to append to an already existing file or to overwrite an already existing file.
%                  Requirements
%                  1. The value must be a 1 or a 0.  A 1 will append to the file, a 0 will not append to the file.
%                  Example:   We want to append to the file.
%                      AppendFile = 1;
%                  
%     Parameter3 - Extraction Mode - This value indicates what records to extract from Matlab and write to file.  See Extract Modes Above for Description.
%                  Requirements
%                  1. The value must be in the range 1..NumModes.
%                  Example:   We want to extract all the records from matlab and write them to file.
%                      ExtractMode = 1;
% 
%     Parameter4 - Extraction Mode Array - This array indicates what records to extract from Matlab and write to file.  Each Selection requires a different array.  They are as follows:
%                  1. Extract All - This mode does not require an Array, please insert a 1 into this variable.
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
%     Parameter5 - Number Of Records In Matlab Arrays - This denotes the size of each array we are dealing with.
%                  Requirements
%                  1. The Array must be an integer.
%                  2. All Matlab Matrices or Arrays must have the same number of records which is equal to this parameter.
%                  Example:   If our Matlab TimeStamp Matrix contains 100 timestamps then we set our variable to 100.
%                        NumRecs = 100;
%
%     Parameter6 - Field Selection Array - This denotes which fields to extract from matlab and write to file, any combination is possible.
%                  Requirements
%                  1. The Array must contain only 1's and 0's.  A 1 will extract the field, a 0 will not extract the field.
%                  2. The Array must be of size NumFields.
%                  Example:   We want to extract fields for Timestamps and Samples Only.
%                      FieldSelection(1) = 1;
%                      FieldSelection(2) = 0;
%                      FieldSelection(3) = 0;
%                      FieldSelection(4) = 0;
%                      FieldSelection(5) = 1;
%                      FieldSelection(6) = 0;
%                  
%     Parameter7-12 - Selected Fields that will be written to File. These variables are dependent on the field selections from Parameter 6.
%                  Requirements
%                  1. The variable must match the field Selection Array.  
%                  Example:   If we look at Parameter 6 example were we were going to extract and write to file the TimeStamp
%                  and Samples fields only, then we would have two parameters.  The order must also follow along with the
%                  appropriate field selection. Reversing the order of the parameters would cause incorrect data to be written to file.
%                      Parameter 7 = Timestamps;
%                      Parameter 8 = Samples;
%                  
%
%
% Complete Example:
%         If we combine all examples from above we come up with a typical call to our functions.
% 
%         Function( Filename, AppendFile, ExtractMode, 1, NumRecs, FieldSelection, Timestamps, Samples );
%         
%         Another Example: We want to write all fields(this includes a header) and use record list extraction mode.
% 
%         FieldSelection(1) = 1;
%         FieldSelection(2) = 1;
%         FieldSelection(3) = 1;
%         FieldSelection(4) = 1;
%         FieldSelection(5) = 1;
%         FieldSelection(6) = 1;
% 
%         ExtractMode = 3;
% 
%         ModeArray(1) = 0;
%         ModeArray(2) = 4;
%         ModeArray(3) = 9;
% 
%         Function( Filename, AppendFile, ExtractMode, ModeArray, NumRecs, FieldSelection, Timestamps, ChannelNumbers, SampleFrequency, NumberValidSamples, Samples, Header );
%         

