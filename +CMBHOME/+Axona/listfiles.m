function fileList=listfiles(dirName)

dirData = dir(dirName);      % # Get the data for the current directory
dirIndex = [dirData.isdir];  % # Find the index for directories
fileList = {dirData(~dirIndex).name}';  % # Get a list of the files

end