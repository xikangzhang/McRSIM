function folderNames = listFolder(dataPath)

folderList = dir(dataPath);
isDir = [folderList(:).isdir];
folderList = folderList(isDir);
folderList(strncmp({folderList.name}, '.', 1)) = [];
folderList(strncmp({folderList.name}, '..', 1)) = [];
folderNames = {folderList.name};

end