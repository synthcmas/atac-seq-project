mat = mat[,names(labels)]
cells_data=colnames(mat)
cellLines_csv <- data.frame(cellLines = labels)
rownames(cellLines_csv) = cells_data
rm(cells_data)
mat = as(mat,'dgCMatrix')
print(class(mat))

cisTopicObject <- createcisTopicObject(mat, project.name=paste(tolower(org),'_data',sep=''))
rm(mat)
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cellLines_csv)
print(class(cisTopicObject))
rm(cellLines_csv)
print('cisTopicObject created!')
