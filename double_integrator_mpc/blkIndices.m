function matrixIndices = blkIndices(blockIndex, blockSize)
numBlocks = length(blockIndex);
matrixIndices = zeros(numBlocks*blockSize,1);
for i=1:numBlocks
  matrixIndices(blockSize*i-(blockSize-1) : blockSize*i) = ...
      [blockSize*blockIndex(i)-(blockSize-1) : blockSize*blockIndex(i)];
end
end