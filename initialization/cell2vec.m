function [ Vec ] = cell2vec( Cell )
% converting cell array to vectror array
for ii=1:numel(Cell)
   Vec(ii) = Cell{ii}; 
end

end

