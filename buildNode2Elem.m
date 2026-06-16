function node2elem = buildNode2Elem(elements, nNodes)

node2elem = cell(nNodes,1);

for k = 1:size(elements,1)

    for j = 1:3
        node = elements(k,j);
        node2elem{node}(end+1) = k;
    end

end

end