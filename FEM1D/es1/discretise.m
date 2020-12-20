function [partition,nodeslist,mesh_size,nDOFs] = discretise(sim_step, length)
    
    nElements = 2^sim_step;
    mesh_size = length/nElements;
    nDOFs = nElements - 1;
    N = nDOFs + 2; % total number of nodes, i.e. inner nodes plus the two endpoints of the domain interval
    nodeslist = zeros(2,N);
    nodeslist(1,:) = linspace(0,length,N);
    nodeslist(2,1) = -1;
    counter = 0;
    for j=2:N-1
        counter = counter + 1;
        nodeslist(2,j) = counter;
    end
    nodeslist(2,N) = -2;
    partition = zeros(nElements,2);
    for j=1:nElements
        partition(j,1) = j;
        partition(j,2) = j+1;
    end
    
end