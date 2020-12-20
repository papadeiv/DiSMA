if mrk_j==mrk_k % diagonal entries
    D(mrk_j,mrk_k) = 2*mu(x)/mesh_size;
    C(mrk_j,mrk_k) = 0;
    R(mrk_j,mrk_k) = (2/3)*sigma(x)*mesh_size;
elseif mrk_k==mrk_j-1 % off-diagonal entries: lower diagonal
    D(mrk_j,mrk_k) = -mu(x)/mesh_size; 
    C(mrk_j,mrk_k) = -0.5*beta(x);
    R(mrk_j,mrk_k) = (mesh_size*sigma(x))/6;
elseif mrk_k==mrk_j+1 % off-diagonal entries: upper diagonal
    D(mrk_j,mrk_k) = -mu(x)/mesh_size; 
    C(mrk_j,mrk_k) = 0.5*beta(x);
    R(mrk_j,mrk_k) = (mesh_size*sigma(x))/6;
end