function [markers] = boundary(nodelist, boundary_values)
    markers = zeros(length(nodelist),1);
    for j=1:length(nodelist)
        if nodelist(j)>0 && mod(nodelist(j),2)~=0
            markers(j) = find(boundary_values == nodelist(j));
        else
            markers(j) = 0;
        end
    end
end