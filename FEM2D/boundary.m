function [markers] = boundary(nodelist, boundary_values)
    markers = zeros(length(nodelist),1);
    for j=1:length(nodelist)
        if nodelist(j)>0 && mod(nodelist(j),2)~=0
            if ismember(nodelist(j), boundary_values)
                markers(j) = find(boundary_values == nodelist(j));
            else
                if j==1
                    markers(j) = 1;
                elseif j==2
                    markers(j) = 4;
                elseif j==3
                    markers(j) = 2;
                else
                    markers(j) = 3;
                end
            end
        else
            markers(j) = 0;
        end
    end
end