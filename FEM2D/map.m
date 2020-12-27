function [vj, vk] = map(indices, B)
    
    global nodes;
    global triangles;

    v1(1) = abs(nodes(triangles(indices(3),1),1));
    vj(1) = nodes(triangles(indices(3),indices(1)),1);
    vk(1) = nodes(triangles(indices(3),indices(2)),1);
    
    vj(2:3) = [nodes(v1(1),3); nodes(v1(1),4)] + B*[nodes(vj(1),3);nodes(vj(1),4)];
    if vk(1)>0
        vk(2:3) = [nodes(v1(1),3); nodes(v1(1),4)] + B*[nodes(vk(1),3);nodes(vk(1),4)];
    else
        vk(2:3) = [nodes(v1(1),3); nodes(v1(1),4)] + B*[nodes(-vk(1),3);nodes(-vk(1),4)];
    end
    
end