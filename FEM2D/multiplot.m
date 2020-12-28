function [] = multiplot(Ns, m, z_h, z)

    global triangles;
    global nodes;
    global Ndof;
    x = nodes(:,3);
    y = nodes(:,4);
    % create subplot window for the iterate discrete solution 
    subplot(2,(Ns+1)/2,m+1);
    % plot the discrete surface for the approximate solution
    grid = triangles(:,1:3);
    if Ndof == 6
        scatter3(x,y,z_h,40,z_h,'MarkerFaceColor','flat','LineWidth',1.5,'MarkerEdgeColor','k')
        %colormap(jet);
    else
        trisurf(grid, x, y, z_h, 'LineWidth', 1.5);
    end
    colormap(summer);
    % plot the triangulation underneath the surface
    global geom;
    draw_grid(geom,1)
    % name the plot according to its mesh iteration
    title('Iteration = ', m+1)
    % change azimuth's and elevation's angle of the camera
    if m == Ns
        view(45,0)
    else
        view(75,30)
    end
    % label the axis
    xlabel('x')
    ylabel('y')
    
%     if m == Ns
%         for idx = 0:m-1
%             figure(1)
%             subplot(2,(Ns+1)/2,idx+1);
%             trisurf(grid, x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'LineStyle', '--');
%             colormap(summer);
%         end
%     end
end