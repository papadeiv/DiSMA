function [] = multiplot(multi, Ns, m, z_h, z)

    global triangles;
    global nodes;
    global Ndof;
    global geom;
    
    x = nodes(:,3);
    y = nodes(:,4);
    grid = triangles(:,1:3);
    
    colormap(summer);
    
    if multi=='Y'
        % create subplot window for the iterate discrete solution
        subplot(2,(Ns+1)/2,m+1);
        % plot the discrete surface for the approximate solution
        if Ndof == 6
            scatter3(x,y,z_h,40,z_h,'MarkerFaceColor','flat','LineWidth',1.5,'MarkerEdgeColor','k')
        else
            trisurf(grid, x, y, z_h, 'LineWidth', 1.5);
        end
        % plot the triangulation underneath the surface
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
    else
        if Ndof == 6
            trisurf(grid, x, y, z_h, 'LineWidth', 1.5);
            draw_grid(geom,1)
            %scatter3(x,y,z_h,40,z_h,'MarkerFaceColor','flat','LineWidth',1.5,'MarkerEdgeColor','k')
        else
            trisurf(grid, x, y, z_h, 'LineWidth', 1.5);
            draw_grid(geom,1)
        end
        view(45,0)
        %view(150,6)
        xlabel('x')
        xlabel('x')
        %trisurf(grid, x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'LineStyle', '--');
    end
    
end