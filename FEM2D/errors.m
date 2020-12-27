function errors(h ,A ,Nh , E_inf ,E_L2, E_H1)
    
    % compute error's convergence order in different norms and for different mesh size measures
    alpha_inf_h = polyfit(log10(h),log10(E_inf),1);
    alpha_inf_h_str = sprintf('$ ||u-u_h||_{L^{\\infty}(\\Omega)} \\approx h^{\\;\\:%f} $', alpha_inf_h(1));
    alpha_inf_A = polyfit(log10(A),log10(E_inf),1);
    alpha_inf_A_str = sprintf('$ ||u-u_h||_{L^{\\infty}(\\Omega)} \\approx A^{\\;\\:%f} $', alpha_inf_A(1));
    alpha_inf_Nh = polyfit(log10(Nh),log10(E_inf),1);
    alpha_inf_Nh_str = sprintf('$ ||u-u_h||_{L^{\\infty}(\\Omega)} \\approx N_h^{\\;\\:%f} $', alpha_inf_Nh(1));
    
    alpha_L2_h=polyfit(log10(h),log10(E_L2),1);
    alpha_L2_h_str = sprintf('$ ||u-u_h||_{L^2(\\Omega)} \\approx h^{\\;\\:%f} $', alpha_L2_h(1));
    alpha_L2_A = polyfit(log10(A),log10(E_L2),1);
    alpha_L2_A_str = sprintf('$ ||u-u_h||_{L^2(\\Omega)} \\approx A^{\\;\\:%f} $', alpha_L2_A(1));
    alpha_L2_Nh = polyfit(log10(Nh),log10(E_L2),1);
    alpha_L2_Nh_str = sprintf('$ ||u-u_h||_{L^2(\\Omega)} \\approx N_h^{\\;\\:%f} $', alpha_L2_Nh(1));
    
    alpha_H1_h=polyfit(log10(h),log10(E_H1),1);
    alpha_H1_h_str = sprintf('$ ||u-u_h||_{H^1(\\Omega)} \\approx h^{\\;\\:%f} $', alpha_H1_h(1));
    alpha_H1_A = polyfit(log10(A),log10(E_H1),1);
    alpha_H1_A_str = sprintf('$ ||u-u_h||_{H^1(\\Omega)} \\approx A^{\\;\\:%f} $', alpha_H1_A(1));
    alpha_H1_Nh = polyfit(log10(Nh),log10(E_H1),1);
    alpha_H1_Nh_str = sprintf('$ ||u-u_h||_{H^1(\\Omega)} \\approx N_h^{\\;\\:%f} $', alpha_H1_Nh(1));
    
    % create subplot window for the error plots
    figure('units','normalized','outerposition',[0 0 1 1])
    
    % plot the error decay law with h measured in different norms 
    subplot(1, 3, 1)
    loglog(h, E_inf, 'k--*', 'LineWidth', 2, 'DisplayName', alpha_inf_h_str)
    hold on
    loglog(h, E_L2, 'b--o', 'LineWidth', 2, 'DisplayName', alpha_L2_h_str)
    hold on
    loglog(h, E_H1, 'r--s', 'LineWidth',2, 'DisplayName', alpha_H1_h_str)
    xlabel('Mesh size (h)')
    ylabel('Discretization error')
    leg = legend('Interpreter', 'latex', 'Interpreter', 'latex', 'Interpreter', 'latex');
    set(leg, 'Location', 'best', 'FontSize', 18)
    
    % plot the error decay law with A measured in different norms 
    subplot(1, 3, 2);
    loglog(A, E_inf, 'k--*', 'LineWidth', 2, 'DisplayName', alpha_inf_A_str)
    hold on
    loglog(A, E_L2, 'b--o', 'LineWidth', 2, 'DisplayName', alpha_L2_A_str)
    hold on
    loglog(A, E_H1, 'r--s', 'LineWidth',2, 'DisplayName', alpha_H1_A_str)
    xlabel('Mesh size (A ~ h^2)')
    ylabel('Discretization error')
    leg = legend('Interpreter', 'latex', 'Interpreter', 'latex', 'Interpreter', 'latex');
    set(leg, 'Location', 'best', 'FontSize', 18)
    
    % plot the error decay law with Nh measured in different norms 
    subplot(1, 3, 3);
    loglog(Nh, E_inf, 'k--*', 'LineWidth', 2, 'DisplayName', alpha_inf_Nh_str)
    hold on
    loglog(Nh, E_L2, 'b--o', 'LineWidth', 2, 'DisplayName', alpha_L2_Nh_str)
    hold on
    loglog(Nh, E_H1, 'r--s', 'LineWidth',2, 'DisplayName', alpha_H1_Nh_str)
    xlabel('Mesh size (N_h ~ N_x\cdot N_y)')
    ylabel('Discretization error')
    leg = legend('Interpreter', 'latex', 'Interpreter', 'latex', 'Interpreter', 'latex');
    set(leg, 'Location', 'best', 'FontSize', 18)

end