function SEIR_DBD_Animation
    % Parameter
    beta = 0.5;
    eta = 0.5;
    delta_eh = 0.2;
    theta_ev = 0.3;
    gamma_h = 0.1;
    mu_h = 0.01;
    mu_v = 0.05;
    Nh = 1000;
    Nv = 5000;
    p = 0.05; % Vaksinasi aktif

    % Waktu simulasi
    tspan = [0 200];
    y0 = [4000; 500; 500; 900; 50; 50; 0];  % [Sv, Ev, Iv, Sh, Eh, Ih, Rh]

    % Simulasi ODE
    [t, y] = ode45(@(t, y) odefun(t, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p), tspan, y0);

    % Ekstrak kompartemen
    Sv = y(:,1); Ev = y(:,2); Iv = y(:,3);
    Sh = y(:,4); Eh = y(:,5); Ih = y(:,6); Rh = y(:,7);

    % Siapkan animasi
    figure;
    for i = 1:length(t)
        % Subplot manusia
        subplot(2,1,1);
        bar([Sh(i), Eh(i), Ih(i), Rh(i)], 'FaceColor', [0.2 0.6 0.9]);
        title(['Populasi Manusia pada Hari ke-', num2str(round(t(i)))]);
        set(gca, 'XTickLabel', {'S_h', 'E_h', 'I_h', 'R_h'});
        ylim([0 Nh]);
        ylabel('Jumlah Individu');
        grid on;

        % Subplot nyamuk
        subplot(2,1,2);
        bar([Sv(i), Ev(i), Iv(i)], 'FaceColor', [0.9 0.4 0.4]);
        title(['Populasi Nyamuk pada Hari ke-', num2str(round(t(i)))]);
        set(gca, 'XTickLabel', {'S_v', 'E_v', 'I_v'});
        ylim([0 Nv]);
        ylabel('Jumlah Individu');
        grid on;

        drawnow; % Update figure
        pause(0.05); % Jeda untuk kecepatan animasi
    end
end

function dydt = odefun(~, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p)
    Sv = y(1); Ev = y(2); Iv = y(3);
    Sh = y(4); Eh = y(5); Ih = y(6); Rh = y(7);

    dSv = mu_v * Nv - (eta * Sv * Ih / Nh) - mu_v * Sv;
    dEv = (eta * Sv * Ih / Nh) - theta_ev * Ev - mu_v * Ev;
    dIv = theta_ev * Ev - mu_v * Iv;

    dSh = mu_h * Nh - (beta * Sh * Iv / Nh) - p * Sh - mu_h * Sh;
    dEh = (beta * Sh * Iv / Nh) - delta_eh * Eh - mu_h * Eh;
    dIh = delta_eh * Eh - gamma_h * Ih - mu_h * Ih;
    dRh = gamma_h * Ih + p * Sh - mu_h * Rh;

    dydt = [dSv; dEv; dIv; dSh; dEh; dIh; dRh];
end
