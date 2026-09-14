function SEIR_DBD_Model_with_Comparison_Plot
    % Parameter tetap
    beta = 0.16667;
    eta = 0.1428;
    delta_eh = 0.375;
    theta_ev = 0.75;
    gamma_h = 0.328833;
    mu_h = 0.0000002;
    mu_v = 0.0323;
    Nh = 9419600;
    Nv = 5000;
    tspan = [0 200];

    % Waktu tetap untuk interpolasi hasil
    t_fix = linspace(tspan(1), tspan(2), 500);  % sederhananya 200 titik waktu tersebut diperbanyak menjadi 500 titik

    % Kondisi awal
    y0 = [94500000; 47250000; 47250000; 8477640; 423882; 423882; 0];

    % Variasi laju vaksinasi
    p_values = [0, 0.01, 0.1, 0.5, 0.7];
    Ih_all = zeros(length(t_fix), length(p_values));  % Preallocate hasil Ih

    % Warna untuk perbandingan Ih
    colors = lines(length(p_values));

    for i = 1:length(p_values)
        p = p_values(i);
        [t, y] = ode45(@(t, y) odefun(t, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p), tspan, y0);

        % Interpolasi Ih ke t_fix
        Ih_interp = interp1(t, y(:, 6), t_fix, 'linear');

        % Simpan hasil Ih
        Ih_all(:, i) = Ih_interp;

        % Tampilkan dalam figure terpisah
        figure;
        plot(t, y, 'LineWidth', 1.5);
        legend('S_v','E_v','I_v','S_h','E_h','I_h','R_h');
        xlabel('Waktu (hari)');
        ylabel('Jumlah individu');
        title(['Simulasi SEIR DBD dengan p = ' num2str(p)]);
        grid on;
    end

    % Plot perbandingan Ih
    figure;
    hold on;
    for i = 1:length(p_values)
        plot(t_fix, Ih_all(:, i), 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', ['p = ' num2str(p_values(i))]);
    end
    hold off;
    legend;
    xlabel('Waktu (hari)');
    ylabel('Jumlah Ih (Terinfeksi)');
    title('Perbandingan Kurva Ih pada Berbagai Nilai p');
    grid on;
end

function dydt = odefun(~, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p)
    Sv = y(1);
    Ev = y(2);
    Iv = y(3);
    Sh = y(4);
    Eh = y(5);
    Ih = y(6);
    Rh = y(7);

    dSv = mu_v * Nv - (eta * Sv * Ih / Nh) - mu_v * Sv; 
    dEv = (eta * Sv * Ih / Nh) - theta_ev * Ev - mu_v * Ev; 
    dIv = theta_ev * Ev - mu_v * Iv; 

    dSh = mu_h * Nh - (beta * Sh * Iv / Nh) - p * Sh - mu_h * Sh;
    dEh = (beta * Sh * Iv / Nh) - delta_eh * Eh - mu_h * Eh;
    dIh = delta_eh * Eh - gamma_h * Ih - mu_h * Ih;
    dRh = gamma_h * Ih + p * Sh - mu_h * Rh;

    dydt = [dSv; dEv; dIv; dSh; dEh; dIh; dRh];
end
