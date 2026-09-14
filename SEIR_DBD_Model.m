function SEIR_DBD_Model
    % Parameter
    beta = 0.16667;     % Tingkat penularan manusia oleh nyamuk
    eta = 0.1428;       % Tingkat penularan nyamuk oleh manusia
    delta_eh = 0.375;   % Laju paparan pada manusia (E_h ke I_h)
    theta_ev = 0.75;    % Laju paparan pada nyamuk (E_v ke I_v)
    gamma_h = 0.328833; % Laju pemulihan manusia
    mu_h = 0.0000002;   % Angka kematian manusia
    mu_v = 0.0323;      % Angka kematian nyamuk
    Nh = 9419600;       % Populasi manusia total
    Nv = 5000;          % Populasi nyamuk total
    p = 0.7;           % Laju vaksinasi

    % Waktu simulasi
    tspan = [0 200];

    % Kondisi awal [Sv, Ev, Iv, Sh, Eh, Ih, Rh]
    y0 = [94500000; 47250000; 47250000; 8477640; 423882; 423882; 0];

    % Pemanggilan ODE
    [t, y] = ode45(@(t, y) odefun(t, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p), tspan, y0);

    % Plot hasil
    figure;
    plot(t, y);
    legend('S_v','E_v','I_v','S_h','E_h','I_h','R_h');
    xlabel('Waktu (hari)');
    ylabel('Jumlah individu');
    title('Model SEIR DBD dengan Vaksinasi');
end

function dydt = odefun(~, y, beta, eta, delta_eh, theta_ev, gamma_h, mu_h, mu_v, Nh, Nv, p)
    Sv = y(1);
    Ev = y(2);
    Iv = y(3);
    Sh = y(4);
    Eh = y(5);
    Ih = y(6);
    Rh = y(7);

    % Persamaan diferensial
    dSv = mu_v * Nv - (eta * Sv * Ih / Nh) - mu_v * Sv; 
    dEv = (eta * Sv * Ih / Nh) - theta_ev * Ev - mu_v * Ev; 
    dIv = theta_ev * Ev - mu_v * Iv; 

    dSh = mu_h * Nh - (beta * Sh * Iv / Nh) - p * Sh - mu_h * Sh;
    dEh = (beta * Sh * Iv / Nh) - delta_eh * Eh - mu_h * Eh;
    dIh = delta_eh * Eh - gamma_h * Ih - mu_h * Ih;
    dRh = gamma_h * Ih + p * Sh - mu_h * Rh;

    dydt = [dSv; dEv; dIv; dSh; dEh; dIh; dRh];
end
