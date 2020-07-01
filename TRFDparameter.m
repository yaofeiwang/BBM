% -------------------------------------------------------------------------
% BBM | Feb 2020 | version 1.0 | 
% -------------------------------------------------------------------------
% Contact: yaofei@mail.ustc.edu.cn, Yaofei Wang
% -------------------------------------------------------------------------
% Non-Additive Cost Functions for JPEG Steganography Based on Block Boundary Maintenance
% -------------------------------------------------------------------------

function [D4, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6] = TRFDparameter(q_tab, W, L)

    T = ones(1, 8);
    T(1) = 2;
    T = 4 * T;

    wr = ones(8, 8);
    wr(1, :) = 1 ./ sqrt(2);
    wc = wr';
    R = q_tab .* wr .* wc ./ 2;

    F = zeros(8, 8);

    for i = 1:8

        for j = 1:8
            F(i, j) = cos(pi .* (i - 1) ./ 16) .* cos(pi .* (j - 1) ./ 16) + cos(15 .* pi .* (i - 1) ./ 16) .* cos(15 .* pi .* (j - 1) ./ 16);
        end

    end

    F_hor2 = circshift(F, -2, 2);
    F_hor4 = circshift(F, -4, 2);
    F_hor6 = circshift(F, -6, 2);

    F_ver2 = circshift(F, -2, 1);
    F_ver4 = circshift(F, -4, 1);
    F_ver6 = circshift(F, -6, 1);

    for i = 1:8
        F_hor2(:, i) = F_hor2(i, i);
        F_hor4(:, i) = F_hor4(i, i);
        F_hor6(:, i) = F_hor6(i, i);
        F_ver2(i, :) = F_ver2(i, i);
        F_ver4(i, :) = F_ver4(i, i);
        F_ver6(i, :) = F_ver6(i, i);
    end

    D4 = zeros(8, 8);

    for i = 1:8

        for j = 1:8
            D4(i, j) = R(i, j) .* R(i, j) .* (F(i, i) .* T(j) + F(j, j) .* T(i));
        end

    end

    T_hor = zeros(8, 8);

    for i = 1:8
        T_hor(i, :) = T(i);
    end

    T_ver = zeros(8, 8);

    for i = 1:8
        T_ver(:, i) = T(i);
    end

    T_hor = repmat(T_hor, W ./ 8, L ./ 8);
    T_ver = repmat(T_ver, W ./ 8, L ./ 8);
    R = repmat(R, W ./ 8, L ./ 8);
    F_hor2 = repmat(F_hor2, W ./ 8, L ./ 8);
    F_hor4 = repmat(F_hor4, W ./ 8, L ./ 8);
    F_hor6 = repmat(F_hor6, W ./ 8, L ./ 8);
    F_ver2 = repmat(F_ver2, W ./ 8, L ./ 8);
    F_ver4 = repmat(F_ver4, W ./ 8, L ./ 8);
    F_ver6 = repmat(F_ver6, W ./ 8, L ./ 8);
    D4 = repmat(D4, W ./ 8, L ./ 8);

end
