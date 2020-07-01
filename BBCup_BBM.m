% -------------------------------------------------------------------------
% BBCup_BBM | Feb 2020 | version 1.0 | 
% -------------------------------------------------------------------------
% Contact: yaofei@mail.ustc.edu.cn, Yaofei Wang
% -------------------------------------------------------------------------
% Non-Additive Cost Functions for JPEG Steganography Based on Block Boundary Maintenance
% -------------------------------------------------------------------------

function BBCup_BBM(payload, cover_img, stego_img, DP, D4, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6, embed_order, costfun)

    default_gray_jpeg_obj = jpeg_read(cover_img);

    H = [0 1 0; 1 0 1; 0 1 0];
    H_hor = zeros(1, 17);
    H_hor([1, end]) = 1;
    H_ver = H_hor';

    C_STRUCT = jpeg_read(cover_img);
    dct_coef = C_STRUCT.coef_arrays{1};
    dct_coef_cover = dct_coef;
    q_tab = C_STRUCT.quant_tables{1};

    [X, Y] = size(dct_coef);
    q_tab1 = repmat(q_tab, X / 8, Y / 8);

    cost_update_choice = 0;

    if costfun == "UERD"
        [rhoP1_ini, rhoM1_ini] = UERD_Cost(cover_img, cost_update_choice);
    elseif costfun == "UNIWARD"
        [rhoP1_ini, rhoM1_ini] = J_UNIWARD_COST(cover_img, cost_update_choice);
    elseif costfun == "HDS"
        [rhoP1_ini, rhoM1_ini] = HDS_Cost(cover_img, cost_update_choice);
    elseif costfun == "RBV"
        [rhoP1_ini, rhoM1_ini] = RBV_Cost(cover_img, cost_update_choice);
    elseif costfun == "BET_HILL"
        [rhoP1_ini, rhoM1_ini] = BET_HILL(cover_img, cost_update_choice);
    end

    [rhoP1_ini4(:, :, 1), rhoP1_ini4(:, :, 2), rhoP1_ini4(:, :, 3), rhoP1_ini4(:, :, 4)] = block8to4(rhoP1_ini);
    [rhoM1_ini4(:, :, 1), rhoM1_ini4(:, :, 2), rhoM1_ini4(:, :, 3), rhoM1_ini4(:, :, 4)] = block8to4(rhoM1_ini);

    [W, L] = size(rhoM1_ini4(:, :, 1));

    for i = 1:numel(DP)
        DP{i} = repmat(DP{i}, W / 8, L / 8);
    end

    S_COEFFS = dct_coef_cover;

    [dct_coef_cover4(:, :, 1), dct_coef_cover4(:, :, 2), dct_coef_cover4(:, :, 3), dct_coef_cover4(:, :, 4)] = block8to4(dct_coef_cover);

    for i = 1:4
        nzAC = nnz(dct_coef_cover4(:, :, i)) - nnz(dct_coef_cover(1:8:end, 1:8:end)) ./ 4;
        message(i) = round(payload * nzAC);
    end

    for i = 1:4
        BBC_DP = zeros(16, 16);

        if i == 1
            BBC_DP(1:8, 1:8) = 1;
        elseif i == 2
            BBC_DP(1:8, 9:16) = 1;
        elseif i == 3
            BBC_DP(9:16, 9:16) = 1;
        elseif i == 4
            BBC_DP(9:16, 1:8) = 1;
        end

        BBC_DP = logical(repmat(BBC_DP, X ./ 16, Y ./ 16));

        diff = S_COEFFS - dct_coef_cover;
        diff_hor = imfilter(diff, H_hor);
        diff_ver = imfilter(diff, H_ver);
        diff_hor(:, 2:2:end) = -1 * diff_hor(:, 2:2:end);
        diff_ver(2:2:end, :) = -1 * diff_ver(2:2:end, :);
        diff = diff_hor + diff_ver;
        [rhoP1, rhoM1] = deal(rhoP1_ini, rhoM1_ini);

        rhoP1(diff > 0) = rhoP1(diff > 0) .* 0.5;
        rhoM1(diff < 0) = rhoM1(diff < 0) .* 0.5;

        [rhoP1_ini4(:, :, 1), rhoP1_ini4(:, :, 2), rhoP1_ini4(:, :, 3), rhoP1_ini4(:, :, 4)] = block8to4(rhoP1);
        [rhoM1_ini4(:, :, 1), rhoM1_ini4(:, :, 2), rhoM1_ini4(:, :, 3), rhoM1_ini4(:, :, 4)] = block8to4(rhoM1);
        [S_COEFFS4(:, :, 1), S_COEFFS4(:, :, 2), S_COEFFS4(:, :, 3), S_COEFFS4(:, :, 4)] = block8to4(S_COEFFS);
        [dct_coef_cover4(:, :, 1), dct_coef_cover4(:, :, 2), dct_coef_cover4(:, :, 3), dct_coef_cover4(:, :, 4)] = block8to4(dct_coef_cover);

        S_COEFFS4(:, :, i) = BBMpart(dct_coef_cover4(:, :, i), rhoP1_ini4(:, :, i), rhoM1_ini4(:, :, i), message(i), DP, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6, embed_order);
        S_COEFFS = block4to1(S_COEFFS4(:, :, 1), S_COEFFS4(:, :, 2), S_COEFFS4(:, :, 3), S_COEFFS4(:, :, 4));

    end

    S_STRUCT = default_gray_jpeg_obj;
    S_STRUCT.coef_arrays{1} = S_COEFFS;
    jpeg_write(S_STRUCT, stego_img);

end
