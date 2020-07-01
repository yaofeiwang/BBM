% -------------------------------------------------------------------------
% BBM | Feb 2020 | version 1.0 |
% -------------------------------------------------------------------------
% Contact: yaofei@mail.ustc.edu.cn, Yaofei Wang
% -------------------------------------------------------------------------
% Non-Additive Cost Functions for JPEG Steganography Based on Block Boundary Maintenance
% -------------------------------------------------------------------------

function S_COEFFS = BBMpart(dct_coef_cover, rhoP1_ini, rhoM1_ini, message, DP, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6, order)
    All_embed_order = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

    left2 = @(x)(circshift(x.data, -2, 2));
    left4 = @(x)(circshift(x.data, -4, 2));
    left6 = @(x)(circshift(x.data, -6, 2));
    up2 = @(x)(circshift(x.data, -2, 1));
    up4 = @(x)(circshift(x.data, -4, 1));
    up6 = @(x)(circshift(x.data, -6, 1));

    S_COEFFS = dct_coef_cover;

    Message_part = BlockMessage(rhoP1_ini, rhoM1_ini, message, DP);

    S_COEFFS(logical(DP{All_embed_order(order, 1)})) = DCTEmbeddingSimulator(dct_coef_cover(logical(DP{All_embed_order(order, 1)})), rhoP1_ini(logical(DP{All_embed_order(order, 1)})), rhoM1_ini(logical(DP{All_embed_order(order, 1)})), Message_part(All_embed_order(order, 1)), All_embed_order(order, 1));

    for i = 2:numel(DP)

        [rhoP1, rhoM1] = deal(rhoP1_ini, rhoM1_ini);
        diff = S_COEFFS - dct_coef_cover;
        diff = diff .* R;
        diff_hor2 = blockproc(diff, [8, 8], left2);
        diff_hor4 = blockproc(diff, [8, 8], left4);
        diff_hor6 = blockproc(diff, [8, 8], left6);
        diff_ver2 = blockproc(diff, [8, 8], up2);
        diff_ver4 = blockproc(diff, [8, 8], up4);
        diff_ver6 = blockproc(diff, [8, 8], up6);
        g_fun = -R .* ((diff_hor2 .* F_hor2 + diff_hor4 .* F_hor4 + diff_hor6 .* F_hor6) .* T_hor + (diff_ver2 .* F_ver2 + diff_ver4 .* F_ver4 + diff_ver6 .* F_ver6) .* T_ver);
        f_fun = 1 ./ ((exp(1)).^abs(g_fun));
        rhoP1(g_fun > 0) = rhoP1(g_fun > 0) .* f_fun(g_fun > 0);
        rhoM1(g_fun < 0) = rhoM1(g_fun < 0) .* f_fun(g_fun < 0);

        S_COEFFS(logical(DP{All_embed_order(order, i)})) = DCTEmbeddingSimulator(dct_coef_cover(logical(DP{All_embed_order(order, i)})), rhoP1(logical(DP{All_embed_order(order, i)})), rhoM1(logical(DP{All_embed_order(order, i)})), Message_part(All_embed_order(order, i)), All_embed_order(order, i));
    end

end
