% -------------------------------------------------------------------------
% BBM | Feb 2020 | version 1.0 | 
% -------------------------------------------------------------------------
% Contact: yaofei@mail.ustc.edu.cn, Yaofei Wang
% -------------------------------------------------------------------------
% Non-Additive Cost Functions for JPEG Steganography Based on Block Boundary Maintenance
% -------------------------------------------------------------------------

function BBM_DP_scale_order(payload, cover_img, stego_img, DP1, T_hor, T_ver, R, F_hor2, F_hor4, F_hor6, F_ver2, F_ver4, F_ver6, costfun, Scale_fun, embed_order)

    left2 = @(x)(circshift(x.data, -2, 2));
    left4 = @(x)(circshift(x.data, -4, 2));
    left6 = @(x)(circshift(x.data, -6, 2));
    up2 = @(x)(circshift(x.data, -2, 1));
    up4 = @(x)(circshift(x.data, -4, 1));
    up6 = @(x)(circshift(x.data, -6, 1));
    All_embed_order = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

    default_gray_jpeg_obj = jpeg_read(cover_img);
    C_STRUCT = jpeg_read(cover_img);
    dct_coef = C_STRUCT.coef_arrays{1};
    dct_coef_cover = dct_coef;

    S_COEFFS = dct_coef;
    nzAC = nnz(dct_coef_cover) - nnz(dct_coef_cover(1:8:end, 1:8:end));
    message = round(payload * nzAC);

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

    Message_part = BlockMessage(rhoP1_ini, rhoM1_ini, message, DP1);

    S_COEFFS(logical(DP1{All_embed_order(embed_order, 1)})) = DCTEmbeddingSimulator(dct_coef_cover(logical(DP1{All_embed_order(embed_order, 1)})), rhoP1_ini(logical(DP1{All_embed_order(embed_order, 1)})), rhoM1_ini(logical(DP1{All_embed_order(embed_order, 1)})), Message_part(All_embed_order(embed_order, 1)), All_embed_order(embed_order, 1));

    for i = 2:numel(DP1)
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

        switch Scale_fun
            case 'lin'
                f_fun = 1 ./ (1 + abs(g_fun));
            case 'exp'
                f_fun = 1 ./ ((exp(1)).^abs(g_fun));
            case 'log'
                f_fun = 1 ./ (log(exp(1) + abs(g_fun)));
            case 'log10'
                f_fun = 1 ./ (log10(10 + abs(g_fun)));
        end

        rhoP1(g_fun > 0) = rhoP1(g_fun > 0) .* f_fun(g_fun > 0);
        rhoM1(g_fun < 0) = rhoM1(g_fun < 0) .* f_fun(g_fun < 0);
        S_COEFFS(logical(DP1{All_embed_order(embed_order, i)})) = DCTEmbeddingSimulator(dct_coef_cover(logical(DP1{All_embed_order(embed_order, i)})), rhoP1(logical(DP1{All_embed_order(embed_order, i)})), rhoM1(logical(DP1{All_embed_order(embed_order, i)})), Message_part(All_embed_order(embed_order, i)), All_embed_order(embed_order, i));
    end

    S_STRUCT = default_gray_jpeg_obj;
    S_STRUCT.coef_arrays{1} = S_COEFFS;
    jpeg_write(S_STRUCT, stego_img);

end
