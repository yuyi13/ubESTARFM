# This is the R code for the ubESTARFM algorithm
# This is a variant modified from the ESTARFM algorithm developed by Zhu et al. (2010), which was originally written in Python
# Contact: Yi Yu (yi.yu1@anu.edu.au)

ubESTARFM = function(w = 25, DN_min = 250, DN_max = 350, patch_long = 200, tmp_path, out_path, method='zero bias',
                    rst_fine1, rst_fine2, rst_coarse1, rst_coarse2, rst_coarse0){

nl = orig_nl = nrow(rst_fine1); ns = orig_ns = ncol(rst_fine1)

n_nl = as.integer(orig_nl / patch_long); if (n_nl * patch_long < orig_nl) n_nl = n_nl + 1
n_ns = as.integer(orig_ns / patch_long); if (n_ns * patch_long < orig_ns) n_ns = n_ns + 1

# index of patches
ind_patch = matrix(NA, n_nl * n_ns, 4)

# get the edge indices for different patches
for (i_ns in 1:n_ns){
    for (i_nl in 1:n_nl){
        ind_patch[n_ns * (i_nl-1) + i_ns, 1] = (i_ns-1) * patch_long + 1
        ind_patch[n_ns * (i_nl-1) + i_ns, 2] = min(c(ns, i_ns * patch_long))
        ind_patch[n_ns * (i_nl-1) + i_ns, 3] = (i_nl-1) * patch_long + 1
        ind_patch[n_ns * (i_nl-1) + i_ns, 4] = min(c(nl, i_nl * patch_long))
    }
}

# row index of images
row_index = matrix(0, nl, ns)
for (t in 1:nl){
    row_index[t, ] = t
}

# column index of images
col_index = matrix(0, nl, ns)
for (t in 1:ns){
    col_index[, t] = t
}

# asign the whole raster values to these matrices
fine1 = as.matrix(rst_fine1); coarse1 = as.matrix(rst_coarse1)
fine2 = as.matrix(rst_fine2); coarse2 = as.matrix(rst_coarse2)
coarse0 = as.matrix(rst_coarse0)

# find interaction of valid pixels of all input images: exclude missing pixels and background
valid_index = matrix(NA, nl, ns)
ind_valid = which(!is.na(fine1) & !is.na(fine2) & !is.na(coarse1) & !is.na(coarse2) & !is.na(coarse0))
valid_index[ind_valid] = 1  # mark good pixels in all images

# compute the distance of each pixel in the window with the target pixel (integrate window)
D_temp1 = w - t(matrix(rep(0:(w*2),w*2+1), w*2+1, w*2+1)); d1 = D_temp1 ^ 2
D_temp2 = w - matrix(rep(0:(w*2),w*2+1), w*2+1, w*2+1); d2 = D_temp2 ^ 2
D_D_all = 1.0 + sqrt(d1 + d2) / w
D_D_all = as.vector(D_D_all)

print(paste0('there are total ', n_nl*n_ns, ' patches'))

foreach (isub = 1:(n_nl * n_ns), .combine=cbind, .packages=('raster')) %dopar% {

    print(paste0('starting to process the ', isub, 'th patch'))

    col1 = ind_patch[isub, 1]
    col2 = ind_patch[isub, 2]
    row1 = ind_patch[isub, 3]
    row2 = ind_patch[isub, 4]

    # the matries on the prediction date
    fine0 = matrix(NA, nl, ns)
    
    # compute the threshold of similar pixel seeking
    similar_th = c()
    similar_th[1] = sd(fine1[row1:row2, col1:col2], na.rm=TRUE)
    similar_th[2] = sd(fine2[row1:row2, col1:col2], na.rm=TRUE)

    for (j in row1:row2){    # retrieve each target pixel
        for (i in col1:col2){

            if (!is.na(valid_index[j, i])){    # do not process the background

                ai = max(1, i - w)
                bi = min(ns, i + w)
                aj = max(1, j - w)
                bj = min(nl, j + w)

                ind_wind_valid = which(as.vector(valid_index[aj:bj, ai:bi]) == 1)
                position_cand = rep(1, (bi-ai+1)*(bj-aj+1))    # place the location of each similar pixel
                row_wind = row_index[aj:bj, ai:bi]    # a window having row index
                col_wind = col_index[aj:bj, ai:bi]    # a window having col index

                # searching for similar pixels
                for (ipair in 1:2){
                    cand_band = rep(0, (bi-ai+1)*(bj-aj+1))
                    if (ipair == 1){
                        S_S = abs(fine1[aj:bj, ai:bi] - fine1[j, i])
                    } else if (ipair == 2){
                        S_S = abs(fine2[aj:bj, ai:bi] - fine2[j, i])
                    }
                    ind_cand = which(S_S < similar_th[ipair])
                    cand_band[ind_cand] = 1
                    position_cand = position_cand * cand_band
                }

                indcand = which(position_cand != 0 & valid_index[aj:bj, ai:bi] == 1) # the selected pixels' indices (index of candidates)
                number_cand = length(indcand)

                if (number_cand > 5){   # for the case we can select enough similar pixels

                    # compute the correlation
                    S_D_cand = rep(0, number_cand)
                    x_cand = as.vector(col_wind)[indcand]
                    y_cand = as.vector(row_wind)[indcand]
                    finecand = matrix(NA, 2, number_cand)
                    coarsecand = matrix(NA, 2, number_cand)

                    # select the defined similar fine and coarse pixels
                    finecand[1,] = as.vector(fine1[aj:bj, ai:bi])[indcand]
                    finecand[2,] = as.vector(fine2[aj:bj, ai:bi])[indcand]
                    coarsecand[1,] = as.vector(coarse1[aj:bj, ai:bi])[indcand]
                    coarsecand[2,] = as.vector(coarse2[aj:bj, ai:bi])[indcand]

                    if (method == 'zero bias'){
                        # ubESTARFM
                        fine1_pixel_value = fine1[j, i] - mean(finecand[1,]) + mean(coarsecand[1,])
                        finecand[1,] = finecand[1,] - mean(finecand[1,]) + mean(coarsecand[1,]) 

                        fine2_pixel_value = fine2[j, i] - mean(finecand[2,]) + mean(coarsecand[2,])
                        finecand[2,] = finecand[2,] - mean(finecand[2,]) + mean(coarsecand[2,])  

                    } else if (method == 'baseline'){
                        # ESTARFM
                        fine1_pixel_value = fine1[j, i]
                        fine2_pixel_value = fine2[j, i]
                    }

                    S_D_cand = 1.0 - 0.5*(abs((finecand[1,] - coarsecand[1,]) / (finecand[1,] + coarsecand[1,])) +
                                          abs((finecand[2,] - coarsecand[2,]) / (finecand[2,] + coarsecand[2,])))

                    ind_nan = which(S_D_cand > 1.0 | S_D_cand < (-1.0))
                    if (length(ind_nan) > 0){
                        S_D_cand[ind_nan] = 0.5    # correct the abnormal value of correlation
                    }

                    # ********************************** #
                    # spatial distance
                    D_D_cand = rep(NA, number_cand)   
                    if ((bi-ai+1)*(bj-aj+1) < (w*2.0+1)*(w*2.0+1)){   # not an integrate window
                        D_D_cand = 1.0 + sqrt((i-x_cand)^2 + (j-y_cand)^2) / w
                    } else {
                        D_D_cand[1:number_cand] = D_D_all[indcand]      # integrate window
                    }
                    C_D = (1.0 - S_D_cand) * D_D_cand + 0.0000001    # combined distance
                    weight = (1.0/C_D) / sum(1.0/C_D)    # calculate the weight of each similar pixel

                    # ********************************** #
                    # compute V
                    # In this version, we just set the conversion coefficient as 1.0
                    V_cand = 1.0

                    # compute the temporal weight
                    difc_pair1 = abs(mean(as.vector(coarse0[aj:bj, ai:bi])[ind_wind_valid] - as.vector(coarse1[aj:bj, ai:bi])[ind_wind_valid])) + 0.01^5
                    difc_pair2 = abs(mean(as.vector(coarse0[aj:bj, ai:bi])[ind_wind_valid] - as.vector(coarse2[aj:bj, ai:bi])[ind_wind_valid])) + 0.01^5
                    T_weight1 = (1.0/difc_pair1) / (1.0/difc_pair1+1.0/difc_pair2)
                    T_weight2 = (1.0/difc_pair2) / (1.0/difc_pair1+1.0/difc_pair2)

                    # predict from pair1
                    coase0_cand = as.vector(coarse0[aj:bj, ai:bi])[indcand]
                    coase1_cand = as.vector(coarse1[aj:bj, ai:bi])[indcand]
                    fine01 = fine1_pixel_value + sum(weight * V_cand * (coase0_cand-coase1_cand))
                    # predict from pair2
                    coase2_cand = as.vector(coarse2[aj:bj, ai:bi])[indcand]
                    fine02 = fine2_pixel_value + sum(weight * V_cand * (coase0_cand-coase2_cand))
                    # the final prediction
                    fine0[j, i] = T_weight1 * fine01 + T_weight2 * fine02
                    # revise the abnormal prediction
                    if (fine0[j, i] <= DN_min | fine0[j, i] >= DN_max){

                        fine01 = sum(weight * finecand[1,])
                        fine02 = sum(weight * finecand[2,])
                        fine0[j, i] = T_weight1 * fine01 + T_weight2 * fine02
                    }
                } else {   # for the case that no enough similar pixel are selected

                    # then we use all the valid pixels within the window
                    if (method == 'zero bias'){
                        # zero bias
                        fine1_pixel_value = fine1[j, i] - mean((fine1[aj:bj, ai:bi])[ind_wind_valid]) + mean((coarse1[aj:bj, ai:bi])[ind_wind_valid])       
                        fine2_pixel_value = fine2[j, i] - mean((fine2[aj:bj, ai:bi])[ind_wind_valid]) + mean((coarse2[aj:bj, ai:bi])[ind_wind_valid])

                    } else if (method == 'baseline'){
                        fine1_pixel_value = fine1[j, i]
                        fine2_pixel_value = fine2[j, i]
                    }

                    # compute the temporal weight
                    difc_pair1 = mean(as.vector(coarse0[aj:bj, ai:bi])[ind_wind_valid] - as.vector(coarse1[aj:bj, ai:bi])[ind_wind_valid]) + 0.01^5
                    difc_pair1_a = abs(difc_pair1)
                    difc_pair2 = mean(as.vector(coarse0[aj:bj, ai:bi])[ind_wind_valid] - as.vector(coarse2[aj:bj, ai:bi])[ind_wind_valid]) + 0.01^5
                    difc_pair2_a = abs(difc_pair2)
                    T_weight1 = (1.0/difc_pair1_a) / (1.0/difc_pair1_a+1.0/difc_pair2_a)
                    T_weight2 = (1.0/difc_pair2_a) / (1.0/difc_pair1_a+1.0/difc_pair2_a)
                    fine0[j, i] = T_weight1 * (fine1_pixel_value + difc_pair1) + T_weight2 * (fine2_pixel_value + difc_pair2)
                }
            }
        }
    }

    pred_rst = raster(fine0, xmn=xmin(rst_fine1), xmx=xmax(rst_fine1), ymn=ymin(rst_fine1), ymx=ymax(rst_fine1), crs=crs(rst_fine1))
    writeRaster(pred_rst, paste0(tmp_path, 'ubESTARFM_', isub, '.tif'), overwrite=TRUE)

}

print('The processing has been done')

fl = grep(list.files(tmp_path, full.names=TRUE), pattern='.xml', invert=TRUE, value=TRUE)

rst_stk = stack(fl)
final_rst = calc(rst_stk, fun = mean, na.rm = TRUE)

writeRaster(final_rst, out_path, overwrite=TRUE)
file.remove(fl)

return(final_rst)
}