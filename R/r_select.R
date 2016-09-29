### Automatic selection for r (based on the moments estimator for the extremal index
### by Ferros and Segers)
r.select = function(x){

    # Compute the estimate from Ferros and Segers
    theta = rep(list(NULL), NCOL(x$varmat))
    for (j in 1:length(theta)){
        ind = which(x$varmat[,j] > x$threshold)
        vec = as.numeric(names(ind))

        tmp1 = c(1, which(diff(vec) > x$days.in.months)+1)
        tmp2 = c(which(diff(vec) > x$days.in.months), length(vec))

        for (i in 1:length(tmp1)){
            tpart = diff(vec[tmp1[i]:tmp2[i]])
            if (length(tpart) > 0){
                if (max(tpart) <= 2){
                    theta[[j]] = c(theta[[j]], min(1,
                        (2*sum(tpart)^2) / (sum(tpart^2)*length(tpart))))
                } else {
                    theta[[j]] = c(theta[[j]], min(1,
                        (2*sum(tpart-1)^2) / (sum((tpart-1)*(tpart-2))*length(tpart))))
                    }
            } else {
                theta[[j]] = c(theta[[j]], 1)
                }
            }
        }

    # Given an r, compute the mle of theta (n_c / n_u)
    # this tests r=0,...,14, make length(theta.r) bigger to test more
    # theta.r is non-increasing with r
    theta.r = double(15)
    r = seq(0, length.out = length(theta.r))
    for (i in 1:length(r)){
        n_c = 0
        n_u = 0
        for (j in 1:length(x$run.index)){
            tmp = decluster(x = x$varmat[,x$run.index[j]],
                times = x$time.dates, u = x$threshold, r = r[i], doplot = FALSE)
            n_c = n_c + tmp$n_clust
            n_u = n_u + tmp$n_u
            }
        theta.r[i] = n_c / n_u
        }

    return (r[which.min(abs(theta.r - mean(unlist(theta))))])
    }

r_select_hier = function(x, p){

    # Compute the estimate from Ferros and Segers
    theta = NULL
    ind = which(x$varmat[,p] > x$threshold[p])
    vec = as.numeric(names(ind))

    tmp1 = c(1, which(diff(vec) > x$days.in.months)+1)
    tmp2 = c(which(diff(vec) > x$days.in.months), length(vec))

    for (i in 1:length(tmp1)){
        tpart = diff(vec[tmp1[i]:tmp2[i]])
        if (length(tpart) > 0){
            if (max(tpart) <= 2){
                theta = c(theta, min(1,
                    (2*sum(tpart)^2) / (sum(tpart^2)*length(tpart))))
            } else {
                theta = c(theta, min(1,
                    (2*sum(tpart-1)^2) / (sum((tpart-1)*(tpart-2))*length(tpart))))
                }
        } else {
            theta = c(theta, 1)
            }
        }

    # Given an r, compute the mle of theta (n_c / n_u)
    # this tests r=0,...,14, make length(theta.r) bigger to test more
    # theta.r is non-increasing with r
    theta.r = double(15)
    r = seq(0, length.out = length(theta.r))
    for (i in 1:length(r)){
        tmp = decluster(x = x$varmat[,p],
            times = x$time.dates, u = x$threshold[p], r = r[i], doplot = FALSE)
        theta.r[i] = tmp$n_clust / tmp$n_u
        }

    return (r[which.min(abs(theta.r - mean(unlist(theta))))])
    }
