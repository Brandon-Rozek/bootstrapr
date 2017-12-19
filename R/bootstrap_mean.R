bootstrap.mean.test = function(x, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, conf.level = 0.95, bootstrap_samples = 1000000) {
    ## Select the alternative chosen
    alternative = match.arg(alternative)

    ## Make sure argument mu is correct
    if(!missing(mu) && (length(mu) != 1 || is.na(mu))) {
        stop("'mu' must be a single number")
    }

    ## Make sure confidence interval is between 0 and 1
    if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1)) {
        stop("'conf.level' must be a single number between 0 and 1")
    }

    if (is.null(y)) {
        dname = deparse(substitute(x))
    } else {
        dname = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }

    ## Take out missing values
    x = x[!is.na(x)]
    if (!is.null(y)) {
        y = y[!is.na(y)]
    }

    ## Check to see if we have enough observations and enough variation in data
    nx = length(x)
    mx = mean(x)
    vx = var(x)
   if (is.null(y)) {
        if (nx < 2) {
            stop("not enough 'x' observations")
        }
	    stderr = sqrt(vx / nx)
        if (stderr < 10 * .Machine$double.eps * abs(mx)) {
            stop("data are essentially constant")
        }
    } else {
	    ny = length(y)
        if(nx < 2) {
            stop("not enough 'x' observations")
        }
	    if(ny < 2) {
            stop("not enough 'y' observations")
        }
	    my = mean(y)
	    vy = var(y)
	    stderrx = sqrt(vx / nx)
	    stderry = sqrt(vy / ny)
	    stderr = sqrt(stderrx^2 + stderry^2)
        if(stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) {
            stop("data are essentially constant")
        }
    }


    ## Preallocate the vectors
    bootstrap_replicants = numeric(bootstrap_samples)
    bootstrap_means = numeric(bootstrap_samples)
    ## If it's a one sample bootstrap test
    if (is.null(y)) {
        observed_difference = mx - mu
        ## Simulate null hypothesis by shifting x to have a new mean of mu
        x_shifted = x - mx + mu
        for (i in 1:bootstrap_samples) {
            ## Resample under the null hypothesis the difference between the means of mu
            bootstrap_xshifted = sample(x_shifted, nx, replace = T)
            bootstrap_replicants[i] = mean(bootstrap_xshifted) - mu
            ## Resample original x for confidence interval 
            bootstrap_x = sample(x, nx, replace = T)
            bootstrap_means[i] = mean(bootstrap_x)
        }
        method = "One Sample Bootstrap Test"
        estimate = setNames(mx, "mean of x")
        
    } else { ## If it's a two sample bootstrap test
        observed_difference = mx - my - mu
        ## Under the null hypothesis the difference in means between x and y is mu
        pooled_mean = mean(c(x, y)) + mu
        x_shifted = x - mx + pooled_mean
        y_shifted = y - my + pooled_mean
        for (i in 1:bootstrap_samples) {
            ## Resample under null hypothesis the difference of means
            bootstrap_xshifted = sample(x_shifted, nx, replace = T)
            bootstrap_yshifted = sample(y_shifted, ny, replace = T)
            bootstrap_replicants[i] = mean(bootstrap_xshifted) - mean(bootstrap_yshifted) - mu
            ## Resample original data for confidence interval
            bootstrap_x = sample(x, nx, replace = T)
            bootstrap_y = sample(x, ny, replace = T)
            bootstrap_means[i] = mean(bootstrap_x) - mean(bootstrap_y)
        }
        method = "Two Sample Bootstrap Test"
        estimate = c(mx,my)
        names(estimate) = c("mean of x","mean of y")
    }

    if (alternative == "two.sided") {
        pval = sum(abs(bootstrap_replicants) >= abs(observed_difference)) / bootstrap_samples
        lower_bound = quantile(bootstrap_means, (1 - conf.level) / 2)
        upper_bound = quantile(bootstrap_means, conf.level + (1 - conf.level) / 2)
    } else if (alternative == "less") {
        pval = sum(bootstrap_replicants <= observed_difference) / bootstrap_samples
        lower_bound = -Inf
        upper_bound = quantile(bootstrap_means, conf.level)
    } else if (alternative == "greater") {
        pval = sum(bootstrap_replicants >= observed_difference) / bootstrap_samples
        lower_bound = quantile(bootstrap_means, 1 - conf.level)
        upper_bound = Inf
    }
    cint = c(lower_bound, upper_bound)
    names(mu) = if(!is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") = conf.level
    rval = list(p.value = pval, conf.int = cint, estimate = estimate, null.value = mu,
	       alternative = alternative, method = method, data.name = dname)
    class(rval) = "htest"
    return(rval)
}