/*import the NLRB data from the csv file*/

data mydata;
    infile 'C:\Users\flash\Downloads\206 TRG Exam\NLRB_data.csv' dlm=',' firstobs=2; /* Specify the delimiter and skip the header row */
    input x y; 
run;

/*fit and LPM Model , for binary outcome variable y*/

proc reg data=mydata;
    model y = x / hccme=0;

    output out=lmp_results predicted = prob_hat; /* Store predicted probabilities in a new dataset */

run;
quit;

proc iml ;
    cal = 0.2903 + 0.01295 * 20.39;
    print cal;
quit;


/*now we are considering model 2*/

proc iml;
    use mydata;
    read all var {x y} into data;
    close mydata;
    
    n=nrow(data);
    x = data[,1];       
    y = data[,2];
    
    /*define log-lkelihood function*/

    start calc_ll(y,x,beta,lambda);
        term1 =1 - beta*exp(-lambda*x);
        term2 = beta*exp(-lambda*x);
        
        term1 = max(term1, 1e-10); /* Avoid log(0) by ensuring term1 is not too small */
        term2 = max(term2, 1e-10); /* Avoid log(0) by ensuring term2 is not too small */
        ll_vec = y#log(term1) + (1-y)#log(term2);
        return(sum(ll_vec));

    finish;


    /*define gradient function*/
    start calc_grad(y, x, beta, lambda);
        exp_term = exp(-lambda * x);
        denom = 1 - beta * exp_term;
        
        /* Using element-wise multiplication (#) for vector operations */
        grad_beta = sum( (1 - y)/beta - (y # exp_term) / denom );
        grad_lambda = sum( (y # beta # x # exp_term) / denom - (1 - y) # x );
        
        return (grad_beta // grad_lambda); /* Return as 2x1 column vector */
    finish;

    /* 4. Define the Hessian Matrix */
    start calc_hess(y, x, beta, lambda);
        exp_term = exp(-lambda * x);
        exp_term2 = exp(-2 * lambda * x);
        denom2 = (1 - beta * exp_term)##2;
        
        h11 = sum( -(1 - y)/(beta**2) - (y # exp_term2) / denom2 );
        h22 = sum( -(y # beta # (x##2) # exp_term) / denom2 );
        h12 = sum( (y # x # exp_term) / denom2 );
        
        /* Construct 2x2 matrix */
        hess = (h11 || h12) // 
               (h12 || h22);
        return (hess);
    finish;

    /* 5. Initialize Parameters and Newton-Raphson Settings */
    beta = 0.5;
    lambda = 0.01;
    theta = beta // lambda;
    
    tol = 1e-8;      /* Convergence tolerance */
    max_iter = 100;  /* Maximum allowable iterations to prevent infinite loops */
    diff = 1;
    iter = 0;        /* Initialize iteration count */

    print "--- Newton-Raphson Iteration History ---";
    
    /* 6. Execute Newton-Raphson Optimization Loop */
    do while (diff > tol & iter < max_iter);
        /* Calculate current Log-Likelihood */
        ll = calc_ll(y, x, beta, lambda);
        
        /* Print current iteration state */
        Iteration = iter;
        Current_Beta = beta;
        Current_Lambda = lambda;
        Log_Likelihood = ll;
        print Iteration Current_Beta Current_Lambda Log_Likelihood;
        
        /* Compute Gradient and Hessian */
        grad = calc_grad(y, x, beta, lambda);
        hess = calc_hess(y, x, beta, lambda);
        
        /* Newton-Raphson Update Rule: theta^(t+1) = theta^(t) - H^(-1) * grad */
        theta_new = theta - inv(hess) * grad;
        
        /* Calculate the maximum absolute difference between old and new parameters */
        diff = max(abs(theta_new - theta));
        
        /* Update parameters for the next iteration */
        theta = theta_new;
        beta = theta[1];
        lambda = theta[2];
        iter = iter + 1;
    end;

    /* 7. Final Output Presentation */
    ll_final = calc_ll(y, x, beta, lambda);
    
    print "--- Optimization Converged ---";
    Total_Iterations = iter;
    print Total_Iterations diff;
    
    Final_Estimates = theta;
    param_names = {"Beta", "Lambda"};
    print Final_Estimates[rowname=param_names label="Final Parameter Estimates"];
    
    print ll_final[label="Final Maximized Log-Likelihood"];
quit;



proc iml;
    /* 1. Load Data */
    use mydata;
    read all var {y} into y;
    read all var {x} into x;
    close mydata;

    n = nrow(y);
    B = 500; /* Number of bootstrap resamples */
    
    /* Matrix to store the 6 target estimates for each of the 500 iterations */
    boot_results = J(B, 6, .); 
    cnames = {"b0", "b1", "pred1", "beta", "lambda", "pred2"};
    
    call randseed(12345); /* Set seed for reproducibility */
    
    /* Bootstrap Loop */
    do i = 1 to B;
    
        /* --- A. Generate Bootstrap Sample --- */
        u = J(n, 1, 0);
        call randgen(u, "Uniform");
        idx = ceil(n * u); /* Generates indices 1 to n with replacement */
        y_b = y[idx];
        x_b = x[idx];
        
        /* --- B. Model 1 (LPM) Estimation via OLS Matrix Algebra --- */
        X_des = J(n, 1, 1) || x_b; /* Design matrix with intercept column */
        beta_ols = inv(X_des` * X_des) * X_des` * y_b;
        
        b0 = beta_ols[1];
        b1 = beta_ols[2];
        pred1 = b0 + b1 * 20.39; /* Predict for x=20.39 */
        
        /* --- C. Model 2 (Non-linear) Estimation via Newton-Raphson --- */
        beta_nl = 0.5; 
        lambda_nl = 0.01;
        theta = beta_nl // lambda_nl;
        tol = 1e-8; 
        diff = 1; 
        iter = 0; 
        max_iter = 50; /* Limit iterations in case a bootstrap sample is ill-conditioned */
        
        do while (diff > tol & iter < max_iter);
            exp_term = exp(-lambda_nl * x_b);
            denom = 1 - beta_nl * exp_term;
            
            /* Guard against division by zero in volatile resamples */
            denom = denom <> 1e-12; 
            
            /* Gradient calculation */
            g1 = sum( (1 - y_b)/beta_nl - (y_b # exp_term) / denom );
            g2 = sum( (y_b # beta_nl # x_b # exp_term) / denom - (1 - y_b) # x_b );
            grad = g1 // g2;
            
            /* Hessian calculation */
            exp_term2 = exp(-2 * lambda_nl * x_b);
            denom2 = denom##2;
            h11 = sum( -(1 - y_b)/(beta_nl**2) - (y_b # exp_term2) / denom2 );
            h22 = sum( -(y_b # beta_nl # (x_b##2) # exp_term) / denom2 );
            h12 = sum( (y_b # x_b # exp_term) / denom2 );
            hess = (h11 || h12) // (h12 || h22);
            
            /* Newton-Raphson Update */
            theta_new = theta - inv(hess) * grad;
            diff = max(abs(theta_new - theta));
            
            theta = theta_new;
            beta_nl = theta[1];
            lambda_nl = theta[2];
            iter = iter + 1;
        end;
        
        pred2 = 1 - beta_nl * exp(-lambda_nl * 20.39); /* Predict for x=20.39 */
        
        /* --- D. Store results for iteration i --- */
        boot_results[i, ] = b0 || b1 || pred1 || beta_nl || lambda_nl || pred2;


    end;

    /* --- E. Calculate Mean, 95% CI, and 95% PI for Predictions --- */
    
    /* Extract the vectors of the 500 predicted probabilities for x=20.39 */
    p1_boot = boot_results[, 3]; /* Model 1 predictions */
    p2_boot = boot_results[, 6]; /* Model 2 predictions */

    /* 1. Calculate Bootstrap Means */
    mean_p1 = mean(p1_boot);
    mean_p2 = mean(p2_boot);

    /* 2. Calculate 95% Confidence Intervals (Empirical Percentile Method) */
    /* Sort the predictions from lowest to highest */
    call sort(p1_boot, 1);
    call sort(p2_boot, 1);
    
    /* Find the indices for the 2.5th and 97.5th percentiles */
    idx_lower = round(B * 0.025);
    idx_upper = round(B * 0.975);
    
    ci_p1_lower = p1_boot[idx_lower];
    ci_p1_upper = p1_boot[idx_upper];
    
    ci_p2_lower = p2_boot[idx_lower];
    ci_p2_upper = p2_boot[idx_upper];

    /* 3. Calculate 95% Prediction Intervals (Normal Approximation) */
    /* Total Variance = Var(Estimate) + Inherent Variance of New Observation Y */
    /* For a Bernoulli outcome, the inherent variance is p(1-p) */
    var_p1_hat = var(p1_boot);
    var_p2_hat = var(p2_boot);
    
    se_pi1 = sqrt(var_p1_hat + mean_p1 * (1 - mean_p1));
    se_pi2 = sqrt(var_p2_hat + mean_p2 * (1 - mean_p2));
    
    /* Compute normal approximation bounds for the PI */
    pi_p1_lower = mean_p1 - 1.96 * se_pi1;
    pi_p1_upper = mean_p1 + 1.96 * se_pi1;
    
    pi_p2_lower = mean_p2 - 1.96 * se_pi2;
    pi_p2_upper = mean_p2 + 1.96 * se_pi2;
    
    /* --- Print the Results --- */
    print "--- Bootstrap Inference for Predictions at x = 20.39 ---";
    
    mod1_res = mean_p1 || ci_p1_lower || ci_p1_upper || pi_p1_lower || pi_p1_upper;
    mod2_res = mean_p2 || ci_p2_lower || ci_p2_upper || pi_p2_lower || pi_p2_upper;
    
    summary_matrix = mod1_res // mod2_res;
    r_names = {"Model 1 (LPM)", "Model 2 (Non-linear)"};
    c_names = {"Mean", "CI_Lower", "CI_Upper", "PI_Lower", "PI_Upper"};
    
    print summary_matrix[rowname=r_names colname=c_names format=8.4 label="95% CI and PI for Predictions"];
    

    /* --- F. Calculate Mean and 95% CI for Parameter Estimates --- */
    
    /* Extract the vectors of the 500 parameter estimates */
    b0_boot = boot_results[, 1];
    b1_boot = boot_results[, 2];
    beta_boot = boot_results[, 4];
    lambda_boot = boot_results[, 5];

    /* 1. Calculate Bootstrap Means */
    mean_b0 = mean(b0_boot);
    mean_b1 = mean(b1_boot);
    mean_beta = mean(beta_boot);
    mean_lambda = mean(lambda_boot);

    /* 2. Calculate 95% Confidence Intervals (Empirical Percentile Method) */
    /* Sort each parameter array from lowest to highest */
    call sort(b0_boot, 1);
    call sort(b1_boot, 1);
    call sort(beta_boot, 1);
    call sort(lambda_boot, 1);
    
    /* Ensure percentile indices are set */
    idx_lower = round(B * 0.025);
    idx_upper = round(B * 0.975);
    
    ci_b0_lower = b0_boot[idx_lower];
    ci_b0_upper = b0_boot[idx_upper];
    
    ci_b1_lower = b1_boot[idx_lower];
    ci_b1_upper = b1_boot[idx_upper];
    
    ci_beta_lower = beta_boot[idx_lower];
    ci_beta_upper = beta_boot[idx_upper];
    
    ci_lambda_lower = lambda_boot[idx_lower];
    ci_lambda_upper = lambda_boot[idx_upper];

    /* --- Print the Results --- */
    print "--- Bootstrap Inference for Parameter Estimates ---";
    
    mod1_b0_res = mean_b0 || ci_b0_lower || ci_b0_upper;
    mod1_b1_res = mean_b1 || ci_b1_lower || ci_b1_upper;
    mod2_beta_res = mean_beta || ci_beta_lower || ci_beta_upper;
    mod2_lambda_res = mean_lambda || ci_lambda_lower || ci_lambda_upper;
    
    param_summary = mod1_b0_res // mod1_b1_res // mod2_beta_res // mod2_lambda_res;
    r_names_param = {"Model 1: Intercept (b0)", "Model 1: Income (b1)", "Model 2: Beta", "Model 2: Lambda"};
    c_names_param = {"Mean", "CI_Lower", "CI_Upper"};
    
    print param_summary[rowname=r_names_param colname=c_names_param format=8.4 label="95% CI for Parameters"];
    /* Export bootstrap results matrix to a SAS dataset for plotting */
    create boot_data from boot_results[colname=cnames];
    append from boot_results;
    close boot_data;
quit;


/* 2. Plotting the Bootstrap Distributions */

title "Model 1: Bootstrap Distribution for Intercept (b0)";
proc sgplot data=boot_data; 
    histogram b0; 
    density b0; 
run;

title "Model 1: Bootstrap Distribution for Income Coefficient (b1)";
proc sgplot data=boot_data; 
    histogram b1; 
    density b1; 
run;

title "Model 1: Bootstrap Predicted Probability at x=20.39";
proc sgplot data=boot_data; 
    histogram pred1; 
    density pred1; 
run;

title "Model 2: Bootstrap Distribution for Parameter Beta";
proc sgplot data=boot_data; 
    histogram beta; 
    density beta; 
run;

title "Model 2: Bootstrap Distribution for Parameter Lambda";
proc sgplot data=boot_data; 
    histogram lambda; 
    density lambda; 
run;

title "Model 2: Bootstrap Predicted Probability at x=20.39";
proc sgplot data=boot_data; 
    histogram pred2; 
    density pred2; 
run;
title; /* Clear title */