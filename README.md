To run estimation, run the script “RunEstimationVN.m”,

Different models an be estimated via the following:


### MNLogit
44: opts.Models={‘MNP’};
46: opts.Prob={‘Logit’};

Output: estimate for sigma (variance)

### MNProbit (Independent Errors)
44: opts.Models={‘MNP’};
46: opts.Prob={‘Probit’};

Output: estimate for sigma (variance)

### DivNorm (simple normalization with normal errors, this is the original Louie et al. model)
44: opts.Models={‘DN’};
46: opts.Prob={‘Probit’};

Output: estimate for sigma and omega

### DivNorm (with beta controlling curvature of normalization)
44: opts.Models={‘DNb’};
46: opts.Prob={‘Probit’};

Output: estimate for sigma and omega and beta

If you want to estimate the model with beta in it (curvature of normalization) I recommend making omega hierarchical (see below).

### DivNorm (with asymmetric recurrent weight)
44: opts.Models={‘DNw’};
46: opts.Prob={‘Probit’};

This was actually our best-fitting model in the paper.

Output: sigma, omega, omega_i-omega (so you have to sum the last two to get the estimate of omega_i)

### Hierarchical Model
line 49: To turn on hierarchical estimation for any parameter, put a 1 in the appropriate spot.

The implementation of this code is a bit rough. For instance, the code is not tested for parameters other than omega, but should work. The output order of the the first K parameters is the same as for the non-hierarchical model. Then for each hierarchical parameter, the second “shape” parameter follows in the same order. So for a hierarchical omega the output will be:

Output: sigma, omega (scale), beta, omega (shape)

This will also take a looong time to run because of the additional integration step needed. But it’s worth it given its ability to address some of the heterogeneity in the data.

### Additional Notes
Though the code will give you standard errors for all parameters, I highly recommend using Likelihood Ratio tests for nested models when possible (i.e. the test on omega is a LR stat distributed chi-square with 1 degree of freedom).

line 69: theta0=[]; %will start the algorithm at random initial points. You can also set them manually.