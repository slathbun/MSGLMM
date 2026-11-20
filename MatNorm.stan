// Matrix normal priors
// define the data. 
// z: response variables 
// x: regression covariates
// loc: x/y location coordinates
// NObs: number of observations -- this is the number of study participants
// NSite: number of sites used to construct the covar matrix. NSite <= NObs.
data 
{
  int<lower=1> NObs;              /* Number of observations */
  int<lower=1> NSite;             /* Number of sites (unique coordinates). NSite <= NObs */
  int<lower=1> P;                 /* number response variables */ 
  int<lower=1> K;                 /* number of coefficients (including intercept) */
  int<lower=1> Q;                 /* number of latent random fields */

  matrix[NObs, K] X;              /* Design matrix */ 
  array[NSite] vector[2] Loc;     /* Coordinate Pairs */  
  array[NObs] int<lower=0> Z1;    /* response 1: MOCA (Binomial, n=30) */
  array[NObs] int<lower=0> Z2;    /* response 2: RecallDelay (Binomial, n=10) */
  array[NObs] int<lower=0> Z3;    /* response 3: Memory Rey Figure */
  array[NObs] int<lower=0> Z4;    /* response 4: BNT */
  array[NObs] int<lower=0> Z5;    /* response 5: FlFastTotal (Poisson) */ 
  array[NObs] real Z6;            /* response 6: Log Executive Task Time */
  array[NObs] int<lower=0> Z7;    /* response 7: Clock Drawing */
  array[NObs] int<lower=0> Z8;    /* response 8: Logical Memory */
  array[NObs] int<lower=0> Site;  /* Lookup table: Loc[i] is site of observation i */

}

transformed data 
{
  /* takes NSite x 1 vector to NObs x 1 vector */
  /* e.g. Site2Obs * epsilon --> epsilon_per_obs */
  matrix[NObs,NSite] SiteToObs; 

  for (i in 1:NObs) 
  {
    for (j in 1:NSite) 
    {
      SiteToObs[i, j] = 0.0;
    }
  }

  for (i in 1:NObs) 
  {
    SiteToObs[i, Site[i]] = 1.0;
  }
}

parameters 
{
  array[Q] real<lower=0> diff;  /* range parameter */
  cov_matrix[P] phi;
  matrix[P,Q] zeta;
  array[P] vector[K] betas;      /* regression coefficients */
  array[Q] vector[NSite] w;
  real<lower=0> vardelta;

}

transformed parameters
{
  array[Q] real<lower=0> range;
  range[1]=diff[1];
  matrix[P,Q] A;     /* coregionalization matrix */
  matrix[P,P] L_phi;
  for(i in 2:Q)
     {
         range[i]=range[i-1]+diff[i];
     }
  L_phi=cholesky_decompose(phi);
  A=L_phi*zeta;
}

model
{
  array[P] vector[NObs] eta;           /* canonical params */
  array[Q] matrix[NSite,NSite] Psi;    /* covariance matrices of latent variables \xi */
  array[Q] matrix[NSite,NSite] L_Psi;  /* Cholesky decomoposition of Psi */
  array[Q] vector[NSite] xi;           /* latent random field */
  array[P] vector[NSite] epsilon;      /* latent random field */
  
  for(k in 1:Q)
  {
       Psi[k] = gp_exponential_cov(Loc, 1.0, range[k]/3.0);
       L_Psi[k]=cholesky_decompose(Psi[k]);
       xi[k]=L_Psi[k]*w[k];
  }

  for(i in 1:P)
  {
     for(j in 1:NSite)
     {
         epsilon[i][j]=0.0;
         for(k in 1:Q)
         {
             epsilon[i][j]=epsilon[i][j]+A[i,k]*xi[k][j];
         }
      }
  }
 
  for(i in 1:P)
  {
      eta[i] =  X*betas[i] + SiteToObs*epsilon[i]; 
  }


 /*    priors    */

  w ~ multi_normal(rep_vector(0,NSite), identity_matrix(NSite));
  diff~inv_gamma(2.0,10.0);
  phi~inv_wishart(P+4,2.0*identity_matrix(P));
  to_vector(zeta)~normal(0.0,0.5);
  betas ~ multi_normal(rep_vector(0,K), 10 * identity_matrix(K));
  vardelta~inv_gamma(3.0,2.0);
  Z1 ~ binomial_logit(30, eta[1]);
  Z2 ~ binomial_logit(10, eta[2]);
  Z3 ~ binomial_logit(24, eta[3]);
  Z4 ~ binomial_logit(30, eta[4]);
  Z5 ~ poisson_log(eta[5]);
  Z6 ~ normal(eta[6],vardelta);
  Z7 ~ binomial_logit(15, eta[7]);
  Z8 ~ binomial_logit(25, eta[8]);
}

generated quantities
{
  matrix[P,P] T;
  array[P] vector[K] expbeta;
  matrix[P,P] corr;
  T=A*A';
   for(i in 1:P)
     {
      for(j in 1:K)
        {
         expbeta[i][j]=exp(betas[i][j]);
        }
      for(j in 1:P)
        {
         corr[i,j]=T[i,j]/sqrt(T[i,i]*T[j,j]);
        }
     }
 }





