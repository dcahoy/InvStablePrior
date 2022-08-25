#' Bayesian inference for the true rate of double exponential (Laplace) distribution.
#'
#'
#' Generates random numbers from the prior and posterior distributions of the inverse stable-double exponential model. The random variates can be used to simulate prior and posterior predictive distributions as well.
#'
#'
#' @param x vector of double exponential (Laplace) data.
#' @param B test size for the adaptive rejection sampling algorithm.
#' @param alp value between 0 and 1 that controls the shape of the inverse stable prior.
#' @param rho positive value that scales the mean of the inverse stable prior.
#'
#'
#'
#' @return list consisting of the vectors of random numbers from the prior and posterior distributions, the accepted sample size, and the acceptance probability of the adaptive rejection sampling procedure (Algorithm 2 of the first reference below).
#'
#'
#'
#' @references
#'
#' Cahoy and Sedransk (2019). \emph{Inverse stable prior for exponential models.}  Journal of Statistical Theory and Practice, 13, Article 29.
#' <doi:10.1007/s42519-018-0027-2>
#'
#' Meerschaert and Straka (2013). \emph{Inverse stable subordinators.} Math. Model. Nat. Phenom.,  8(2), 1-16.
#' <doi:10.1051/mmnp/20138201>
#'
#' Mainardi, Mura, and Pagnini (2010). \emph{The M-Wright Function in Time-Fractional Diffusion Processes: A Tutorial Survey}. Int. J. Differ. Equ., Volume 2010.
#' <doi:10.1155/2010/104505>
#'
#'
#'
#' @examples
#'
#' alp=0.95
#' require(nimble)
#' dat=rdexp(30, location = 0,  rate = 2)
#' rho=1/sd(dat)
#'
#'
#' #b=n
#' #a=sum(abs(dat) )
#' #rho=optimize(function(r){exp(-b)*(b/a)^b - (r^b)*exp(-a*r)}, c(0,20),  tol=10^(-50)  )$min
#'
#' out= isdexp(dat, B=1000000, alp , rho)
#' #prior samples
#' thetprior=unlist(out[2])
#' summary(thetprior)
#'
#' #posterior samples
#' thet=unlist(out[1])
#'
#' #95% Credible intervals
#' quantile (thet, c(0.025,0.975) )
#' summary(thet)
#'
#' #The accepted sample size:
#' unlist(out[3])
#'
#' #The acceptance probability:
#' unlist(out[4])
#'
#' #Plotting with normalization to have a maximum of 1
#' #for comparing prior and posterior
#' out2=density(thet)
#' ymaxpost=max(out2$y)
#' out3=density(thetprior)
#' ymaxprior=max(out3$y)
#' plot(out2$x,out2$y/ymaxpost, xlim=c(0,5), col="blue", type="l",
#'  xlab="theta", ylab="density", lwd=2, frame.plot=FALSE)
#' lines(out3$x,out3$y/ymaxprior,lwd=2,col="red")
#' #points(dat,rep(0,length(dat)), pch='*')
#'
#'
#' #Generating 1000 random numbers from the Inverse Stable (alpha=0.4,rho=5) prior
#' U1 = runif(1000)
#' U2 = runif(1000)
#' alp=0.4
#' rho=5
#' stab = ( ( sin(alp*pi*U1)*sin((1-alp)*pi*U1)^(1/alp-1) )
#' / (  ( sin(pi*U1)^(1/alp) )*abs(log(U2))^(1/alp-1))  )
#' #Inverse stable random numbers are below:
#' #rho*stab^(-alp)
#'
#'
#' @import stats nimble
#'
#'
#' @export
isdexp = function(x, B, alp, rho){
  B=B
  x=x
  alp=alp
  rho=rho
  U1 = runif(B)
  U2 = runif(B)
  rstab= ( ( sin(alp*pi*U1)*sin((1-alp)*pi*U1)^(1/alp-1) )/ (  ( sin(pi*U1)^(1/alp) )*abs(log(U2))^(1/alp-1))  )
  rstab=rstab[is.finite(rstab)]
  rstab=na.omit(rstab)
  n=length(x)
  b=n
  a=sum(abs(x) )
  v= rho*rstab^(-alp)
  n2=length(rstab)
  u=runif(n2)
  test=    -a*v+b+ b*log(a*v/b)
  test2= v*(  log(u) <test )
  thet=na.omit(test2[test2!=0])
  return(list(thet, v[1:length(thet)],length(thet),length(thet)/n2))
}
