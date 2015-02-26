\name{prediction}
\alias{prediction}

\title{ Prediction probabilities for Cox proportionnal hazard, Shared and Joint frailty models }

\description{

\bold{For Cox proportionnal hazard model}

A predictive probability of event between t and horizon time t+w, with w the window of prediction.

\deqn{P(t,t+w)=\frac{S_i(t)-S_i(t+w)}{S_i(t)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{\exp(\beta'Z_i)}}

\bold{For Gamma Shared Frailty model for clustered (not recurrent) events}

Two kinds of predictive probabilities can be calculated:

- a conditional predictive probability of event between t and horizon time t+w, i.e. given a specific group

\deqn{P^{cond}(t,t+w)=\frac{S_{ij}(t)-S_{ij}(t+w)}{S_{ij}(t)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{u_i\exp(\beta'Z_{ij})}}

- a marginal predictive probability of event between t and horizon time t+w, i.e. averaged over the population

\deqn{P^{marg}(t,t+w)=1-\left(\frac{1+\theta H_0(t)\exp(\beta'Z_{ij})}{1+\theta H_0(t+w)\exp(\beta'Z_{ij})}\right)^{1/ \theta}}

\bold{For Gaussian Shared Frailty model for clustered (not recurrent) events}

Two kinds of predictive probabilities can be calculated:

- a conditional predictive probability of event between t and horizon time t+w, i.e. given a specific group
    
\deqn{P^{cond}(t,t+w)=\frac{S_{ij}(t)-S_{ij}(t+w)}{S_{ij}(t)}=1-\left(\frac{S_0(t+w)}{S_0(t)}\right)^{\exp(\eta_i+\beta'Z_{ij})}}
     
- a marginal predictive probability of event between t and horizon time t+w, i.e. averaged over the population
       
\deqn{P^{marg}(t,t+w)=\frac{\int_{-\infty}^{+\infty}(S_{ij}(t)-S_{ij}(t+w))g(\eta)d\eta}{\int_{-\infty}^{+\infty}S_{ij}(t)g(\eta)d\eta}}

\bold{For Joint Frailty model}

You can predict risk of death knowing patients' characteristics i.e. predicting the probability of death in a specific time window given the history of the patient before the time of prediction t. The history of the patient is the information on covariates before time t, but also the number of recurrences and the time of occurences. Three types of marginal probabilities are computed:

- a prediction of death between t and t+w given that the patient had exactly J recurrences before t

\deqn{P^1(t,t+w)=P(D_i \le t+w|D_i>t,H_i^{J,1})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)](u_i)^JS_{i(J+1)}^R(t)g(u_i)du_i}{\int_0^\infty S_i^D(t)(u_i)^JS_{i(J+1)}^R(t)g(u_i)du_i}}

- a prediction of death between t and t+w given that the patient had at least J recurrences before t

\deqn{P^2(t,t+w)=P(D_i \le t+w|D_i>t,H_i^{J,2})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)](u_i)^JS_{iJ}^R(X_{iJ})g(u_i)du_i}{\int_0^\infty S_i^D(t)(u_i)^JS_{iJ}^R(X_{iJ})g(u_i)du_i}}

- a prediction of death between t and t+w considering the recurrence history only in the parameters estimation. It corresponds to the average probability of death between t and t+w for a patient with these given characteristics.

\deqn{P^3(t,t+w)=P(D_i \le t+w|D_i>t)=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+w)]g(u_i)du_i}{\int_0^\infty S_i^D(t)g(u_i)du_i}}

It is possible to compute all these predictions in two ways on a scale of times :
- either you want a cumulative probability of developping the event between t and t+w (with t fixed, but with a varying window of prediction w);
- either you want at a specific time the probability to develop the event in the next w (ie, for a varying prediction time t, but for a fixed window of prediction). See Details.

With Gaussian frailties (\eqn{\eta}), the same expressions are used but with \eqn{{u_i}^{J}} replaced by \eqn{\exp(J\eta_i)} and  \eqn{g(\eta)} corresponds to the gaussian distribution.
}

\usage{

prediction(fit, data, t, window, group, MC.sample=0)

}

\arguments{
  \item{fit}{A frailtyPenal or jointPenal object.}
  \item{data}{Dataframe for the prediction. See Details.}
  \item{t}{Time or vector of times for prediction.}
  \item{window}{Window or vector of windows for prediction.}
  \item{group}{Only for shared frailty model, the group on which you want to make the conditional prediction (its frailty value will be used). If you specify a group then a conditional prediction will be computed, otherwise it will be a marginal prediction.}
  \item{MC.sample}{Number of samples used to calculate confidence bands with a Monte-Carlo method (with a maximum of 1000 samples). If MC.sample=0 (default value), no confidence intervals are calculated.}
}

\value{

The following components are included in a 'predFrailty' object obtained by using prediction function for Cox proportionnal hazard and shared frailty model.

\item{npred}{Number of individual predictions}
\item{x.time}{A vector of prediction times of interest (used for plotting predictions): vector of prediction times t if fixed window. Otherwise vector of prediction times t+w}
\item{window}{Prediction window or vector of prediction windows}
\item{pred}{Predictions estimated for each profile}
\item{icproba}{Logical value. Were confidence intervals estimated ?}
\item{predLow}{Lower limit of Monte-Carlo confidence interval for each prediction}
\item{predHigh}{Upper limit of Monte-Carlo confidence interval for each prediction}
\item{type}{Type of prediction probability (marginal or conditional)}
\item{group}{For conditional probability, the group on which you make predictions}

The following components are included in a 'predJoint' object obtained by using prediction function for joint frailty model.

\item{npred}{Number of individual predictions}
\item{x.time}{A vector of prediction times of interest (used for plotting predictions): vector of prediction times t if fixed window. Otherwise vector of prediction times t+w}
\item{window}{Prediction window or vector of prediction windows}
\item{group}{Id of each patient}
\item{pred1}{Estimation of probability of type 1: exactly j recurrences}
\item{pred2}{Estimation of probability of type 2: at least j recurrences}
\item{pred3}{Estimation of probability of type 3}
\item{icproba}{Logical value. Were confidence intervals estimated ?}
\item{predlow1}{Lower limit of Monte-Carlo confidence interval for probability of type 1}
\item{predhigh1}{Upper limit of Monte-Carlo confidence interval for probability of type 1}
\item{predlow2}{Lower limit of Monte-Carlo confidence interval for probability of type 2}
\item{predhigh2}{Upper limit of Monte-Carlo confidence interval for probability of type 2}
\item{predlow3}{Lower limit of Monte-Carlo confidence interval for probability of type 3}
\item{predhigh3}{Upper limit of Monte-Carlo confidence interval for probability of type 3}

}

\details{
To compute predictions with a prediction time t fixed and a variable window:
\preformatted{prediction(fit, datapred, t=10, window=seq(1,10,by=1))}
Otherwise, you can have a variable prediction time and a fixed window.
\preformatted{prediction(fit, datapred, t=seq(10,20,by=1), window=5)}
Or fix both prediction time t and window.
\preformatted{prediction(fit, datapred, t=10, window=5)}

The dataframe building is an important step. It will contain profiles of patient on which you want to do predictions.
To make predictions on a Cox proportionnal hazard or a shared frailty model, only covariates need to be included. You have to distinguish between numerical and categorical variables (factors). If we fit a shared frailty model with two covariates sex (factor) and age (numeric), here is the associated dataframe for three profiles of prediction.

\preformatted{
datapred <- data.frame(sex=0,age=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex)<- c(1,2)
datapred[1,] <- c(1,40) # man, 40 years old
datapred[2,] <- c(2,45) # woman, 45 years old
datapred[3,] <- c(1,60) # man, 60 years old
}

To use the prediction function on a joint frailty model, the construction will be a little bit different. In this case, the prediction for the terminal event takes into account covariates but also history of relapse times for a patient. You have to create a dataframe with the relapse times, the indicator of event, the cluster variable and the covariates. Relapses occuring after the prediction time may be included but will be ignored for the prediction. A joint model with calendar-timescale need to be fitted with Surv(start,stop,event), relapse times correspond to the "stop" variable and indicators of event correspond to the "event" variable (if event=0, the relapse will not be taken into account). For patients without relapses, all the values of "event" variable should be set to 0. Finally, the same cluster variable name needs to be in the joint model and in the dataframe for predictions ("id" in the following example). For instance, we observe relapses of a disease and fit a joint model adjusted for two covariates sex (1:male 2:female) and chemo (treatment by chemotherapy 1:no 2:yes).
We describe 3 different profiles of prediction all treated by chemotherapy: 1) a man with four relapses at 100, 200, 300 and 400 days, 2) a man with only one relapse at 1000 days, 3) a woman without relapse.

\preformatted{
datapred <- data.frame(time=0,event=0,id=0,sex=0,chemo=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex) <- c(1,2)
datapred$chemo <- as.factor(datapred$chemo)
levels(datapred$chemo) <- c(1,2)
datapred[1,] <- c(100,1,1,1,2) # first relapse of the patient 1
datapred[2,] <- c(200,1,1,1,2) # second relapse of the patient 1
datapred[3,] <- c(300,1,1,1,2) # third relapse of the patient 1
datapred[4,] <- c(400,1,1,1,2) # fourth relapse of the patient 1
datapred[5,] <- c(1000,1,2,1,2) # one relapse at 1000 days for patient 2
datapred[6,] <- c(100,0,3,2,2) # patient 3 did not relapse
}

The data can also be the dataset used to fit the joint model. In this case, you will obtain as many prediction rows as patients.

}

\references{
A. Mauguen, B. Rachet, S. Mathoulin-Pelissier, G. MacGrogan, A. Laurent, V. Rondeau (2013). Dynamic prediction of risk of death using history of cancer recurrences in joint frailty models. \emph{Statistics in Medicine}, \bold{32(30)}, 5366-80.

V. Rondeau, A. Laurent, A. Mauguen, P. Joly, C. Helmer. Dynamic prediction models for clustered and interval-censored outcomes: investigating the intra-couple correlation in the risk of dementia. To appear.
}

\examples{

\dontrun{

#####################################################
#### prediction on a COX or SHARED frailty model ####
#####################################################

data(readmission)
#-- here is a generated cluster (31 clusters of 13 subjects)
readmission <- transform(readmission,group=id\%\%31+1)

#-- we compute predictions of death
#-- we extract last row of each subject for the time of death
readmission <- aggregate(readmission,by=list(readmission$id),
                         FUN=function(x){x[length(x)]})[,-1]

##-- predictions on a Cox proportionnal hazard model --##
cox <- frailtyPenal(Surv(t.stop,death)~sex+dukes,
n.knots=10,kappa=10000,data=readmission)

#-- construction of the dataframe for predictions
datapred <- data.frame(sex=0,dukes=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex)<- c(1,2)
datapred$dukes <- as.factor(datapred$dukes)
levels(datapred$dukes)<- c(1,2,3)
datapred[1,] <- c(1,2) # man, dukes 2
datapred[2,] <- c(2,3) # woman, dukes 3

#-- prediction of death for two patients between 100 and 100+w,
#-- with w in (50,100,...,1900)
pred.cox <- prediction(cox,datapred,100,seq(50,1900,50))
plot(pred.cox)

#-- prediction of death for two patients between t and t+400,
#-- with t in (100,150,...,1500)
pred.cox2 <- prediction(cox,datapred,seq(100,1500,50),400)
plot(pred.cox2)

##-- predictions on a shared frailty model for clustered data --##
sha <- frailtyPenal(Surv(t.stop,death)~cluster(group)+sex+dukes,
n.knots=10,kappa=10000,data=readmission)

#-- marginal prediction
pred.sha.marg <- prediction(sha,datapred,100,seq(50,1900,50))
plot(pred.sha.marg)

#-- conditional prediction, given a specific cluster (group=5)
pred.sha.cond <- prediction(sha,datapred,100,seq(50,1900,50),group=5)
plot(pred.sha.cond)

#####################################################
######## prediction on a JOINT frailty model ########
#####################################################

data(readmission)

##-- predictions of death on a joint model --##
joi <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)
+sex+dukes+terminal(death),formula.terminalEvent=~sex
+dukes,data=readmission,n.knots=10,kappa=c(100,100),recurrentAG=TRUE)

#-- construction of the dataframe for predictions
datapredj <- data.frame(t.stop=0,event=0,id=0,sex=0,dukes=0)
datapredj$sex <- as.factor(datapredj$sex)
levels(datapredj$sex) <- c(1,2)
datapredj$dukes <- as.factor(datapredj$dukes)
levels(datapredj$dukes) <- c(1,2,3)
datapredj[1,] <- c(100,1,1,1,2)
datapredj[2,] <- c(200,1,1,1,2)
datapredj[3,] <- c(300,1,1,1,2)
datapredj[4,] <- c(400,1,1,1,2)
datapredj[5,] <- c(380,1,2,1,2)

#-- prediction of death between 100 and 100+500 given relapses 
pred.joint0 <- prediction(joi,datapredj,100,500)
print(pred.joint0)

#-- prediction of death between 100 and 100+w given relapses (with confidence intervals)
pred.joint <- prediction(joi,datapredj,100,seq(50,1500,50),MC.sample=100)
plot(pred.joint,conf.bands=TRUE)
# each y-value of the plot corresponds to the prediction between [100,x]

#-- prediction of death between t and t+500 given relapses 
pred.joint2 <- prediction(joi,datapredj,seq(100,1000,50),500)
plot(pred.joint2)
# each y-value of the plot corresponds to the prediction between [x,x+500], or in the next 500

}

}

\keyword{misc}