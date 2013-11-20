
\name{prediction}
\alias{prediction}

\title{ Prediction probabilities for Cox proportionnal hazard, Shared and Joint frailty model }

\description{

\bold{For Cox proportionnal hazard model}

A predictive probability of event between t and horizon time t+h.

\deqn{P(t,t+h)=\frac{S_i(t)-S_i(t+h)}{S_i(t)}=1-\left(\frac{S_0(t+h)}{S_0(t)}\right)^{\exp(\beta'Z_i)}}

\bold{For Shared Frailty model for clustered (not recurrent) events}

Two kinds of predictive probabilities can be calculated:

- a conditional predictive probability of event between t and horizon time t+h, i.e. given a specific group

\deqn{P^{cond}(t,t+h)=\frac{S_{ij}(t)-S_{ij}(t+h)}{S_{ij}(t)}=1-\left(\frac{S_0(t+h)}{S_0(t)}\right)^{\omega_i\exp(\beta'Z_{ij})}}

- a marginal predictive probability of event between t and horizon time t+h, i.e. averaged over the population

\deqn{P^{marg}(t,t+h)=1-\left(\frac{1+\theta H_0(t)\exp(\beta'Z_{ij}(t))}{1+\theta H_0(t+h)\exp(\beta'Z_{ij}(t))}\right)^{1/ \theta}}

\bold{For Joint Frailty model}

You can predict risk of death knowing patients' prognosis i.e. predicting the probability of death in a specific time window given the history of the patient before the time of prediction t. The history of the patient is the information on covariates before time t, but also the number of recurrences and the time of occurences. Three types of probability are computed:

- a probability of death between t and t+h given that the patient has exactly J recurrences before t

\deqn{P^1(t,t+h)=P(D_i \le t+h|D_i>t,H_i^{J,1})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+h)](u_i)^JS_{i(J+1)}^R(t)g(u_i)du_i}{\int_0^\infty S_i^D(t)(u_i)^JS_{i(J+1)}^R(t)g(u_i)du_i}}

- a probability of death between t and t+h given that the patient has at least J recurrences before t

\deqn{P^2(t,t+h)=P(D_i \le t+h|D_i>t,H_i^{J,2})=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+h)](u_i)^JS_{iJ}^R(X_{iJ})g(u_i)du_i}{\int_0^\infty S_i^D(t)(u_i)^JS_{iJ}^R(X_{iJ})g(u_i)du_i}}

- a probability of death between t and t+h considering the recurrence history only in the parameters estimation. It corresponds to the average probability of death between t and t+h for a patient with these given characteristics.

\deqn{P^3(t,t+h)=P(D_i \le t+h|D_i>t)=\frac{\int_0^\infty[S_i^D(t)-S_i^D(t+h)]g(u_i)du_i}{\int_0^\infty S_i^D(t)g(u_i)du_i}}

}

\usage{

prediction(fit, data, predTime, horizon, predMax, type="marginal",
           group, MC.sample=0)
             
}

\arguments{
  \item{fit}{A frailtyPenal or jointPenal object}
  \item{data}{Dataframe for the prediction. See Details.}
  \item{predTime}{First time of predictions}
  \item{horizon}{Length of the split}
  \item{predMax}{Last time of predictions}
  \item{type}{Only for shared frailty model, type of the prediction "marginal" or "conditional". Need to specify argument 'group' if type="conditional". Marginal by default}
  \item{group}{Only for shared frailty model, the group on which you want to make the conditional prediction (its frailty value will be used)}
  \item{MC.sample}{Number of samples used to calculate confidence bands with a Monte-Carlo method. If MC.sample=0 (default value), no confidence bands calculated.}
}

\value{

The following components are included in a 'predFrailty' object obtained by using prediction function for Cox proportionnal hazard and shared frailty model.

\item{time}{Vector of prediction times}
\item{pred}{Predictions estimated for each profile}
\item{icproba}{Logical value. Were confidence bands estimated ?}
\item{predLow}{Lower limit of confidence interval for each prediction}
\item{predHigh}{Upper limit of confidence interval for each prediction}
\item{type}{Type of prediction probability (marginal or conditional)}
\item{horizon}{Prediction horizon}
\item{group}{For conditional probability, the group on which you make predictions}

The following components are included in a 'predJoint' object obtained by using prediction function for joint frailty model.

\item{npred}{Number of individual predictions}
\item{time}{Vector of prediction times}
\item{group}{Id of each patient}
\item{proba1}{Estimation of probability of type 1: exactly j recurrences}
\item{proba2}{Estimation of probability of type 2: at least j recurrences}
\item{proba3}{Estimation of probability of type 3}
\item{icproba}{Logical value. Were confidence bands estimated ?}
\item{probalow1}{Lower limit of Monte-Carlo confidence interval for probability of type 1}
\item{probahigh1}{Upper limit of Monte-Carlo confidence interval for probability of type 1}
\item{probalow2}{Lower limit of Monte-Carlo confidence interval for probability of type 2}
\item{probahigh2}{Upper limit of Monte-Carlo confidence interval for probability of type 2}
\item{probalow3}{Lower limit of Monte-Carlo confidence interval for probability of type 3}
\item{probahigh3}{Upper limit of Monte-Carlo confidence interval for probability of type 3}

}

\details{
The dataframe building is an important step. It will contain profiles of patient on which you want to do predictions.
To make predictions on a Cox proportionnal hazard or a shared frailty model, only covariates need to be included. You have to distinguish between numerical and categorical variables (factors). If we fit a shared frailty model with two covariates sex (factor 1:male 2:female) and age (numeric), here is the associated dataframe for the three profiles of prediction.

\preformatted{
datapred <- data.frame(sex=0,age=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex)<- c(1,2)
datapred[1,] <- c(1,40) # man 40 years old
datapred[2,] <- c(2,45) # woman 45 years old
datapred[3,] <- c(1,60) # man 60 years old
}

If you want to use the prediction function on a joint frailty model, the construction will be a little bit different. In this case, the prediction for the terminal event takes into account covariates but also history of relapse times for a patient. You have to create as much rows in the dataframe as there are relapses before the prediction time for a specific patient. If you did your joint model in gap-time using Surv(time,event), relapse times correspond to the "time" variable. In case of calendar-time with Surv(start,stop,event), relapse times correspond to the "stop" variable. For patients without relapses, "time" or "stop" variable should be equal to 0. Finally, the same cluster variable name needs to be in the joint model and in the dataframe for predictions. For instance, we observe relapses of a disease and fit a joint model adjusted for two covariates sex (1:male 2:female) and chemo (treatment by chemotherapy 1:no 2:yes). We describe 3 different profiles of prediction: 1) a man with four relapses at 100, 200, 300 and 400 days, 2) a man with only one relapse at 1000 days, 3) a woman without relapse, all treated by chemotherapy.

\preformatted{
datapred <- data.frame(time=0,id=0,sex=0,chemo=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex) <- c(1,2)
datapred$chemo <- as.factor(datapred$chemo)
levels(datapred$chemo) <- c(1,2)
datapred[1,] <- c(100,1,1,2) # first relapse of the patient 1
datapred[2,] <- c(200,1,1,2) # second relapse of the patient 1
datapred[3,] <- c(300,1,1,2) # third relapse of the patient 1
datapred[4,] <- c(400,1,1,2) # fourth relapse of the patient 1
datapred[5,] <- c(1000,2,1,2) # one relapse at 1000 days for patient 2
datapred[6,] <- c(0,3,2,2) # patient 3 did not relapse
}

}

\references{
A. Mauguen, B. Rachet, S. Mathoulin-Pelissier, G. MacGrogan, A. Laurent, V. Rondeau (2013). Dynamic prediction of risk of death using history of cancer recurrences in joint frailty models. Statistics in Medicine.

V. Rondeau, A. Laurent, A. Mauguen, P. Joly, C. Helmer. Simple or joint frailty models for clustered and interval-censored outcome with a dependant terminal event: investigating the intra-couple correlation in the risk of dementia and death. To appear.
}

\examples{

\dontrun{

data(readmission)

##-- predictions on a Cox proportionnal hazard model --##
cox <- frailtyPenal(Surv(time,event)~sex+dukes,
n.knots=10,kappa1=10000,data=readmission,Frailty=FALSE)

#-- construction of the dataframe for predictions
datapred <- data.frame(sex=0,dukes=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex)<- c(1,2)
datapred$dukes <- as.factor(datapred$dukes)
levels(datapred$dukes)<- c(1,2,3)
datapred[1,] <- c(1,2)
datapred[2,] <- c(2,3)

#-- prediction of the event for two patients, between time 100 and time 2000
pred.cox <- prediction(cox,datapred,100,50,2000)
plot(pred.cox)

##-- predictions on a shared frailty model --##
sha <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes,
n.knots=10,kappa1=10000,data=readmission,Frailty=TRUE)

#-- marginal prediction
pred.sha.marg <- prediction(sha,datapred,100,50,2000,type="marginal")
plot(pred.sha.marg)

#-- conditional prediction, given a specific cluster (group=5)
pred.sha.cond <- prediction(sha,datapred,100,50,2000,type="conditional",group=5)
plot(pred.sha.cond)

##-- predictions on a joint frailty model --##
joi <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)
+sex+dukes+terminal(death),formula.terminalEvent=~sex
+dukes,data=readmission,n.knots=10,kappa1=9.55e9,
kappa2=1.41e12,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE)

#-- creation of the dataframe for predictions
datapredj <- data.frame(t.stop=0,id=0,sex=0,dukes=0)
datapredj$sex <- as.factor(datapredj$sex)
levels(datapredj$sex) <- c(1,2)
datapredj$dukes <- as.factor(datapredj$dukes)
levels(datapredj$dukes) <- c(1,2,3)
datapredj[1,] <- c(100,1,1,2)
datapredj[2,] <- c(200,1,1,2)
datapredj[3,] <- c(300,1,1,2)
datapredj[4,] <- c(400,1,1,2)
datapredj[5,] <- c(380,2,1,2)

#-- only recurrences observed before the prediction
#-- time should be included in the 'datapred'
pred.joint <- prediction(joi,datapredj,500,50,2000)
plot(pred.joint)

}

}

\keyword{misc}
