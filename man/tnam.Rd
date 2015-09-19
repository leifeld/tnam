\name{tnam}
\alias{tnam}
\alias{tnamdata}
\title{Fit (temporal) network autocorrelation models}
\description{Fit (temporal) network autocorrelation models.}
\usage{
tnam(formula, family = gaussian, re.node = FALSE, 
    re.time = FALSE, time.linear = FALSE, time.quadratic = FALSE, 
    center.y = FALSE, na.action = na.omit, ...)

tnamdata(formula, center.y = FALSE)
}
\details{
The \code{tnam} function serves to estimate temporal or cross-sectional network autocorrelation models. Model terms such as spatial lags, temporal lags, spatio-temporal lags, centrality etc. can be specified in the \code{formula} argument. Details on the model terms can be found on the \code{\link{tnam-terms}} help page.

The \code{tnamdata} function accepts a formula (like in the \code{tnam} function) and returns a data frame with the response variable and the covariates for estimation with any estimation function. \code{tnam} first calls \code{tnamdata} internally and then hands over the resulting data structure to a \code{glm}, \code{lmer}, or \code{nlmer} call. If models such as tobit, multinomial or multilevel models should be estimated, one can leave out the estimation step and feed the results of \code{tnamdata} manually into any type of model.
}
\arguments{
\item{formula}{ A formula where the left-hand side specifies either a vector containing the outcome variable (for a cross-sectional model) or a list of such vectors (for modeling the outcome at multiple time steps) or a data frame with one time step per column (also for longitudinal models of behavior). The right-hand side of the formula consists of tnam-specific model terms like \code{netlag}, \code{structsim} and other terms which are described on the help page of \link{tnam-terms}. }
\item{family}{ The link function for fitting the generalized linear model or the mixed effects model, for example \code{gaussian} or \code{binomial}. The options are the same as in the \code{\link[stats]{glm}} and \code{\link[lme4]{glmer}} functions. For details on the \code{family} argument, see the \code{\link[stats]{family}} help page. }
\item{re.node}{ If multiple time steps are present: should a random effect for the nodes be added to the model? This results in the estimation of a mixed effects model. }
\item{re.time}{ If multiple time steps are present: should a random effect for the time steps be added to the model? This results in the estimation of a mixed effects model. }
\item{time.linear}{ If multiple time steps are present: should a linear effect for time be added to the model? This can be estimated in the standard GLM framework. }
\item{time.quadratic}{ If multiple time steps are present: should a squared effect for time be added to the model? This can be estimated in the standard GLM framework. }
\item{center.y}{ Center the dependent variable by subtracting the mean from the actual value within each time step? }
\item{na.action}{ How should missing values be treated? By default, they are omitted. See the \code{\link{na.omit}} help page for details. }
\item{...}{ Further arguments that should be passed to the \code{\link{glm}}, \code{\link{lmer}}, or \code{\link{glmer}} function, which is used under the hood for estimating the model. }
}
\examples{
# The following example models delinquency among adolescents at 
# multiple time steps as a function of (1) their nodal attributes 
# like sex or religion, (2) their peers' delinquency levels, (3) 
# their own and their peers' past delinquency behavior, and (4) 
# their structural position in the network. See ?knecht for 
# details on the dataset. Before estimating the model, all data 
# should be labeled with the names of the nodes such that tnam 
# is able to merge the information on multiple nodes across time 
# points.

library("tnam")
data("knecht")

# prepare the dependent variable y
delinquency <- as.data.frame(delinquency)
rownames(delinquency) <- letters

# replace structural zeros (denoted as 10) and add row labels
friendship[[3]][friendship[[3]] == 10] <- NA
friendship[[4]][friendship[[4]] == 10] <- NA
for (i in 1:length(friendship)) {
  rownames(friendship[[i]]) <- letters
}

# prepare the covariates sex and religion
sex <- demographics$sex
names(sex) <- letters
sex <- list(t1 = sex, t2 = sex, t3 = sex, t4 = sex)
religion <- demographics$religion
names(religion) <- letters
religion <- list(t1 = religion, t2 = religion, t3 = religion, 
    t4 = religion)

# Estimate the model. The first term is the sex of the respondent, 
# the second term is the religion of the respondent, the third 
# term is the previous delinquency behavior of the respondent, 
# the fourth term is the delinquency behavior of direct friends, 
# the fifth term is the delinquency behavior of indirect friends 
# at a path distance of 2, the sixth effect is the past delinquency 
# of direct friends, the seventh term indicates whether the 
# respondent has any contacts at all, and the last term captures 
# the effect of the betweenness centrality of the respondent on 
# his or her behavior. Apparently, previous behavior, being an 
# isolate, and religion seem to have an effect on delinquency in 
# this dataset. There is also a slight positive trend over time, 
# and direct friends exert a minor effect (not significant).
# Note that a linear model may not be the best specification for 
# modeling the ordered categorical delinquency variable, but it 
# suffice here for illustration purposes.

model1 <- tnam(
    delinquency ~ 
    covariate(sex, coefname = "sex") + 
    covariate(religion, coefname = "religion") + 
    covariate(delinquency, lag = 1, exponent = 1) + 
    netlag(delinquency, friendship) + 
    netlag(delinquency, friendship, pathdist = 2, decay = 1) + 
    netlag(delinquency, friendship, lag = 1) + 
    degreedummy(friendship, deg = 0, reverse = TRUE) + 
    centrality(friendship, type = "betweenness"), 
    re.node = TRUE, time.linear = TRUE
)
summary(model1)

# for nice table output, use the texreg package
library("texreg")
screenreg(model1)
}
\seealso{
\link{tnam-package} \link{tnam-terms} \link[xergm.common]{knecht}
}
