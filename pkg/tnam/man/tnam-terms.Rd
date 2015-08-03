\name{tnam-terms}
\alias{tnam-terms}
\alias{tnam.terms}
\alias{terms-tnam}
\alias{terms.tnam}
\alias{netlag}
\alias{weightlag}
\alias{structsim}
\alias{covariate}
\alias{centrality}
\alias{degreedummy}
\alias{clustering}
\alias{cliquelag}
\alias{interact}
\alias{attribsim}
\docType{methods}
\title{Terms used in (Temporal) Network Autocorrelation Models (tnam)}
\description{
  The function \code{\link{tnam}} is used to fit (temporal) network 
  autocorrelation models.
  
  The function \code{\link{tnamdata}} can be used alternatively to 
  create a data frame containing all the data ready for estimation. 
  This may be useful when a non-standard model should be estimated, 
  like a tobit model or a model with zero inflation, for example.
  
  Both functions accept a formula containing several model terms. 
  The model terms are themselves functions which can be called 
  separately. For example, one model term is called \code{netlag}. 
  This model term can be part of the formula handed over to the 
  \code{tnam} function, or \code{netlag} can be called directly 
  in order to create a single variable.
  
  This help page describes the different model terms available in 
  (temporal) network autocorrelation models. See the 
  \code{\link{tnam}} help page for details on the model.
}

\section{Model terms for \code{\link{tnam}}}{
  \describe{
    
    \item{attribsim}{\emph{Spatial lag based on attribute 
    similarity}
    The \code{attribsim} model term computes a similarity matrix 
    based on the \code{attribute} argument and uses this similarity 
    matrix to construct a spatial lag by multiplying the similarity 
    matrix and the outcome vector \code{y}. The intuition behind 
    this model term is that node i's behavior may be influenced 
    by node j's behavior if nodes i and j are similar on another 
    dimension. For example, if i and j both smoke while k does not 
    smoke, j's alcohol consumption may affect i's alcohol 
    consumption to a larger extent than node k's alcohol 
    consumption. In this example, the \code{y} outcome variable 
    is alcohol consumption and the \code{attribute} argument is 
    smoking. If \code{match = FALSE}, the absolute similarity 
    between i and j is computed by subtracting j's attribute value 
    from i's attribute value and taking the absolute value to 
    construct the similarity matrix. If \code{match = TRUE}, the 
    function computes a matrix containing values of 1 if i and j 
    have the same attribute value and 0 otherwise. A scenario where 
    the \code{attribsim} model term makes sense is degree 
    assortativity: if i and j have the same degree centrality, they 
    may be inclined to learn from each other's behavior, even in 
    the absence of a direct connection between them.}
    
    \item{centrality}{\emph{Node centrality}
    The \code{centrality} model term computes a centrality index 
    for the nodes in a network or matrix. This can capture important 
    structural effects because being central often implies certain 
    constraints or opportunities more peripheral nodes do not have. 
    For example, central nodes in a network of employees might be 
    able to perform better.}
    
    \item{cliquelag}{\emph{Spatial lag of k-clique co-members}
    The \code{cliquelag} model term computes a clique co-membership 
    matrix and multiplies this matrix with the outcome variable. The 
    intuition behind this is that in some settings individuals may 
    be influenced to a particularly strong extent by peers in the 
    same cliques. A clique is defined as a maximal connected 
    subgraph of size k. For example, a deviant behavior of a person 
    may be conditioned by the deviant behavior of the person's 
    friends -- but only if these friends are tied to each other as 
    well so that a clique among these persons exists. A minimal and 
    a maximal \code{k} may be defined, where \code{k} is the size 
    of the cliques. In the clique co-membership matrix, all cliques 
    with \code{k.min <= k <= k.max} are included.}
    
    \item{clustering}{\emph{Local clustering coefficient}, or 
    \emph{transitivity}
    The \code{clustering} model term computes the local clustering 
    coefficient, which is also known as transitivity. This index 
    has high values if the direct neighborhood of a node is densely 
    interconnected. For example, if one's friends are friends with 
    each other, this may have repercussions on ego's behavior.}
    
    \item{covariate}{\emph{Exogenous nodal covariate}
    The \code{covariate} model term adds an exogenous nodal 
    covariate to the model. For example, when performance of 
    employees is modeled, a covariate could be seniority of these 
    employees. It is possible to add lagged covariates to model 
    the effect of past nodal attributes on current behavior. 
    Similarly, this model term can be used to add autoregressive 
    terms, that is, the effect of previous behavior on current 
    behavior.}
    
    \item{degreedummy}{\emph{Dummy variable for degree centrality 
    values}
    The \code{degreedummy} model term controls for specific degree 
    centralities or ranges of degree centrality. For example, 
    do nodes with a degree of 0 (isolates) show different behavior 
    than nodes who are connected? Or do nodes with a degree 
    centrality larger than three exert different behavior?}
    
    \item{interact}{\emph{Interactions between other model terms}
    The \code{interact} model term adds an interaction effect 
    between two other model terms by multiplying the result 
    vectors of these two model terms. When using interaction 
    terms, centering the result is recommended. Note that only 
    the interaction term is created; the main effects must be 
    introduced to the model using the other model terms.}
    
    \item{netlag}{\emph{Spatial network lag}
    The \code{netlag} model term captures the autocorrelation 
    inherent in networks. For example, when political actors 
    are members of a policy network, their success of achieving 
    policy outcomes is not independent from each other. Most 
    likely, being connected to policy winners increases the 
    success rate. In many settings, indirect effects may be 
    important as well: how does the behavior of my friends' 
    friends affect my own behavior? In some contexts, 
    spatio-temporal lags are useful: how does the past behavior 
    of my friends affect my current behavior? The \code{netlag} 
    model term is designed for binary networks because things 
    like indirect effects, restriction to reciprocal dyads, 
    decay of indirect relations etc. is possible. For weighted 
    networks, the \code{weightlag} term is recommended. If no 
    other arguments are specified and loops are absent and a 
    binary matrix is used, both model terms produce the same 
    results.}
    
    \item{structsim}{\emph{Structural similarity}
    The \code{structsim} model term computes the structural 
    similarity with other nodes in the network and multiplies 
    this similarity matrix with the outcome variable. The 
    intuition is that behavior is sometimes affected by 
    comparison with structurally similar nodes. For example, 
    a worker may be impressed by the performance of other 
    workers who are embedded in the same team or who report 
    to the same bosses. As with the other model terms, temporal 
    lags are possible.}
    
    \item{weightlag}{\emph{Weighted spatial lag}
    The \code{weightlag} model term captures spatial 
    autocorrelation in weighted networks. For example, the 
    GDP per capita of a country may be affected by the 
    GDP of proximate other countries or by the GDP of trade 
    partners. In these cases, indirect contacts etc. do not 
    make any sense, therefore the distinction between the 
    \code{weightlag} and the \code{netlag} model term. The 
    weight matrix is multiplied by the outcome variable, 
    possibly after row or column normalization.}
  }
}
\usage{
attribsim(y, attribute, match = FALSE, lag = 0, 
    normalization = c("no", "row", "column"), center = FALSE, 
    coefname = NULL)

centrality(networks, type = c("indegree", "outdegree", "freeman", 
    "betweenness", "flow", "closeness", "eigenvector", 
    "information", "load", "bonpow"), directed = TRUE, lag = 0, 
    rescale = FALSE, center = FALSE, coefname = NULL, ...)

cliquelag(y, networks, k.min = 2, k.max = Inf, directed = TRUE, 
    lag = 0, normalization = c("no", "row", "column"), 
    center = FALSE, coefname = NULL)

clustering(networks, directed = TRUE, lag = 0, center = FALSE, 
    coefname = NULL, ...)

covariate(y, lag = 0, exponent = 1, center = FALSE, 
    coefname = NULL)

degreedummy(networks, deg = 0, type = c("indegree", "outdegree", 
    "freeman"), reverse = FALSE, directed = TRUE, lag = 0, 
    center = FALSE, coefname = NULL, ...)

interact(x, y, lag = 0, center = FALSE, coefname = NULL)

netlag(y, networks, lag = 0, pathdist = 1, decay = pathdist^-1, 
    normalization = c("no", "row", "column", "complete"), 
    reciprocal = FALSE, center = FALSE, coefname = NULL, ...)

structsim(y, networks, lag = 0, method = c("euclidean", 
    "minkowski", "jaccard", "binary", "hamming"), center = FALSE, 
    coefname = NULL, ...)

weightlag(y, networks, lag = 0, normalization = c("no", "row", 
    "column"), center = FALSE, coefname = NULL)

}
\arguments{
  \item{attribute}{A vector, list of vectors or data frame with the same dimensions as \code{y}. Based on this attribute, the similarity between nodes i and j will be calculated, and the resulting similarity matrix is used to weight the \code{y} variable. }
  \item{center}{Should the model term be centered? That is, should the mean of the variable be subtracted from the actual value at each time step?}
  \item{coefname}{An additional name that is used as part of the coefficient label for easier identification in the summary output of the model.}
  \item{decay}{For each value in \code{pathdist}, the \code{decay} argument specifies the relative importance. By default, a geometric decay is used, that is, the behavior of nodes at path distance 2 is counted only half as much as the behavior of adjacent nodes. Alternatively, if both are equally important, it is possible to write \code{pathdist = c(1, 2)} and \code{decay = c(1, 1)}.}
  \item{deg}{The degree (e.g., \code{deg = 2}) or degree range (e.g., \code{deg = 1:3}).}
  \item{directed}{Is the input matrix or network a directed network?}
  \item{exponent}{The exponent of a covariate. For example, \code{exponent = 2} creates a squared variable. This may be helpful for modeling non-linear effects or for modeling a quadratic behavior shape.}
  \item{k.max}{Maximal clique size.}
  \item{k.min}{Minimal clique size.}
  \item{lag}{The temporal lag. The default value 0 means there is no lag. A value of 1 would specify a single-period lag, that is, current behavior is modeled conditional on previous influence. A value of 2 would specify a two-period lag, that is, current behavior is modeled conditional on pre-previous influence, etc.}
  \item{match}{If \code{match = FALSE}, a similarity matrix is computed by subtracting node j's attribute value from node i's attribute value, standardizing the resulting distance between 0 and 1, and converting it into a similarity by subtracting it from 1. This similarity matrix is used as a weight matrix to compute a spatial lag. If \code{match = TRUE} is specified, the weight matrix contains values of 1 whenever node i and j have the same attribute value and 0 otherwise.}
  \item{method}{The distance function used for computing structural similarity. Possible values are \code{"euclidean"}, \code{"minkowski"}, \code{"jaccard"}, \code{"binary"}, and \code{"hamming"}.}
  \item{networks}{The network(s) for computing the peer influence, also known as the weight matrix. This can be a matrix or a network object (for a single time step) or a list of matrices or network objects (for multiple time steps).}
  \item{normalization}{Possible values: \code{"no"} for switching off normalization, \code{"row"} for row normalization of the weight matrix, \code{"column"} for column normalization of the weight matrix, and \code{"complete"} for complete normalization.
  
  If \code{"no"} is selected, this corresponds to the total similarity or the sum of all influences of tied alters.
  
  If \code{"row"} is selected, this corresponds to the average alter effect. If \code{i} is the row node and \code{j} is the column node, row normalization computes the peer influence of \code{j} on \code{i} as a fraction of \code{i}'s overall number of outgoing ties (i.e., \code{i}'s row sum or outdegree centrality). The theoretical intuition is that other nodes' influences on \code{i} are not cumulative; \code{i} rather perceives the average influence of his or her peers. Note that row normalization does not necessarily entail that the values are standardized between 0 and 1.
  
  If \code{"column"} is selected, this corresponds to the peer influence of node \code{j} on node \code{i} as a fraction of the number of incoming ties \code{j} has (i.e., \code{j}'s column sum or indegree centrality). This captures the theoretical effect that \code{j} may distribute his or her influence among many nodes, in which case \code{j}'s influence on \code{i} is relatively weak. Thus the "exerted influence of \code{j} on \code{i} decreases with the number of actors \code{j} influences" (Leenders 2002). Note that column normalization does not necessarily entail that the values are standardized between 0 and 1.
  
  If \code{"complete"} is selected, this captures the peer influence of node \code{j} on node \code{i} as a fraction of all nodes' cumulative outcome values (including non-tied dyads; except \code{i}'s own outcome value). In other words, complete normalization corresponds to the actual exposure of \code{i} to the influence of his or her tied alters \code{j} over the the exposure \code{i} could receive if \code{i} were tied to all other nodes in the network. If a decay \code{<= 1} is used and if only direct friends are considered, this effectively standardizes the influence scores between \code{0} and \code{1}.}
  \item{pathdist}{An integer or a vector of integers. For example, if \code{pathdist = 1} is used, this computes the sum of the behavior of adjacent nodes. If \code{pathdist = 2} is specified, this computes the effect of indirect paths of length 2 ("friends of friends"). If \code{pathdist = 1:2} is set, both directly connected nodes' behavior and the behavior of nodes at a path distance of 2 from the focal node are counted. Arbitrary (sets of) path distances can be used. See also the \code{decay} argument.}
  \item{reciprocal}{If \code{reciprocal = TRUE} is specified, only the behavior of nodes to which a reciprocal relation exists is counted (that is, a link in both directions).}
  \item{rescale}{Should the centrality index be rescaled between 0 and 1?}
  \item{reverse}{Reverse the selection of degrees. For example, when \code{deg = 0} and \code{reverse = FALSE} are specified, resulting values of 1 indicate that a node has no connections, whereas the combination \code{deg = 0} and \code{reverse = TRUE} results in the value 1 representing nodes which have a degree of at least 1.}
  \item{type}{The type of centrality measure. Possible values are \code{"indegree"}, \code{"outdegree"}, \code{"freeman"}, \code{"betweenness"}, \code{"flow"}, \code{"closeness"}, \code{"eigenvector"}, \code{"information"}, \code{"load"}, and \code{"bonpow"}.}
  \item{x}{A variable that should be interacted with \code{y}. Either a vector or a list of vectors or another model term (this is the preferred way).}
  \item{y}{The outcome or behavior variable. Either a vector (for a single time step) or a list of vectors with named elements in each vector (for multiple time steps) or a data frame with row names where each column is one time step (for multiple time steps).}
  \item{...}{Additional arguments to be handed over to subroutines.}
}
\references{
Leenders, Roger Th. A. J. (2002):
Modeling Social Influence through Network Autocorrelation: 
Constructing the Weight Matrix. \emph{Social Networks} 24: 21--47.
\url{http://dx.doi.org/10.1016/S0378-8733(01)00049-1}.

Daraganova, Galina and Garry Robins (2013):
Autologistic Actor Attribute Models. In: Lusher, Dean, Johan 
Koskinen and Garry Robins, "Exponential Random Graph Models for 
Social Networks: Theory, Methods, and Applications", 
Cambridge University Press, chapter 9: 102--114.

Hays, Jude C., Aya Kachi and Robert J. Franzese Jr. (2010):
A Spatial Model Incorporating Dynamic, Endogenous Network 
Interdependence: A Political Science Application.
\emph{Statistical Methodology} 7: 406--428.
\url{http://dx.doi.org/10.1016/j.stamet.2009.11.005}
}
\seealso{
\link{tnam-package} \link{tnam} \link{tnamdata} \link[xergm.common]{knecht}
}
