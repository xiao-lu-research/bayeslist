#' The 'bayeslist' package.
#'
#' @description Estimate list experiment data applying full Bayesian approaches.
#'
#' @aliases bayeslist-package
"_PACKAGE"
#' @useDynLib bayeslist, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. https://mc-stan.org
#'
NULL

#' The Sri Lanka List Experiment on Wartime Sexual Violence
#'
#' This dataset, which includes male respondents from Tamil, is a subset of the list experiment administered in Sri Lanka on wartime sexual violence.The main question reads as follows:
#' Now we would like to ask you some more questions about what happened during the war.
#' Please refer to the following list and tell me how many of these experiences happened to you during the war.
#' Please don’t tell me which specific statements you believe to be true, only how many:
#' (1) I won money in a lottery or competition;
#' (2) I was involved in an accident;
#' (3) I received help from a stranger;
#' (4) (I was personally sexually assaulted.)
#' The forth item in bracket is the sensitive item.
#' In addition to the above list, there are also two direct questions asking about sexual abuse.
#'
#'
#'
#' @name srilanka
#' @docType data
#' @format A data frame containing the following 9 variables for 247 observations.
#'
#' \tabular{lll}{
#' sexaussault \tab integer \tab Reported item count for list experiment question. \cr
#' sexaussault_d \tab integer \tab First direct item. \cr
#' sexaussault_d2 \tab integer \tab Second direct item. \cr
#' treatment \tab integer \tab Indicator for list experiment treatment group. \cr
#' age \tab numeric \tab Age. \cr
#' edu \tab integer \tab Education. \cr
#' eastern \tab integer \tab Whether the respondent comes from eastern Tamil. \cr
#' assist.army \tab integer \tab Whether the respondent has assisted rebel groups. \cr
#' displace \tab integer \tab Displacement. \cr
#' }
#'
#' @references Traunmüller, R., Kijewski, S., & Freitag, M. (2019). The silent victims of sexual violence during war: Evidence from a list experiment in Sri Lanka. Journal of conflict resolution, 63(9), 2015-2042. \doi{10.1177/0022002719828053}
#' @keywords datasets
NULL



#' The 2017 London List Experiment
#'
#' This dataset is the 2017 London list experiment on voter turnout fielded via online YouGov survey of a sample of 3189 Greater Londoners.The main question reads as follows:
#' The next question deals with the recent general election on 8th June.
#' Here is a list of four (five) things that some people did and some people did not do during the election campaign or on Election Day.
#' Please say how many of these things you did. Here are the four (five) things:
#' (1) Discussed the election with family and friends;
#' (2) (Voted in the election);
#' (3) Criticised a politician on social media;
#' (4) Avoided watching the leaders debates;
#' (5) Put up a poster for a political party in my window or garden.
#' How many of these things did you do?
#' The second item in bracket is the sensitive item.
#' In addition to the above list, there is a direct question asking about turnout:
#' Talking with people about the recent general election on 8th June, we have found that a lot of people didn't manage to vote. How about you, did you manage to vote in the general election?
#'
#'
#'
#' @name london
#' @docType data
#' @format A data frame containing the following 18 variables for 3189 observations.
#'
#' \tabular{lll}{
#' ID \tab integer \tab Respondent ID number. \cr
#' age \tab integer \tab Respondent age in years. \cr
#' agegrp \tab factor \tab Respondent age group. \cr
#' gender \tab factor \tab YouGov panel measure of gender. \cr
#' social_grade \tab factor \tab YouGov panel measure of respondent social grade. \cr
#' qual \tab factor \tab Measure of highest educational qualification from YouGov panel. \cr
#' validationfactor \tab factor \tab Detailed measure of turnout validation outcome for respondent. \cr
#' validturnout \tab integer \tab Summary measure of true respondent turnout. \cr
#' direct \tab integer \tab Response to direct turnout question asked of list experiment control group. \cr
#' baselineTurnout \tab integer \tab Response to baseline direct turnout question after the election. \cr
#' listTreat \tab integer \tab Indicator for list experiment treatment group. \cr
#' listCount \tab integer \tab Reported item count for list experiment question. \cr
#' qtime \tab numeric \tab Time taken to answer list experiment question, in seconds. \cr
#' recallfirst \tab character \tab Respondent recall of first item from list question. Open text response. \cr
#' recalllast \tab character \tab Respondent recall of last item from list question. Open text response. \cr
#' recallfirst.hand.correct \tab factor \tab Did respondent correctly recall first list experiment item? \cr
#' recalllast.hand.correct \tab factor \tab Did respondent correctly recall last list experiment item? \cr
#' comfort \tab numeric \tab How comfortable do you feel revealing whether you did/did not vote in last election? \cr
#' }
#'
#' @references Kuhn, P. M., & Vivyan, N. (2021). The misreporting trade-off between list experiments and direct questions in practice: Partition validation evidence from two countries. Political Analysis, 1-22.
#' @source The full data set is available at \doi{10.7910/DVN/W90Q7B})
#' @keywords datasets
NULL


