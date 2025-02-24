#' U.S. real GDP percent change and GDP implicit price deflator percent change.
#'
#' A dataset containing a quarterly U.S. time series with two components:
#' the percentage change of real GDP and the percentage change of GDP implicit price deflator,
#' covering the period from 1959Q1 - 2019Q4.
#'
#' @format A numeric matrix of class \code{'ts'} with 244 rows and 2 columns with one time series in each column:
#' \describe{
#'   \item{First column (GDP):}{The quarterly percent change of real U.S. GDP, from 1959Q1 to 2019Q4, \url{https://fred.stlouisfed.org/series/GDPC1}.}
#'   \item{Second column (GDPDEF):}{The quarterly percent change of U.S. GDP implicit price deflator, from 1959Q1 to 2019Q4, \url{https://fred.stlouisfed.org/series/GDPDEF}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database
"gdpdef"


#' A quarterly U.S. data covering the period from 1954Q3 to 2021Q4 (270 observations) and consisting three variables:
#' cyclical component of the log of real GDP, the log-difference of GDP implicit price deflator, and an interest rate variable.
#' The interest rate variable is the effective federal funds rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016)
#' shadow rate, which is not constrained by the zero lower bound and also quantifies unconventional monetary policy measures.
#' The log-differences of the GDP deflator and producer price index are multiplied by hundred.
#'
#' The cyclical component of the log of real GDP was obtained by applying a one-sided Hodrick-Prescott (HP) filter with the
#' standard smoothing parameter lambda=1600. The one-sided filter was obtained from the two-sided HP filter by applying the
#' filter up to horizon t, taking the last observation, and repeating this procedure for the full sample t=1,...,T.
#' In order to allow the series to start from any phase of the cycle, we applied the one-sided filter to the full available
#' sample from 1947Q1 to 2021Q1 before extracting our sample period from it. We computed the two-sided HP filters with the R
#' package lpirfs (Adämmer, 2021)
#'
#' @format A numeric matrix of class \code{'ts'} with 270 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{First column (GDP):}{The cyclical component of the log of real GDP, \url{https://fred.stlouisfed.org/series/GDPC1}.}
#'   \item{Second column (GDPDEF):}{The log-difference of GDP implicit price deflator, \url{https://fred.stlouisfed.org/series/GDPDEF}.}
#'   \item{Third column (RATE):}{The Federal funds rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016) shadow rate,
#'    \url{https://fred.stlouisfed.org/series/FEDFUNDS}, \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
#' @references
#'  \itemize{
#'    \item Adämmer P. 2021. lprfs: Local Projections Impulse Response Functions. R package version: 0.2.0,
#'      \url{https://CRAN.R-project.org/package=lpirfs}.
#'    \item Wu J. and Xia F. 2016. Measuring the macroeconomic impact of monetary policy at the zero lower bound.
#'      \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
"usamone"


#' A monthly  U.S. data covering the period from 1961I to 2022III (735 observations) and consisting four variables.
#' First, The Actuaries Climate Index (ACI), which is a measure of the frequency of severe weather and the extend changes in sea levels.
#' Second, the monthly GDP growth rate constructed by the Federal Reserve Bank of Chicago from a collapsed dynamic factor analysis of
#' a panel of 500 monthly measures of real economic activity and quarterly real GDP growth. Third, the monthly growth rate of the
#' consumer price index (CPI). Third, an interest rate variable, which is the effective federal funds rate that is replaced by the
#' the Wu and Xia (2016) shadow rate during zero-lower-bound periods. The Wu and Xia (2016) shadow rate is not bounded by the zero
#' lower bound and also quantifies unconventional monetary policy measures, while it closely follows the federal funds rate when the
#' zero lower bound does not bind.
#'
#' @format A numeric matrix of class \code{'ts'} with 735 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{First column (GDP):}{The cyclical component of the log of real GDP, \url{https://fred.stlouisfed.org/series/GDPC1}.}
#'   \item{Second column (GDPDEF):}{The log-difference of GDP implicit price deflator, \url{https://fred.stlouisfed.org/series/GDPDEF}.}
#'   \item{Third column (RATE):}{The Federal funds rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016) shadow rate,
#'    \url{https://fred.stlouisfed.org/series/FEDFUNDS}, \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
#' @references
#'  \itemize{
#'    \item American Academy of Actuaries, Canadian Institute of Actuaries, Casualty Actuarial Society,
#'     and Society of Actuaries, 2023. Actuaries Climate Index. https://actuariesclimateindex.org.
#'   \item Federal Reserve Bank of Chicago, 2023. Monthly GDP Growth Rate Data. \url{https://www.chicagofed.org/publications/bbki/index}.
#'    \item Wu J. and Xia F. 2016. Measuring the macroeconomic impact of monetary policy at the zero lower bound.
#'      \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
"acidata"


#' A monthly  U.S. data covering the period from 1987:4 to 2024:12 (453 observations) and consisting six variables.
#' First, the climate policy uncertainty index (CPUI) (Gavridiilis, 2021), which is a news based measure of climate policy uncertainty.
#' Second, the economic policy uncertainty index (EPUI), which is a news based measure of economic policy uncertainty, including also
#' components quantifying the present value of future scheduled tax code expirations and disagreement among professional forecasters
#' over future goverment purchases and consumer prices.
#' Third, the log-difference of real indsitrial production index (IPI).
#' Fourth, the log-difference of the consumer price index (CPI).
#' Fifth, the log-difference of the producer price index (PPI).
#' Sixth, an interest rate variable, which is the effective federal funds rate that is replaced by the
#' the Wu and Xia (2016) shadow rate during zero-lower-bound periods. The Wu and Xia (2016) shadow rate is not bounded by the zero
#' lower bound and also quantifies unconventional monetary policy measures, while it closely follows the federal funds rate when the
#' zero lower bound does not bind. This is the dataset used in Virolainen (2025)
#'
#' @format A numeric matrix of class \code{'ts'} with 443 rows and 4 columns with one time series in each column:
#' \describe{
#'   \item{First column (CPUI):}{The climate policy uncertainty index, \url{https://www.policyuncertainty.com/climate_uncertainty.html}.}
#'   \item{Second column (EPUI):}{The economic policy uncertainty index, \url{https://www.policyuncertainty.com/us_monthly.html}.}
#'   \item{Third column (IPI):}{The log-difference of real indsitrial production index, \url{https://fred.stlouisfed.org/series/INDPRO}.}
#'   \item{Fourth column (CPI):}{The log-difference of the consumer price index, \url{https://fred.stlouisfed.org/series/CPIAUCSL}.}
#'   \item{Fifth column (PPI):}{The log-difference of the producer price index, \url{https://fred.stlouisfed.org/series/PPIACO}.}
#'   \item{Sixth column (RATE):}{The Federal funds rate from 1954Q3 to 2008Q2 and after that the Wu and Xia (2016) shadow rate,
#'    \url{https://fred.stlouisfed.org/series/FEDFUNDS}, \url{https://www.atlantafed.org/cqer/research/wu-xia-shadow-federal-funds-rate}.}
#' }
#'
#' @source The Federal Reserve Bank of St. Louis database and the Federal Reserve Bank of Atlanta's website
#' @references
#'  \itemize{
#'    \item K. Gavriilidis, 2021. Measuring climate policy uncertainty. \url{https://www.ssrn.com/abstract=3847388}.
#'    \item Federal Reserve Bank of Chicago. 2023. Monthly GDP Growth Rate Data. \url{https://www.chicagofed.org/publications/bbki/index}.
#'    \item Virolainen S. 2025. Identification by non-Gaussianity in structural threshold and smooth transition vector
#'      autoregressive models. Unpublished working paper, available as arXiv.2404.19707.
#'    \item Wu J. and Xia F. 2016. Measuring the macroeconomic impact of monetary policy at the zero lower bound.
#'      \emph{Journal of Money, Credit and Banking}, 48(2-3): 253-291.
#'  }
"usacpu"
