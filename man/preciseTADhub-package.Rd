\name{preciseTADhub-package}

\alias{preciseTADhub-package}
\alias{preciseTADhub}

\title{Pre-trained models obtained using preciseTAD as \code{list} objects.}

\description{
    preciseTADhub is package that give users access to pre-trained random forest
    models that can be leveraged to predict TAD and/or chromatin loop boundaries
    using the preciseTAD R package. These data have been parsed into \code{list}
    objects and RDS files and are available in \code{ExperimentHub}.
}

\details{
    See the vignette for examples of using these data in predicting precise
    boundary location at base-level resolution.

    browseVignettes("preciseTADhub")

    Details of how these data were created are in the scripts/ directory
    of the source package.
}

\examples{
\dontrun{
library(ExperimentHub)
eh <- ExperimentHub()
myfiles <- query(eh, "preciseTADhub")
CHR1_GM12878_5kb_Arrowhead <- myfiles[[1]]
}}

\keyword{utilities}
