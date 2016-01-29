shinyUI(fluidPage(
  titlePanel(title="BFtestvar"),
  "Bayes Factors for TESTing VARiances. This", a("Shiny", href = "http://shiny.rstudio.com/"),
  "application computes the adjusted fractional Bayes factor presented in Anonymous (2015).",
  em("Bayesian evaluation of equality and inequality constrained hypotheses on variances."),
  "Manuscript submitted for publication. The R source code is available on", a("GitHub.", href =
  "https://github.com/fboeingmessing/BFtestvar"),
  br(),
  br(),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Mandatory input",
                 br(),
                 textInput("s2", label = "Sample variances", value = "1, 1.5, 2.25"),
                 textInput("n", label = "Sample sizes", value = "20, 20, 20"),
                 textInput("hypotheses", label = "Hypotheses",
                           value = "1<2<3; not (1<2<3 or 3<2<1); 1<2=3; 1<(2,3)")
        ),

        tabPanel("Optional input",
                 br(),
                 h5(strong("Logarithm of Bayes factors")),
                 checkboxInput("log.BF", "Show logarithm of Bayes factors", value = FALSE),
                 textInput("prior.probabilities", label = "Prior probabilities of the hypotheses",
                           value = NA),
                 textInput("b", label = "Fractions", value = NA),
                 numericInput("nsim",
                              label = "Number of draws from the posterior distribution",
                              value = 100000),
                 numericInput("seed", label = "Seed", value = NA)
        ),
        
        tabPanel("Help",
                 br(),
                 p(strong("Mandatory input")),
                 p(em("* Sample variances:"), "Separate sample variances by a comma (,). Make sure to
                   use a point (.) as the decimal mark."),
                 p(em("* Sample sizes:"), "Separate sample sizes by a comma (,)."),
                 p(em("* Hypotheses:"), "Use numbers 1,...,J to specify at least two hypotheses, where
                   J is the number of groups. The number 1 refers to the population variance of group 1
                   and so on. For an overview of the types of hypotheses we may specify see
                   Böing-Messing et al. (2015). Separate hypotheses by a semicolon (;)."),
                 br(),
                 p(strong("Optional input")),
                 p(em("* Logarithm of Bayes factors:"), "Determines whether the app shows Bayes factors
                   or log Bayes factors. By default the app shows Bayes factors."),
                 p(em("* Prior probabilities of the hypotheses:"), 'Are used to compute the posterior
                   probabilities of the hypotheses. Separate prior probabilities by a comma (e.g. "1/2,
                   1/2"). Fractions are allowed. Note that the prior probabilities need to sum to 1. If
                   no prior probabilities are specified (i.e. if the field is empty), then the app
                   assumes equal prior probabilities.'),
                 # p(HTML(paste("2/n", tags$sub("j"), sep = ""))),
                 p(em("* Fractions:"), "Are used to compute the Bayes factors and the posterior
                   probabilities of the hypotheses. Must be numbers between", HTML(paste("1/n",
                   tags$sub("j"), sep = "")), "and 1, where", HTML(paste("n", tags$sub("j"),
                   sep = "")), "is the sample size in group j. If no fractions are specified (i.e. if
                   the field is empty), then the app uses fractions of", HTML(paste("2/n",
                   tags$sub("j"), ".", sep = "")), "For details see Böing-Messing et al. (2015)."),
                 # http://shiny.rstudio.com/articles/tag-glossary.html
                 p(em("* Number of draws from the posterior distribution:"), "Computing the marginal
                   likelihood under an inequality constrained hypothesis involves sampling
                   from the posterior distribution of the group variances. The default number of draws
                   from the posterior is 100,000. For most applications this is sufficient. Larger
                   numbers result in more precise approximations of the marginal likelihood. In case
                   the number of inequality constraints under a certain hypothesis is large, one might
                   consider increasing the number of draws from the posterior."),
                 p(em("* Seed:"), "Seed for the Monte Carlo procedure described above. Is useful for
                   exactly replicating results. If no seed is specified (i.e. if the field is empty),
                   then the app uses the standard R procedure for automatically generating seeds."),
                 br(),
                 p(strong("Output")),
                 p(em("* Bayes factors:"), "For T > 1 hypotheses there are T*T Bayes factors. These are
                   all shown in the Bayes factors table. The table can be read as follows. For example,
                   the cell in row 2 and column 1 contains the Bayes factor", HTML(paste("B",
                   tags$sub(21), ",", sep = "")), "that is, the Bayes factor testing hypothesis 2
                   against hypothesis 1. Likewise, the cell in row 1 and column 4 contains the Bayes
                   factor", HTML(paste("B", tags$sub(14), ".", sep = "")), "The hypotheses are numbered
                   in the order that they are specified in the mandatory input. Note that",
                   HTML(paste("B", tags$sub(12), " = 1/B", tags$sub(21), ".", sep = "")), "Thus, the
                   Bayes factors above the diagonal are given by the inverse of the corresponding Bayes
                   factors below the diagonal. The Bayes factors on the diagonal are 1 because they
                   test a hypothesis against itself."),
                 p(em("* Posterior probabilities of the hypotheses:"), "Are computed using the marginal
                   likelihoods and the prior probabilities of the hypotheses. The posterior
                   probabilities necessarily sum to 1. The hypotheses are labeled in the order that
                   they are specified in the mandatory input (e.g. H1 refers to the first hypothesis
                   specified).")
        )
      )
    ),
  
    mainPanel(
      strong("Bayes factors"),
      tableOutput("bayesfactors"),
   
      strong("Posterior probabilities of the hypotheses"),
      tableOutput("posteriorprobabilities")
    )
  )
))