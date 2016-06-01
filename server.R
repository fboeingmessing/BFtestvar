source("Functions_13.R")

shinyServer(function(input, output) {
  
  output$bayesfactors <- renderTable({      
    s2 <- as.numeric(unlist(strsplit(input$s2, split = ",")))
    n <- as.numeric(unlist(strsplit(input$n, split = ",")))
    hypotheses <- unlist(strsplit(gsub(" ", "", input$hypotheses), split = ";"))
    log.BF <- input$log.BF
    if (is.na(input$prior.probabilities)) {
      prior.probabilities <- NA
    } else {
      prior.probabilities <- strsplit(gsub(" ", "", unlist(strsplit(input$prior.probabilities,
                                                                    split = ","))), split = "/")
      prior.probabilities <- sapply(prior.probabilities, function(x) if(length(x)==1) {as.numeric(x)}
                                    else {as.numeric(x[1])/as.numeric(x[2])})
    }
    if (is.na(input$b)) {
      b <- "default"
    } else {
      b <- strsplit(gsub(" ", "", unlist(strsplit(input$b, split = ","))), split = "/")
      b <- sapply(b, function(x) if(length(x)==1) {as.numeric(x)} else
                                                  {as.numeric(x[1])/as.numeric(x[2])})
    }
    nsim <- if (is.na(input$nsim)) {100000} else {input$nsim}
    seed <- if (is.na(input$seed)) {NA} else {input$seed}

    results <- shiny.function(s2, n, hypotheses, log.BF, prior.probabilities, b, nsim, seed)[[1]]
    },
    digits = 3, sanitize.rownames.function = function(x) paste('<b>',x,'</b>', sep ='')
  )
  
  output$posteriorprobabilities <- renderTable({      
    s2 <- as.numeric(unlist(strsplit(input$s2, split = ",")))
    n <- as.numeric(unlist(strsplit(input$n, split = ",")))
    hypotheses <- unlist(strsplit(gsub(" ", "", input$hypotheses), split = ";"))
    log.BF <- input$log.BF
    if (is.na(input$prior.probabilities)) {
      prior.probabilities <- NA
    } else {
      prior.probabilities <- strsplit(gsub(" ", "", unlist(strsplit(input$prior.probabilities,
                                                                    split = ","))), split = "/")
      prior.probabilities <- sapply(prior.probabilities, function(x) if(length(x)==1) {as.numeric(x)}
                                    else {as.numeric(x[1])/as.numeric(x[2])})
    }
    if (is.na(input$b)) {
      b <- "default"
    } else {
      b <- strsplit(gsub(" ", "", unlist(strsplit(input$b, split = ","))), split = "/")
      b <- sapply(b, function(x) if(length(x)==1) {as.numeric(x)} else
      {as.numeric(x[1])/as.numeric(x[2])})
    }
    nsim <- if (is.na(input$nsim)) {100000} else {input$nsim}
    seed <- if (is.na(input$seed)) {NA} else {input$seed}
    
    results <- shiny.function(s2, n, hypotheses, log.BF, prior.probabilities, b, nsim, seed)[[2]]
    },
    digits = 3, include.rownames = FALSE
  )

})