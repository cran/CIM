#' Compositional Impact of Migration
#'
#' produce statistical indicators of the impact of migration on the socio-demographic
#' composition of an area. Three measures can be used: ratios, percentages and
#'  the Duncan index of dissimilarity. The input data files are assumed to be in an
#' origin-destination matrix format, with each cell representing a flow count
#' between an origin and a destination area. Columns are expected to represent origins,
#' and rows are expected to represent destinations. The first row and column are assumed to
#' contain labels for each area. See Rodríguez-Vignoli and Rowe (2018)
#' for technical details.
#'
#' @param ... 2 or more data frames, each containing an origin-destination migration matrix by
#' population attribute (i.e. age, sex, education, ethnicity, etc.).
#' Columns are expected to represent origins,
#' and rows are expected to represent destionations. The first row and column are assumed to
#' contain labels for each area.
#'
#' @param calculation a character, indicating the migration impact indicator selected to measure the socio-demographic
#' composition of an area. Users can type one of three options:
#' "ratio", "percentage" or "duncan".
#'
#' @param numerator a number, indicating the index number of the data frame to be used as
#' the numerator for the calculation. Type 1 to use the first data frame included in the function.
#' Type 2 to use the second data frame included in the function, and so on.
#'
#' @param denominator a number, indicating the index number of the data frame to be used as
#' the denominator for the calculation. Type 1 to use the first data frame included in the function.
#' Type 2 to use the second data frame included in the function, and so on.
#' Note the numerator data frame must differ from the denominator data frame.
#'
#' @param DuncanAll, a logical argument. If calculation = "Duncan", this logical argument must be specified.
#' The Duncan index measures the dissimilarity in the spatial
#' distribution of a chosen group (first data frame in the function) against a reference category
#' as specified by the "DuncanAll" argument.
#' If TRUE, the reference category is the sum of all data frames,
#' except for the first data frame included in the function (i.e. chosen group).
#' If FALSE, a specific data frame must be specified to be the reference group. See Duncan and Duncan (1955)
#' for details on the Duncan index, and Rodríguez-Vignoli and Rowe (2017a, b) for an empirical application
#' of the CIM using the Duncan index.
#'
#' @param rest, a logical argument. If calculation = "Duncan", this argument must be specified.
#' It enables a special calculation of the CIM, for a particular area (e.g. the Greater London Metropolitan Area), and
#' the rest of spatial units comprising a country. To correctly compute the CMI, these spatial units
#' need to be amalgamated and included as a single column/row in the matrix - labelled "Rest of the country"
#' (e.g. Rest of the UK). If TRUE, the column/row of the "Rest of the country" is considered for the calculation
#' and is excluded from the denominator of the duncan index.
#' If FALSE, the "Rest of the country" column/row is included in the denominator, producing the wrong results.
#'
#' @return an object containing:
#'
#' @return for the "ratio" and "percentage" calculation options:
#' @return num_results: a data frame containing nine area-level indicators: the Factual Value (FV), Counterfactual Value (CFV),
#' Compositional Impact of Migration (CIM), Compositional Impact of Migration Percentage Change (CIM_PC),
#' Diagonal Cell Indicator (DIAG), Compositional Impact of Migration for Inflows (CIM_I),
#' Compositional Impact of Migration for Outflows (CIM_O), CIM_I as a percentage of CMI (CIM_I_PC), and
#' CIM_O as a percentage of CMI (CIM_O_PC)
#'
#' @return for the "duncan" calculation option:
#'
#' @return duncan_results: a data frame, containing
#' the Factual Value of the Area-Specific Share (ASFVShare_cg), and the Counterfactual Value of the Area-Specific Share (ASCFVShare_cg) for the chosen group;
#' the Factual Value of the Area-Specific Share (ASFVShare_ref) and the Counterfactual Value of the Area-Specific Share (ASCFVShare_ref) for the reference group;
#' the Area-Specific Share Factual Value Difference between the ASFVShare_cg and ASFVShare_ref (ASShareFV_diff);
#' and the Area-Specific Share Counterfactual Value Difference between the ASCFVShare_cg and ASCFVShare_ref (ASShareCFV_diff).
#' The chosen group corresponds to the first data frame in the function. See above the argument "DuncanAll"
#' to specify the reference category.
#'
#' @return duncan_index: a numeric value, indicating the Duncan Index of dissimilarity for the chosen group.
#'
#' @examples
#' ## Read in the two data.frames included in the package
#' m <- male
#' f <- female
#'
#' ## Run the function using "ratio" calculation
#' CIM.ratio <- CIM(m, f, calculation = "ratio", numerator = 1, denominator = 2)
#' ## Print the resulted data.frame
#' CIM.ratio
#'
#' ## Run the function using "percentage" calculation
#' CIM.percentage <- CIM(m, f, calculation = "percentage", numerator = 1, denominator = 2)
#' ## See the resulted data.frame
#' CIM.percentage
#'
#' ## For the Duncan index, we compute impact of internal migration on the spatial pattern of
#' ## residential age segregation of people age 65 and over in the
#' ## local authority districts of Greater London using 2011 census data.
#' ## Chosen group: people aged 65 and over.
#' ## Reference category: the rest of age groups.
#' ## For this example, this group is people aged pop1-14, 15-29, 30-14 and 45-64).
#' CIM.duncan <- CIM(pop65over, pop1_14, pop15_29, pop30_44, pop45_64,
#' calculation = "duncan", numerator = 1, DuncanAll= TRUE)
#' CIM.duncan$duncan_results
#' CIM.duncan$duncan_index
#'
#'
#' @references
#'
#' Duncan, O.D. and Duncan, B., 1955. A methodological analysis of segregation indexes.
#' American sociological review, 20(2), pp.210-217.
#'
#' Rodríguez-Vignoli, J.R. and Rowe, F., 2017a. ¿Contribuye la migración interna a
#' reducir la segregación residencial?: el caso de Santiago de Chile 1977-2002.
#' Revista Latinoamericana de Población, (21), pp.7-46.
#'
#' Rodríguez-Vignoli, J.R. and Rowe, F., 2017b. The Changing Impacts of Internal Migration
#' on Residential Socio-Economic Segregation in the Greater Santiago. 28th International
#' Population Conference of the International Union for the Scientific Study of Population (IUSSP),
#'  Cape Town, South Africa.
#'
#' Rodríguez-Vignoli, J. and Rowe, F., 2018. How is internal migration reshaping
#' metropolitan populations in Latin America? A new method and new evidence.
#' Population studies, 72(2), pp.253-273. doi.org/10.1080/00324728.2017.1416155
#'
#'
#' @importFrom  utils head
#'
#' @export

CIM <- function(..., calculation, numerator, denominator, DuncanAll= TRUE, rest = TRUE) {
  # Put all the data.frames in a list
  args <- list(...)
  # Create an empty list which will hold the new calculated data.frames
  b <- list()
  # Create a counter
  count_tables <- 1
  for (a in args) {
    a <- cbind(a, totalRow= rowSums(a))
    a <- rbind(a, totalCol = colSums(a))
    b[[count_tables]] <- a
    count_tables = count_tables+1
  }

  if(calculation == "ratio"){
    # Calculate the ratio between the specified numerator and denominator
    c <- b[[numerator]]/b[[denominator]] *100
    # Transpose the last row
    c <- cbind(c, t(c[nrow(c),]))
    # Rename the last two columns
    colnames(c)[ncol(c)-1] <- "FV"
    colnames(c)[ncol(c)] <- "CFV"
    # Calculate the CIM index (FV-CFV)
    c$CIM <- c$FV - c$CFV
    # Calculate the CIM inflows index (FV-MII[ii])
    for (i in 1:nrow(c)) {
      c$CIM_I[i] <- c$FV[i] - c[i,i]
    }
    # Calculate the CIM outflows index (MII[ii]- CFV)
    for (i in 1:nrow(c)) {
      c$CIM_O[i] <- c[i,i] -  c$CFV[i]
    }
    # Extract the diagonal of the dataframe
    # Extract the diagonal of the dataframe
    for (i in 1:nrow(c)) {
      c$DIAG[i] <- c[i,i]
    }
    # Calculate Percentage change CIM index ([FV-CFV] / CFV)
    c$CIM_PC <- c$CIM / c$CFV *100
    # Calculate percentage of CMI_I over CMI
    c$CIM_I_PC <- c$CIM_I / c$CIM *100
    # Calculate percentage of CMI_O over CMI
    c$CIM_O_PC <- c$CIM_O / c$CIM *100
    # replace NA values with zeroes
    c[is.na(c)] <- 0
    # Return the columns specified as the result of this option of the function
    res <- list(num_results = c[c("FV", "CFV", "CIM", "CIM_PC","DIAG",
                                  "CIM_I", "CIM_O", "CIM_I_PC", "CIM_O_PC")])
    return(res)
  }

  else if (calculation == "percentage") {
    # Sum all the data.frames
    m_sum <- Reduce('+', b)
    # Calculate the percentage of the specified numerator
    c <- b[[numerator]]/m_sum *100
    # Transpose the last row
    c <- cbind(c, t(c[nrow(c),]))
    # Rename the last two columns
    colnames(c)[ncol(c)-1] <- "FV"
    colnames(c)[ncol(c)] <- "CFV"
    # Calculate the CIM index (FV-CFV)
    c$CIM <- c$FV - c$CFV
    # Calculate the CIM inflows index (FV-MII[ii])
    for (i in 1:nrow(c)) {
      c$CIM_I[i] <- c$FV[i] - c[i,i]
    }
    # Calculate the CIM outflows index (MII[ii]- CFV)
    for (i in 1:nrow(c)) {
      c$CIM_O[i] <- c[i,i] -  c$CFV[i]
    }
    # Calculate Percentage change CIM index ([FV-CFV] / CFV)
    c$CIM_PC <- c$CIM / c$CFV *100
    # Return the columns specified as the result of this option of the function
    res <- list(num_results = c[c("FV", "CFV", "CIM", "CIM_I", "CIM_O", "CIM_PC")])
    return(res)
  }

  else if (calculation == "duncan"){
    if (DuncanAll){
      # Sum all the data.frames except the specified numerator
      sum_matr <- Reduce('+', b[-c(numerator)])
      # Transpose the last row
      sum_matr <- cbind(sum_matr, t(sum_matr[nrow(sum_matr),]))
      c <- b[[numerator]]
      # Transpose the last row
      c <- cbind(c, t(c[nrow(c),]))
      if (rest) {
        # Calculate Factual Value and Counterfactual Value for numerator
        # Ignore the last row
        c$ASFVShare_cg <- c$totalRow / sum(head(c$totalRow, -2))
        c$ASCFVShare_cg <- c$totalCol / sum(head(c$totalCol, -2))
        # Calculate Factual Value and Counterfactual Value for the rest of the data.frames
        # Ignore the last row
        sum_matr$ASFVShare_ref <- sum_matr$totalRow / sum(head(sum_matr$totalRow, -2))
        sum_matr$ASCFVShare_ref <- sum_matr$totalCol / sum(head(sum_matr$totalCol, -2))
        # Calculate Duncan FVs and CFVs for each row
        c$ASShareFV_diff <- abs(c$ASFVShare_cg - sum_matr$ASFVShare_ref)
        c$ASShareCFV_diff <- abs(c$ASCFVShare_cg - sum_matr$ASCFVShare_ref)
        # Sum Duncan FVs and CFVs
        sumDuncanFV <- sum(head(c$ASShareFV_diff, -2))
        sumDuncanCFV <- sum(head(c$ASShareCFV_diff, -2))
        # Duncan Index
        DuncanIndex <- (sumDuncanFV/2) - (sumDuncanCFV/2)
        # remove the rest of the country row which is always the  second row from the bottom of the dataframe
        c <- c[-c(nrow(c)-1),]
        sum_matr <- sum_matr[-c(nrow(sum_matr)-1),]
        # Return the columns specified as the result of this option of the function
        res <- list(duncan_results = cbind(c[c("ASFVShare_cg", "ASCFVShare_cg")],
                                           sum_matr[c("ASFVShare_ref", "ASCFVShare_ref")],
                                           c[c("ASShareFV_diff", "ASShareCFV_diff")]),
                    duncan_index = DuncanIndex)
        return(res)
      }
      # Calculate Factual Value and Counterfactual Value for numerator
      # Ignore the last row
      c$ASFVShare_cg <- c$totalRow / sum(c$totalRow[1:nrow(c)-1])
      c$ASCFVShare_cg <- c$totalCol / sum(c$totalCol[1:nrow(c)-1])
      # Calculate Factual Value and Counterfactual Value for the rest of the data.frames
      # Ignore the last row
      sum_matr$ASFVShare_ref <- sum_matr$totalRow / sum(sum_matr$totalRow[1:nrow(sum_matr)-1])
      sum_matr$ASCFVShare_ref <- sum_matr$totalCol / sum(sum_matr$totalCol[1:nrow(sum_matr)-1])
      # Calculate Duncan FVs and CFVs for each row
      c$ASShareFV_diff <- abs(c$ASFVShare_cg - sum_matr$ASFVShare_ref)
      c$ASShareCFV_diff <- abs(c$ASCFVShare_cg - sum_matr$ASCFVShare_ref)
      # Sum Duncan FVs and CFVs
      sumDuncanFV <- sum(c$ASShareFV_diff)
      sumDuncanCFV <- sum(c$ASShareCFV_diff)
      # Duncan Index
      DuncanIndex <- (sumDuncanFV/2) - (sumDuncanCFV/2)
      # remove the rest of the country row which is always the  second row from the bottom of the dataframe
      c <- c[-c(nrow(c)-1),]
      sum_matr <- sum_matr[-c(nrow(sum_matr)-1),]
      # Return the columns specified as the result of this option of the function
      res <- list(duncan_results = cbind(c[c("ASFVShare_cg", "ASCFVShare_cg")],
                                         sum_matr[c("ASFVShare_ref", "ASCFVShare_ref")],
                                         c[c("ASShareFV_diff", "ASShareCFV_diff")]),
                  duncan_index = DuncanIndex)
      return(res)

    }
    if (rest){
      ratios <- list()
      i <- 1
      for (a in b) {
        a <- cbind(a, t(a[nrow(a),]))
        # ignore the last row
        a$ASFVShare <- a$totalRow / sum(head(a$totalRow, -2))
        a$ASCFVShare <- a$totalCol / sum(head(a$totalCol, -2))
        ratios[[i]] <- a
        i <- i+1
      }
      c <- ratios[[numerator]]
      d <- ratios[[denominator]]
      # Calculate Duncan FVs and CFVs for each row
      c$ASShareFV_diff <- abs(c$ASFVShare - d$ASFVShare)
      c$ASShareCFV_diff <- abs(c$ASCFVShare - d$ASCFVShare)
      # Sum Duncan FVs and CFVs
      sumDuncanFV <- sum(head(c$ASShareFV_diff, -2))
      sumDuncanCFV <- sum(head(c$ASShareCFV_diff, -2))
      # Duncan Index
      DuncanIndex <- (sumDuncanFV/2) - (sumDuncanCFV/2)
      # rename columns to reflect which is the chisen group and which is the reference
      names(c)[names(c) == "ASFVShare"] <- "ASFVShare_cg"
      names(c)[names(c) == "ASCFVShare"] <- "ASCFVShare_cg"
      names(d)[names(d) == "ASFVShare"] <- "ASFVShare_ref"
      names(d)[names(d) == "ASCFVShare"] <- "ASCFVShare_ref"
      # remove the rest of the country row which is always the  second row from the bottom of the dataframe
      c <- c[-c(nrow(c)-1),]
      d <- d[-c(nrow(d)-1),]
      # Return the columns specified as the result of this option of the function
      res <- list(duncan_results = cbind(c[c("ASFVShare_cg", "ASCFVShare_cg")],
                                         d[c("ASFVShare_ref", "ASCFVShare_ref")],
                                         c[c("ASShareFV_diff", "ASShareCFV_diff")]),
                  duncan_index = DuncanIndex)
      return(res)
    }
    ratios <- list()
    i <- 1
    for (a in b) {
      a <- cbind(a, t(a[nrow(a),]))
      # ignore the last row
      a$ASFVShare <- a$totalRow / sum(a$totalRow[1:nrow(a)-1])
      a$ASCFVShare <- a$totalCol / sum(a$totalCol[1:nrow(a)-1])
      ratios[[i]] <- a
      i <- i+1
    }
    c <- ratios[[numerator]]
    d <- ratios[[denominator]]
    # Calculate Duncan FVs and CFVs for each row
    c$ASShareFV_diff <- abs(c$ASFVShare - d$ASFVShare)
    c$ASShareCFV_diff <- abs(c$ASCFVShare - d$ASCFVShare)
    # Sum Duncan FVs and CFVs
    sumDuncanFV <- sum(c$ASShareFV_diff)
    sumDuncanCFV <- sum(c$ASShareCFV_diff)
    # Duncan Index
    DuncanIndex <- (sumDuncanFV/2) - (sumDuncanCFV/2)
    # rename columns to reflect which is the chisen group and which is the reference
    names(c)[names(c) == "ASFVShare"] <- "ASFVShare_cg"
    names(c)[names(c) == "ASCFVShare"] <- "ASCFVShare_cg"
    names(d)[names(d) == "ASFVShare"] <- "ASFVShare_ref"
    names(d)[names(d) == "ASCFVShare"] <- "ASCFVShare_ref"
    # remove the rest of the country row which is always the  second row from the bottom of the dataframe
    c <- c[-c(nrow(c)-1),]
    d <- d[-c(nrow(d)-1),]
    # Return the columns specified as the result of this option of the function
    res <- list(duncan_results = cbind(c[c("ASFVShare_cg", "ASCFVShare_cg")],
                                       d[c("ASFVShare_ref", "ASCFVShare_ref")],
                                       c[c("ASShareFV_diff", "ASShareCFV_diff")]),
                duncan_index = DuncanIndex)
    return(res)
  }
}
