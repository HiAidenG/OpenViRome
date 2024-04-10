# SQLInterface.R

# All the functions in this module are internal
# and should not be called directly by the user.

#' @title constructQuery
#' 
#' @description Constructs a SQL query based on the provided filters.
#' 
#' @param SQLtable The name of the SQL table to query.
#' @param filters A named list of filters to apply to the query.
#' 
#' @return A SQL query string.
#' 
#' @keywords internal
constructQuery <- function(SQLtable = NULL, filters = NULL) {
  if (is.null(SQLtable) || is.null(filters)) {
    stop("Please provide both the SQL table and filters.")
  }

  # Prepare the SELECT clause
  selectClause <- constructSelect(filters = filters)

  # Prepare the WHERE clause
  whereClause <- constructWhere(filters = filters)

  # Construct the full query
  query <- paste(selectClause, "FROM", SQLtable, whereClause)
  
  print(query)

  return(query)
}

#' @title executeQuery
#' 
#' @description Executes a SQL query on the provided connection.
#' 
#' @param conn A DBI connection object.
#' 
#' @param query A SQL query string.
#' 
#' @return The result of the query.
#' 
#' @keywords internal
executeQuery <- function(conn = NULL, query = NULL) {
  if (is.null(query)) {
    stop("Please provide a query to execute.")
  }

  result <- DBI::dbGetQuery(conn, query)

  return(result)
}

#' @title constructSelect
#' 
#' @description Constructs the SELECT clause of a SQL query.
#' 
#' @param filters A named list of filters to apply to the query.
#' 
#' @return A SELECT clause string.
#' 
#' @keywords internal
constructSelect <- function(filters = NULL) {
  if (is.null(filters)) {
    stop("Please provide filters.")
  }

  # All names in the filters are columns to select, irrespective of their values
  columns <- names(filters)

  # Construct the SELECT clause
  selectClause <- sprintf("SELECT %s", paste(columns, collapse = ", "))

  return(selectClause)
}

#' @title constructWhere
#' 
#' @description Constructs the WHERE clause of a SQL query.
#' 
#' @param filters A named list of filters to apply to the query.
#' 
#' @return A WHERE clause string.
#' 
#' @keywords internal
constructWhere <- function(filters = NULL) {
  if (is.null(filters)) {
    stop("Please provide filters.")
  }

  # If a filter's value is NA, ignore it
  where <- filters[!is.na(filters)]

  # First check if any inequality parameters are present
  supportedInequalityParams <- retrieveParamInfo(paramType = 'inequalityQuerying')
  if (any(names(where) %in% names(supportedInequalityParams))) {
    requestedInequalities <- where[names(where) %in% names(supportedInequalityParams)]

    inequalityClause <- processInequalityArgs(requestedInequalities = requestedInequalities)

    where <- where[!(names(where) %in% names(supportedInequalityParams))]

  }
  # Construct the rest of the WHERE clause
  whereClause <- processWhereArgs(whereArgs = where, index = 1)

  # Add the inequality clause, if present
  if (exists("inequalityClause")) {
    whereClause <- paste(whereClause, "AND", inequalityClause)
  }

  return(paste("WHERE", whereClause))
}

#' @title processInequalityArgs
#' 
#' @description Processes the inequality parameters in the WHERE clause.
#' 
#' @param requestedInequalities A named list of inequality parameters.
#' 
#' @return An inequality clause string.
#' 
#' @keywords internal
processInequalityArgs <- function(requestedInequalities = NULL) {
  # Extract the inequality parameters
  supportedInequalityParams <- retrieveParamInfo(paramType = 'inequalityQuerying')

  # Construct the inequality conditions
  inequalityConditions <- sapply(names(requestedInequalities), function(param) {
    inequality <- supportedInequalityParams[[param]]
    sprintf("%s %s %s", param, inequality, requestedInequalities[[param]])
  })

  return(paste(inequalityConditions, collapse = " AND "))
}


#' @title processWhereArgs
#' 
#' @description Proccesses the WHERE clause arguments. Helper for 
#'  constructWhere.
#' 
#' @param whereArgs A named list of filters to apply to the query.
#' 
#' @param index The index of the current filter.
#' 
#' @return A WHERE clause string.
#' 
#' @keywords internal
processWhereArgs <- function(whereArgs = NULL, index = 1) {
  # Base case: if there are no more filters, return an empty string
  if (index > length(whereArgs)) {
    return("")
  }

  # Extract the current filter's name and values
  filterNames <- names(whereArgs)[index]
  filterValues <- whereArgs[[index]]

  condition <- sprintf("%s IN ('%s')", filterNames, paste(filterValues, collapse = "', '"))

  # Recursively process the next filters
  nextConditions <- processWhereArgs(whereArgs = whereArgs, index = index + 1)

  # Combine the current condition with the next conditions using AND, if necessary
  if (nchar(nextConditions) > 0) {
    return(sprintf("%s AND %s", condition, nextConditions))
  } else {
    return(condition)
  }
}

# [END]