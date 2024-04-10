# DBConnection.R

#' @title connectToDatabase
#' 
#' @description Function to connect to the Serratus database.
#' 
#' @return A connection object to the Serratus database if successful, 
#' NULL otherwise.
#' 
#' @examples
#' conn <- connectToDatabase()
#' 
#' @importFrom DBI dbConnect
#' @importFrom RPostgres Postgres
#' 
#' @export
connectToDatabase <- function() {
  print("Connecting to Serratus…")
  tryCatch({
    conn <- DBI::dbConnect(RPostgres::Postgres(),
                           dbname = "summary",
                           host = "serratus-aurora-20210406.cluster-ro-ccz9y6yshbls.us-east-1.rds.amazonaws.com",
                           port = 5432,
                           user = "public_reader",
                           password = "serratus")
    print("Database Connected!")
    return(conn)
  }, error = function(e) {
    message("Unable to connect to Database: ", e$message)
    return(NULL)
  })
}

# [END]
