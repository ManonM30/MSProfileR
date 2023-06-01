species_identification <- function(DB_group, query_IM, top_peaks, top_hits, tolerance) {
  # Process the input files and retrieve the required data
  data <- process_file(DB_group, query_IM)
  DB_peaks <- get_best_mass_per_spectrum(data[[1]], top_peaks)
  query_peaks <- get_best_mass_per_spectrum(data[[2]], top_peaks)
  peaks_matching <- compare_peaks(data[[2]], DB_peaks, query_peaks, tolerance)
  hits_data <- classify_hits(peaks_matching, top_hits)
  
  # Assign global variables
  assign("DB_IM", data[[1]], envir = .GlobalEnv)
  assign("query_IM", data[[2]], envir = .GlobalEnv)
  assign("DB_peaks", DB_peaks, envir = .GlobalEnv)
  assign("query_peaks", query_peaks, envir = .GlobalEnv)
  assign("peak_matching", peaks_matching, envir = .GlobalEnv)
  assign("hits_table", hits_data, envir = .GlobalEnv)
  
  return(data.frame(hits_data))
}

process_file <- function(link_reference, link_query) {
  # Read the reference and query files
  reference <- read.csv(link_reference, header = FALSE, sep = ";")
  query <- read.csv(link_query, header = FALSE, sep = ";")
  
  colnames(query) <- query[1, ]
  rownames(query) <- query[, 1]
  query <- query[-1, -1]
  query <- format(query, scientific = FALSE)
  
  names <- reference[, 1]
  
  new_ref <- matrix(NA, nrow = nrow(reference), ncol = ncol(query), dimnames = list(names, colnames(query)))
  
  for (i in 1:length(names)) {
    # Find the corresponding row index
    row_index <- match(names[i], rownames(query))
    
    # Copy the corresponding row into the new table
    new_ref[i, ] <- as.character(query[row_index, ])
  }
  new_ref <- format(new_ref, scientific = FALSE)
  query <- edit_file_data(query)
  new_ref <- edit_file_data(new_ref)
  list(new_ref, query)
}

edit_file_data <- function(data) {
  # Transpose the data to get the correct format
  data <- t(data)
  
  # Get the current row names
  row_names <- rownames(data)
  
  # Create a new column containing the row names
  new_column <- matrix(row_names, nrow = nrow(data), ncol = 1)
  
  # Combine the new column with the existing data
  data <- cbind(new_column, data)
  
  # Rename the first column as "mass"
  colnames(data)[1] <- "mass"
  
  # Remove row names
  rownames(data) <- NULL
  
  # Convert to data frames
  data <- as.data.frame(data)
  
  # Return the combined data as a list
  return(data)
}

get_best_mass_per_spectrum <- function(data, peaks) {
  best_intensities_table <- apply(data, 2, function(x) head(sort(x, decreasing = TRUE), peaks)) # Retrieve the top intensities for each spectrum!
  
  mass_table <- matrix(NA, nrow = peaks, ncol = ncol(data))
  colnames(mass_table) <- colnames(data)
  
  for (i in 1:ncol(data)) {
    # Retrieve column i from both tables as vectors
    col_data <- as.vector(data[, i])
    col_table <- as.vector(best_intensities_table[, i])
    
    # Convert the vectors to numeric type
    col_data <- as.numeric(col_data)
    col_table <- as.numeric(col_table)
    
    # Get the indices of values in col_data that are present in col_table
    indices <- which(col_data %in% col_table)
    
    indices_sorted <- indices[order(col_data[indices], decreasing = TRUE)]  # Sort the indices based on the order of sorted values
    values_sorted <- col_data[indices_sorted]  # Sort the corresponding values
    first_column_sorted <- data[indices_sorted, 1]  # Retrieve the corresponding first column == mass column
    
    mass_table[, i] <- first_column_sorted
  }
  
  return(mass_table)
}

compare_peaks <- function(query, DB_peaks_table, query_peaks_table, tol) {
  # Initialization of the peak match matrix
  peak_match_matrix <- matrix(NA, nrow = ncol(query_peaks_table), ncol = ncol(DB_peaks_table),
                              dimnames = list(colnames(query_peaks_table), colnames(DB_peaks_table)))
  
  # Column-wise comparison
  for (i in 1:ncol(query_peaks_table)) {
    for (j in 1:ncol(DB_peaks_table)) {
      # Row-wise comparison
      match_count <- 0
      for (k in 1:nrow(query_peaks_table)) {
        for (l in 1:nrow(DB_peaks_table)) {
          if (abs(as.numeric(query_peaks_table[k, i]) - as.numeric(DB_peaks_table[l, j])) < tol) {
            match_count <- match_count + 1
          }
        }
      }
      peak_match_matrix[i, j] <- match_count
    }
  }
  peak_match_matrix <- peak_match_matrix[-1, -1]
  return(peak_match_matrix)
}

classify_hits <- function(matching_table, top_hits) {
  nb_rows <- nrow(matching_table)
  
  # Create column names for "score of matching" and "species matching"
  colnames_score <- paste0("score of matching_", rep(1:top_hits))
  colnames_species <- paste0("species matching_", rep(1:top_hits))
  
  # Create column names in the desired order
  colnames_combined <- c(rbind(colnames_score, colnames_species))
  colnames_combined <- c(colnames_combined[1:top_hits], colnames_combined[(top_hits + 1):(2 * top_hits)])
  
  results <- matrix(NA, nrow = nb_rows, ncol = top_hits * 2,
                    dimnames = list(rownames(matching_table), colnames_combined))
  
  for (i in 1:nb_rows) {
    query <- rownames(matching_table)[i]
    values_sorted <- matching_table[i, ][order(matching_table[i, ], decreasing = TRUE)][1:top_hits]
    colnames_sorted <- colnames(matching_table)[order(matching_table[i, ], decreasing = TRUE)][1:top_hits]
    
    # Extract specific parts from the column names
    regex <- ".*_(.*?)_.*"
    extracted_parts <- sub(regex, "\\1", colnames_sorted)
    
    # Store the values and column names in the results table
    for (j in 1:top_hits) {
      results[i, (j * 2) - 1] <- values_sorted[j]
      results[i, j * 2] <- paste(colnames_sorted[j], extracted_parts[j], sep = ":")
    }
  }
  
  query <- rownames(results)
  rownames(results) <- NULL
  results <- cbind(query, results)
  return(results)
}
