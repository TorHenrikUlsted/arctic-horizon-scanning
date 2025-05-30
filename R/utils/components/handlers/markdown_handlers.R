model_to_md <- function(model) {
  output <- cli::cli_format({
    # Turn off CLI colors temporarily
    old <- options(cli.num_colors = 1)
    on.exit(options(old))
    
    # Execute the function
    output <- capture.output(summary(model), type = "output")
  })
  
  # Initialize empty string for markdown text
  md_text <- ""
  
  # Process the output line by line
  in_coefficients <- FALSE
  current_section <- NULL
  
  for (line in output) {
    # Remove quotes and leading/trailing whitespace
    line <- gsub('"', "", trimws(line))
    
    # Skip empty lines and divider lines
    if (grepl("^\\s*$", line) || grepl("^[-*]{10,}", line)) next
    
    if (grepl("link function:", line)) {
      # Start new section
      section_name <- sub(" link function:.*", "", line)
      md_text <- paste0(md_text, "\n## ", section_name, " Parameters\n\n")
      # Add coefficient table header
      md_text <- paste0(
        md_text,
        "| Parameter | Estimate | Std. Error | t value | p-value |\n",
        "|-----------|-----------|------------|----------|----------|\n"
      )
      in_coefficients <- TRUE
    }
    # Handle coefficient lines
    else if (in_coefficients && grepl("^\\(Intercept\\)|^median|^[a-zA-Z]", line)) {
      # Check if this is actually a statistics line
      if (grepl("^No\\. of observations|^Degrees of Freedom|^Global Deviance|^AIC|^SBC", line)) {
        in_coefficients <- FALSE
        if (!grepl("Model Statistics", md_text)) {
          md_text <- paste0(md_text, "\n## Model Statistics\n\n")
        }
        clean_line <- gsub("\\s+", " ", line)
        md_text <- paste0(md_text, "- ", clean_line, "\n")
        next
      }
      
      # Skip significance code lines
      if (grepl("^Signif\\.", line)) next
      
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) >= 5) {
        param <- parts[1]
        est <- parts[2]
        se <- parts[3]
        t_val <- parts[4]
        
        # Handle p-value and stars
        p_val <- parts[5]
        if (length(parts) > 5) {
          if (grepl("\\*", p_val)) {
            p_val <- gsub("\\*", "\\\\*", p_val)
          } else {
            stars <- gsub("\\*", "\\\\*", parts[6])
            p_val <- paste0(p_val, " ", stars)
          }
        }
        
        md_text <- paste0(
          md_text,
          "| ", param, " | ",
          est, " | ",
          se, " | ",
          t_val, " | ",
          p_val, " |\n"
        )
      }
    } else if (grepl("^No\\. of observations|^Degrees of Freedom|^Global Deviance|^AIC|^SBC|^Residual Deg\\.|^at cycle:", line)) {
      if (!grepl("Model Statistics", md_text)) {
        md_text <- paste0(md_text, "\n## Model Statistics\n\n")
      }
      # Clean up the statistics line formatting
      clean_line <- gsub("\\s+", " ", line) # Replace multiple spaces with single space
      md_text <- paste0(md_text, "- ", clean_line, "\n")
    }
  }
  
  return(md_text)
}

#' Convert CLI styling to Markdown
#' @description Utility functions to convert CLI output to Markdown format
cli_to_markdown <- function(expr) {
  if (is.character(expr)) {
    output <- expr
  } else {
    output <- cli::cli_format({
      # Turn off CLI colors temporarily
      old <- options(cli.num_colors = 1)
      on.exit(options(old))
      
      # Execute the function
      output <- capture.output(expr(), type = "message")
    })
  }
  
  # Remove the extra quotes
  output <- gsub('^"|"$', '', output)
  
  # Only proceed if we have output
  if (length(output) == 0) {
    return("No markdown output")  # Return empty string if no output
  }
  
  # Remove empty strings at the end
  output <- output[nchar(trimws(output)) > 0]
  
  # Process each line
  processed <- sapply(output, function(line) {
    # Convert CLI formatting to Markdown
    line <- gsub("^[-] ", "- ", line)  # List items
    line <- gsub("^[*] ", "* ", line)  # Bullet points
    line <- gsub("^[#]+ (.+)$", "\\1", line)  # Remove cli header formatting
    
    # Convert CLI symbols to Markdown
    line <- gsub("\u2714", "✓", line)  # Tick
    line <- gsub("\u2718", "✗", line)  # Cross
    line <- gsub("\u25CF", "•", line)  # Bullet
    line <- gsub("\u26A0", "⚠", line)  # Warning
    
    # Convert CLI styling to Markdown using fixed patterns
    line <- gsub("\\{.emph ([^}]+)\\}", "*\\1*", line)  # Emphasis
    line <- gsub("\\{.strong ([^}]+)\\}", "**\\1**", line)  # Strong
    line <- gsub("\\{.code ([^}]+)\\}", "`\\1`", line)  # Code
    line <- gsub("\\{.field ([^}]+)\\}", "**\\1**", line)  # Field
    
    # Handle cli alert symbols
    line <- gsub("^ℹ ", "*Info:* ", line)  # Info alerts
    line <- gsub("^✔ ", "*Success:* ", line)  # Success alerts
    line <- gsub("^⚠ ", "*Warning:* ", line)  # Warning alerts
    line <- gsub("^✖ ", "*Error:* ", line)  # Error alerts
    
    # Remove [1] prefixes that R adds to output
    line <- gsub("^\\[\\d+\\] ", "", line)
    
    return(line)
  })
  
  # Post-process the entire output
  md_text <- paste(processed, collapse = "\n\n")
  
  # Convert CLI headers to Markdown headers
  md_text <- gsub("^### (.+)$", "### \\1", md_text, perl = TRUE)
  md_text <- gsub("^## (.+)$", "## \\1", md_text, perl = TRUE)
  md_text <- gsub("^# (.+)$", "# \\1", md_text, perl = TRUE)
  
  # Remove any trailing newlines
  md_text <- trimws(md_text)
  
  return(md_text)
}

# Helper function for severity colors in markdown
severity_to_markdown <- function(severity, text) {
  switch(severity,
         "PASS" = sprintf("✓ %s", text),
         "WARNING" = sprintf("⚠ %s", text),
         "FAIL" = sprintf("✗ %s", text),
         text
  )
}

# Helper function for converting CLI alerts to markdown
cli_alert_to_markdown <- function(type, text) {
  prefix <- switch(type,
                   "success" = "✓",
                   "warning" = "⚠",
                   "error" = "✗",
                   "info" = "ℹ",
                   ">"
  )
  sprintf("%s %s", prefix, text)
}


mdwrite <- function(source, text = NULL, data = NULL, image = NULL, image.out = "./outputs/images/image", image.which = NULL, device = "jpeg", open = "a", veb = TRUE) {
  if (!veb) {
    return(invisible())
  }
  
  create_file_if(source, keep = TRUE)
  
  if (is.data.table(data) || is.data.frame(data)) {
    data <- kable(data, format = "markdown")
  }
  
  if (grepl(";", text)) {
    split_str <- str_split(text, ";")[[1]]
    h_num <- split_str[[1]]
    h_text <- split_str[[2]]
  }
  
  if (!is.null(image)) {
    if (grepl("\\.", basename(image.out))) {
      image_base <- sub("\\..*$", "", basename(image.out))
      image.out <- paste0(dirname(image.out), "/", image_base)
    }
    
    plot_title <- basename(image.out)
    image.out <- paste0(image.out, ".", device)
    
    catn("Writing image to:", colcat(image.out, color = "output"))
    
    do.call(device, list(filename = image.out, width = 500, height = 500))
    
    if (!is.null(image.which)) {
      plot(image, which = image.which, pch = "*", cex = 2, main = h_text)
    } else {
      plot(image, pch = "*", cex = 2, main = h_text)
    }
    
    dev.off()
  }
  
  try(con <- file(source, open = open))
  sink(con, type = "output")
  
  if (!is.null(text)) {
    if (grepl(";", text)) catn(paste0(strrep("#", h_num), " ", h_text))
    if (!grepl(";", text)) catn(text)
  }
  if (!is.null(data)) print(data)
  if (!is.null(image)) catn(paste0("![", h_text, "]", "(", "../images/", basename(image.out), ")"))
  catn("  ")
  
  sink(type = "output")
  close(con)
}

write_wrangled_md <- function(dt.list, name, column = "scientificName") {
  fl <- length(unique(dt.list$formatted, by = column)[[column]])
  
  counts <- lapply(dt.list, function(dt) {
    if (!is.null(dt)) length(unique(dt, by = column)[[column]])
  })
  
  sum <- sum(unlist(counts[names(counts) != "formatted"]), na.rm = TRUE)
  
  md_dt <- data.table(formatted = fl, lost = (fl - sum))
  
  for (key in names(dt.list)) {
    if (key != "formatted" && !is.null(dt.list[[key]])) {
      md_dt[[key]] <- counts[[key]]
    }
  }
  
  setcolorder(md_dt, c(setdiff(names(md_dt), "lost"), "lost"))
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("2;", name),
    data = md_dt
  )
}