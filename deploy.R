################################################################################
#  deploy.R  —  Crichton Group Analysis Suite
#  Run this script from RStudio to deploy / update the app on shinyapps.io
#
#  FIRST-TIME SETUP (once per machine):
#    1. Create a free account at https://www.shinyapps.io
#    2. Go to: Account → Tokens → Add Token → Show → copy the setAccountInfo()
#       command and paste it into your RStudio console, then run it.
#    3. Fill in APP_NAME and (optionally) APP_TITLE below.
#    4. Source this file:  source("deploy.R")
#
#  SUBSEQUENT UPDATES:
#    Just source this file again — it re-deploys with your latest code.
################################################################################

# -- CONFIGURATION (edit these) ------------------------------------------------

APP_NAME  <- "analysis-suite"      # URL slug: <your-username>.shinyapps.io/<APP_NAME>
APP_TITLE <- "Crichton Group Analysis Suite"   # Display name in the dashboard

# ------------------------------------------------------------------------------

cat("\n======================================================\n")
cat(" Crichton Group Analysis Suite — Deployment Script\n")
cat("======================================================\n\n")


# -- 1. Check rsconnect is installed -------------------------------------------

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  cat("Installing rsconnect...\n")
  install.packages("rsconnect")
}
library(rsconnect)


# -- 2. Check shinyapps.io account is linked -----------------------------------

accounts <- rsconnect::accounts()
if (nrow(accounts) == 0) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  ⚠  No shinyapps.io account linked yet.\n\n")
  cat("  To link your account:\n")
  cat("    1. Log in at https://www.shinyapps.io\n")
  cat("    2. Click your username → Account → Tokens\n")
  cat("    3. Click 'Add Token' → 'Show' → copy the command\n")
  cat("    4. Paste and run it in the RStudio console\n")
  cat("    5. Then source this script again\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  stop("Account setup required — see instructions above.", call. = FALSE)
}

cat(sprintf("✓  Linked account: %s (%s)\n", accounts$name[1], accounts$server[1]))


# -- 3. Confirm we are running from the correct directory ----------------------

script_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()
)
setwd(script_dir)
cat(sprintf("✓  Working directory: %s\n", getwd()))


# -- 4. Check all required files are present -----------------------------------

required_files <- c(
  "app.R",
  "softmax_bca_improved.R",
  "tm_analysis_functions.R",
  "plot_akta_improved.R"
)

missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  ⚠  Missing files in this folder:\n")
  for (f in missing) cat(sprintf("      ✗  %s\n", f))
  cat("\n  Make sure all files listed above are in the same\n")
  cat("  folder as this deploy.R script, then try again.\n")
  cat("  (Rename app_WORKING.R → app.R if needed.)\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  stop("Missing required files — see list above.", call. = FALSE)
}

cat("✓  All required files found\n")


# -- 5. List files to be deployed ----------------------------------------------

deploy_files <- c(required_files)

# Include any additional helper scripts in the folder (e.g. future modules)
extra <- setdiff(
  list.files(pattern = "\\.R$", ignore.case = TRUE),
  c(deploy_files, "deploy.R")
)
if (length(extra) > 0) {
  deploy_files <- c(deploy_files, extra)
  cat(sprintf("✓  Extra R files included: %s\n", paste(extra, collapse = ", ")))
}

cat(sprintf("\n  Files to deploy (%d):\n", length(deploy_files)))
for (f in deploy_files) {
  size_kb <- round(file.size(f) / 1024, 1)
  cat(sprintf("    • %-40s  %6.1f KB\n", f, size_kb))
}


# -- 6. Confirm before deploying -----------------------------------------------

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("  Deploying to:  https://%s.shinyapps.io/%s\n",
            accounts$name[1], APP_NAME))
cat("  This will take a few minutes on first deploy.\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Prompt for confirmation if running interactively
if (interactive()) {
  answer <- readline("  Proceed? [y/N]: ")
  if (!tolower(trimws(answer)) %in% c("y", "yes")) {
    cat("  Deployment cancelled.\n")
    invisible(NULL)
    stop("Cancelled by user.", call. = FALSE)
  }
}


# -- 7. Deploy -----------------------------------------------------------------

cat("\n  Deploying...\n\n")

tryCatch({
  rsconnect::deployApp(
    appDir      = getwd(),
    appFiles    = deploy_files,
    appName     = APP_NAME,
    appTitle    = APP_TITLE,
    forceUpdate = TRUE,
    launch.browser = TRUE   # opens the live URL when done
  )

  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  ✓  Deployment successful!\n")
  cat(sprintf("  URL: https://%s.shinyapps.io/%s\n",
              accounts$name[1], APP_NAME))
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

}, error = function(e) {
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  ✗  Deployment failed.\n\n")
  cat("  Error message:\n")
  cat(sprintf("    %s\n", conditionMessage(e)))
  cat("\n  Common causes:\n")
  cat("    • Not connected to the internet\n")
  cat("    • Token expired — re-run the setAccountInfo() command\n")
  cat("      from shinyapps.io → Account → Tokens\n")
  cat("    • Free tier monthly hours exhausted (resets on 1st)\n")
  cat("    • A package failed to install on the remote server —\n")
  cat("      check the shinyapps.io dashboard logs for details\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  stop(conditionMessage(e), call. = FALSE)
})
