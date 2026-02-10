#!/usr/bin/env Rscript
# Set CRAN repo to Posit Public Package Manager (P3M) binary repos for supported distros
.set_repository_options <- function(repo = c("bioconductor", "cran"), bioc_version = NULL, latest_cran = FALSE) {
  repo <- match.arg(repo)

  # Get OS information
  os_release <- tryCatch(readLines("/etc/os-release", warn = FALSE), error = function(e) NULL)
  if (is.null(os_release)) return(invisible(NULL))

  get_field <- function(key) {
    m <- grep(paste0('^', key, '='), os_release, value = TRUE)
    if (length(m) == 0) return(NA_character_)
    val <- sub(paste0('^', key, '=(.*)$'), '\\1', m[1])
    val <- gsub('"', '', val)
    trimws(tolower(val))
  }

  id <- get_field('ID')
  id_like <- get_field('ID_LIKE')
  vid <- get_field('VERSION_ID')
  ver_major <- if (!is.na(vid) && nzchar(vid)) as.integer(sub('\\..*$', '', vid)) else NA_integer_

  # Determine distro codename
  distro_codename <- NA_character_
  if (!is.na(id) && id == 'ubuntu') {
    if (!is.na(ver_major) && ver_major == 22) distro_codename <- 'jammy'
    if (!is.na(ver_major) && ver_major == 24) distro_codename <- 'noble'
  }
  if (is.na(distro_codename) && !is.na(id) && id == 'debian') {
    if (!is.na(ver_major) && ver_major == 12) distro_codename <- 'bookworm'
    if (!is.na(ver_major) && ver_major == 13) distro_codename <- 'trixie'
  }
  if (is.na(distro_codename)) {
    rhel_ids <- c('rhel','redhat','rocky','almalinux')
    if ((!is.na(id) && id %in% rhel_ids) || (!is.na(id_like) && grepl('rhel', id_like))) {
      if (!is.na(ver_major)) distro_codename <- paste0('rhel', ver_major)
    }
  }

  if (is.na(distro_codename)) {
    warning("Unsupported or undetected Linux distribution - no repository changes applied.")
    return(invisible(NULL))
  }

  # Mapping tables
  bioc_to_p3m_date <- c(
    "3.5"  = "2017-10-30",
    "3.6"  = "2018-04-30",
    "3.7"  = "2018-10-29",
    "3.8"  = "2019-05-02",
    "3.9"  = "2019-10-28",
    "3.10" = "2020-04-27",
    "3.11" = "2020-10-28",
    "3.16" = "2023-04-26",
    "3.17" = "2023-10-25",
    "3.18" = "2024-05-01",
    "3.19" = "2024-10-30",
    "3.20" = "2025-04-16",
    "3.21" = "2025-10-30"
  )

  bioc_to_r_local <- c(
    "3.5"  = "3.4",
    "3.6"  = "3.4",
    "3.7"  = "3.5",
    "3.8"  = "3.5",
    "3.9"  = "3.6",
    "3.10" = "3.6",
    "3.11" = "4.0",
    "3.12" = "4.0",
    "3.13" = "4.1",
    "3.14" = "4.1",
    "3.15" = "4.2",
    "3.16" = "4.2",
    "3.17" = "4.3",
    "3.18" = "4.3",
    "3.19" = "4.4",
    "3.20" = "4.4",
    "3.21" = "4.5",
    "3.22" = "4.5",
    "3.23" = "4.6"
  )

  # Get current R version (using base R only - utils::head not available in .Rprofile)
  cur_r <- paste(strsplit(as.character(getRversion()), "\\.")[[1]][1:2], collapse = '.')

  # Helper: find latest Bioconductor version for an R version
  bioc_for_r_local <- function(rver) {
    matches <- names(bioc_to_r_local)[bioc_to_r_local == rver]
    if (!length(matches)) return(NA_character_)
    parts_list <- lapply(matches, function(b) as.integer(strsplit(b, "\\.")[[1]]))
    maj <- sapply(parts_list, `[`, 1)
    minor <- sapply(parts_list, function(x) if (length(x) >= 2) x[2] else 0)
    ord <- order(maj, minor)
    matches[ord[length(ord)]]
  }

  # Main logic
  snapshot_date <- NULL

  if (repo == 'cran') {
    # CRAN logic
    if (latest_cran) {
      # Use latest when explicitly requested
      repo_url <- sprintf('https://packagemanager.posit.co/cran/__linux__/%s/latest', distro_codename)
    } else {
      # Choose snapshot based on current R version
      inferred_bioc <- bioc_for_r_local(cur_r)
      if (!is.na(inferred_bioc) && inferred_bioc %in% names(bioc_to_p3m_date)) {
        snapshot_date <- bioc_to_p3m_date[[inferred_bioc]]
        repo_url <- sprintf('https://packagemanager.posit.co/cran/__linux__/%s/%s', distro_codename, snapshot_date)
      } else {
        # No snapshot available, use latest
        repo_url <- sprintf('https://packagemanager.posit.co/cran/__linux__/%s/latest', distro_codename)
      }
    }

  } else {
    # Bioconductor logic (latest_cran is ignored for bioconductor)
    options(BioC_mirror = 'https://packagemanager.posit.co/bioconductor/latest')
    options(BIOCONDUCTOR_CONFIG_FILE = 'https://packagemanager.posit.co/bioconductor/latest/config.yaml')

    # Determine which Bioconductor version to use
    if (!is.null(bioc_version) && nzchar(bioc_version)) {
      # User specified a version - validate it exists
      if (!(bioc_version %in% names(bioc_to_r_local))) {
        message(sprintf("Unknown Bioconductor version '%s' - no repository changes applied.", bioc_version))
        return(invisible(NULL))
      }

      # Check if it matches current R version
      expected_r <- bioc_to_r_local[[bioc_version]]
      if (expected_r != cur_r) {
        message(sprintf("Warning: Bioconductor %s requires R %s but you have R %s - no repository changes applied.",
                        bioc_version, expected_r, cur_r))
        return(invisible(NULL))
      }
    } else {
      # Infer Bioconductor version from current R version
      bioc_version <- bioc_for_r_local(cur_r)
      if (is.na(bioc_version)) {
        message(sprintf("No Bioconductor version found for R %s - no repository changes applied.", cur_r))
        return(invisible(NULL))
      }
    }

    # Set environment variable for validated version
    Sys.setenv(R_BIOC_VERSION = bioc_version)

    # Choose snapshot date or latest
    if (bioc_version %in% names(bioc_to_p3m_date)) {
      snapshot_date <- bioc_to_p3m_date[[bioc_version]]
      repo_url <- sprintf('https://packagemanager.posit.co/cran/__linux__/%s/%s', distro_codename, snapshot_date)
    } else {
      # Bioc version found but no snapshot - use latest
      repo_url <- sprintf('https://packagemanager.posit.co/cran/__linux__/%s/latest', distro_codename)
    }
  }

  # Set repository
  options(repos = c(CRAN = repo_url))

  # Display informative message
  bioc_env <- Sys.getenv('R_BIOC_VERSION')
  if (repo == 'bioconductor' && nzchar(bioc_env)) {
    message(sprintf("R %s: Bioconductor %s selected; CRAN repo: %s", cur_r, bioc_env, repo_url))
  } else {
    message(sprintf("R %s: CRAN repo: %s", cur_r, repo_url))
  }

  invisible(repo_url)
}

.set_repository_options()
