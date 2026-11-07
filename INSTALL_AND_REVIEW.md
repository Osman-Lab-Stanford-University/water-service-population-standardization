# Updating the repository on GitHub.com without Terminal

This v1.1.0 overlay adds 14 new files and replaces 17 existing files. The remaining 59 original files are unchanged. It is intended to update an existing v1.0.0 repository; do not upload this ZIP as a file.

1. Download this archive and double-click to extract it in Finder.
2. Open your existing GitHub repository in a web browser, select `Add file` → `Upload files` from the root directory.
3. Open the extracted `water-service-population-standardization/` folder in Finder. Drag **its contents**, not the parent folder itself, into the GitHub upload area. Select the contents including the `.github` folder (press Command+Shift+Period to show hidden files in Finder). The original repository files that are not supplied by the overlay should remain unchanged.
4. Make sure GitHub's upload preview shows the correct root-relative paths. New files should appear under `code/original/E27/`, `code/original/E28/`, `data/aggregates/`, `data/metadata/`, `docs/`, and `tests/`. Existing files such as `README.md`, `SHA256SUMS.txt` and the verification scripts should be replaced. Do not create a nested `water-service-population-standardization/` directory within the repository.
5. Choose **create a new branch** rather than committing directly to the default branch; use the message `Add E27–E28 outcome-model sensitivity (v1.1.0)` and propose changes. Review the pull request's `Files changed` tab for all 31 additions/updates, then merge after author/coauthor approval and passing checks. (If Github's interface rejects replacing existing files in a bulk upload, upload the 14 new files and use the built-in pencil editor to replace the 17 changed text files; do not delete unrelated files.)
6. Under the GitHub `Actions` tab, check that the repository validation workflow passes. The E27/E28 empirical scripts cannot be fully rerun from this public package because participant data and the original scorer/decoder are restricted.

Never upload the original `Karur Data.xlsx`, stage ZIPs, `PRIVATE_VERIFICATION`, `LOGS_PRIVATE`, decoder, scoring crosswalk, individual predictions, or survey identifiers. The public package has aggregate data and public source programs only.

For an initial upload of a repository rather than an update, use the separate full repository ZIP, extract it and upload the contents of its root folder into an empty GitHub repository.
