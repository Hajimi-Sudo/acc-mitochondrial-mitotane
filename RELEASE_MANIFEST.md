# Public code release manifest

Recommended repository name: `acc-mitochondrial-mitotane`.

## Include

- The Python and shell scripts in this directory
- `README.md`
- `requirements.txt`
- The project `README.md`
- `AGENTS.md` and the `00_contract/` research contracts when execution
  provenance is needed

## Exclude

- `data/`, raw GEO or TCGA files, and downloaded archives
- Patient-level clinical records and sample-level controlled-access data
- `results/`, `remote_evidence/`, and `publication_package/` outputs
- LaTeX build files such as `.aux`, `.bbl`, `.log`, `.out`, and `.pdf`
- `code/__pycache__/` and compiled Python files
- SSH keys, passwords, API keys, tokens, `.env` files, and local configuration

The scripts are remote-only by design. A public checkout documents the analysis
and can be executed after the user configures an isolated remote project root
with the required public datasets. It must not be interpreted as permission to
download restricted data or to run the analysis on the local manuscript
machine.

The public release contains executable analysis code and runtime requirements.
The repository does not include author metadata, submission metadata, or
manuscript source files.
