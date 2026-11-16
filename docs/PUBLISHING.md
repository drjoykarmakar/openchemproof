# Publish the reviewed repository

GitHub publishing was unavailable in the build environment: the integration was administratively
disabled and the CLI had no authenticated account. No remote repository, PR or release was created.

From a machine where you are authenticated as `drjoykarmakar`, review the data licenses and then:

```bash
cd openchemproof
git init -b main
git add .
git commit -m "Initialize OpenChemProof v0.1 research foundation"
gh repo create drjoykarmakar/openchemproof --public --source=. --remote=origin --push
```

If that repository already exists, do not replace it; inspect its contents and add changes through
a branch. Ensure your Git author identity is configured locally before committing. The archive
omits `.git`, environment directories, API keys, bulk NMR downloads and transient caches.

The included CI workflow can run after upload. Local passing checks are not remote CI results.
Add the repository URL to CITATION.cff only after creation. Obtain a release DOI through a real
archival service before adding one; no DOI has been reserved or invented.
