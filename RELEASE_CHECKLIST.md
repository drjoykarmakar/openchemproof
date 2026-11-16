# v0.1 delivery evidence

- Python 3.12 dependency installation and reproducible uv lockfile.
- Strict Pydantic v2 schemas, RDKit structure/MS calculations, physical optical correction,
  solvent-aware NMR checks and explicit inconclusive outcomes.
- Five MCP tools and two resources exercised through FastMCP's in-memory client.
- Actual stdio and Streamable HTTP subprocess smoke tests.
- Synthetic fixture generation: 500 cases, 100 in each requested category, seed 42.
- Oracle status agreement: 500/500; expected diagnostic flags also agree.
- Real public-source collection: 500 records; 100 per source category; source-specific licenses,
  raw selected evidence, scripts, hashes and missingness retained.
- Baseline/augmented API harness tested with mocked responses and real in-memory MCP calls.
  Paid/live LLM evaluation was not run; no LLM scores are claimed.
- Wheel/sdist build and local CLI examples exercised.
- Source and data scope documented in README, scientific scope, dataset card and source notes.

The final command-level validation results are in `validation.json`. Source downloads and local
checks are completed; publication and scientific validation remain distinct next steps.

GitHub publishing was unavailable because its integration was disabled and no CLI authentication
was present. No remote repository, PR, remote CI result, release DOI or published paper is claimed.
