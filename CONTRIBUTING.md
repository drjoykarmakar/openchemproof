# Contributing

Use Python 3.12+ and `uv sync --locked --all-extras`. Run `uv run ruff check src benchmarks tests`,
`uv run ruff format --check src benchmarks tests`, `uv run mypy`, and `uv run pytest -q`.
Regenerate the fixture dataset and run `uv run chem-artifact-eval --oracle` after changing rules.

Each scientific rule needs a stated domain, primary reference, independent numerical example,
boundary/abstention tests, and a versioned policy change. Do not relabel fixtures merely to match
new predictions. State which observations are measured, source asserted, derived, or simulated.
Never replace missing spectra with calculated evidence and label it experimental.

Keep literature/source/scaffold groups together when designing held-out evaluations. Preserve
native source IDs, original values, source licenses, transformations and checksums. Data source
licenses in `external_data/` remain separate from the Apache-2.0 software license.

Do not commit API keys, experimental records you lack permission to share, local audit logs,
or bulk download caches. MCP HTTP defaults to loopback. Add authenticated deployment outside
this prototype before deliberately exposing it to remote clients.
