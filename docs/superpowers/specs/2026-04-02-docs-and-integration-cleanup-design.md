# Docs And Integration Cleanup Design

## Purpose

This design defines the next cleanup slice after the `docs/superpowers` index work.

The goal is to bundle:

- historical docs path portability cleanup
- targeted codebase integration cleanup around the recently implemented analysis and plotting work

This is still a cleanup slice, not a feature-building slice. It should improve clarity, portability, and discoverability of work that already landed.

## Core Direction

Cleanup should be docs-led.

The order is:

1. normalize historical docs so they are portable in the repository
2. perform a narrow integration cleanup around the completed code/docs surface

This keeps the cleanup low-risk and prevents it from turning into an unfocused refactor.

## Scope

This bundled cleanup should include two phases.

### Phase A: Historical docs normalization

Normalize machine-specific path leakage in tracked historical docs under:

- `docs/superpowers/specs/`
- `docs/superpowers/plans/`

This means:

- replace machine-local absolute paths with repo-relative or plain repo paths where appropriate
- preserve the original meaning of the documents
- avoid changing the substantive design or implementation content

This phase is about portability and maintainability of the internal planning trail, not rewriting history.

### Phase B: Codebase integration cleanup

Perform a targeted cleanup around the already-implemented analysis and plotting work.

This should focus on:

- public surface discoverability
- docs cross-linking
- lightweight consistency fixes around recently added helpers

It should not introduce new analysis features or broader architecture changes.

## Out Of Scope

This cleanup slice should not:

- redesign existing features
- add new analysis workflows
- rename large parts of the codebase for style reasons alone
- refactor internals that are not causing integration friction
- rewrite historical docs for prose polish beyond what is needed for path portability

## Phase A Design: Historical Docs Portability

The tracked historical plan/spec docs currently contain many machine-local absolute paths.

That is acceptable during active execution, but once the docs are committed as repo artifacts they should be portable enough to read in:

- GitHub
- local clones on other machines
- future branches

### Recommended path policy

For historical docs:

- use repo-relative markdown links where links are useful
- use plain repo-relative paths in code blocks and command examples when a clickable link is not needed
- avoid machine-local absolute paths unless a document explicitly needs to show an environment-specific example

### Editing rule

This phase should make the smallest possible content changes:

- path normalization only
- no semantic rewrites unless a path reference becomes misleading when converted

## Phase B Design: Integration Cleanup

The codebase side should prioritize discoverability over internal churn.

### Primary targets

1. public surface and discoverability
   - ensure docs point users to the new helpers and workflows clearly
   - ensure exports/import surfaces are coherent where users would reasonably expect them
   - ensure recent helper behavior is described accurately in the main documentation trail

2. lightweight consistency cleanup
   - only when needed to reduce confusion around completed work
   - examples: small naming inconsistencies, stale references, or thin integration gaps across docs/tests

### Not the focus

This is not the pass for:

- broad internal refactors
- large module reorganizations
- cleanup motivated only by aesthetic preference

## Recommended Files To Inspect

The implementation plan should likely touch or inspect:

- `docs/superpowers/specs/`
- `docs/superpowers/plans/`
- [README.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md)
- [docs/region-contrasts.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md)
- [docs/shared-clustering.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/shared-clustering.md)
- [docs/region-discovery.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-discovery.md)
- [docs/global-analysis.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/global-analysis.md)
- [dimelo/__init__.py](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/__init__.py)
- [dimelo/plotting.py](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py)

The implementation plan does not need to modify all of them. This list exists to constrain where the integration cleanup should look.

## Status / Navigation Interaction

The new central index in `docs/superpowers/README.md` should remain the status source of truth.

This cleanup may update it if:

- normalized paths require link fixes
- status wording needs correction after integration cleanup

But it should not replace the central-index model introduced in the prior cleanup slice.

## Verification

The implementation plan should include two kinds of verification.

### Docs verification

- no remaining machine-local absolute paths in tracked `docs/superpowers` historical docs, unless intentionally preserved
- no broken links introduced by normalization
- no remaining untracked plan/spec artifacts related to this cleanup

### Code/docs integration verification

- any touched docs reference real current helpers/workflows
- any touched exports or helper references resolve in the repo
- if tests are touched, run the relevant focused test subset

## Recommended Build Order

1. Normalize path references in tracked `docs/superpowers/specs` and `docs/superpowers/plans`.
2. Verify no broken relative links were introduced.
3. Inspect the public docs surface for stale or missing integration points around completed work.
4. Apply only the smallest integration cleanup needed for discoverability and consistency.
5. Re-verify docs references and any touched code/test surfaces.

