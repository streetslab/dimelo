# Superpowers Docs Cleanup Design

## Purpose

This design defines the first cleanup pass for the internal planning trail under `docs/superpowers`.

The goal is to make the accumulated specs and plans:

- tracked in git
- easy to navigate
- clear about what has been implemented versus what remains follow-on work

This is a docs-first cleanup slice. It does not yet change package code or user-facing APIs.

## Problem

The repository now contains a substantial internal design and planning trail:

- specs under `docs/superpowers/specs`
- plans under `docs/superpowers/plans`

Some of the plan files exist on disk but are still untracked. The docs set is also missing a simple navigation layer, so it is difficult to tell:

- which specs belong to which plans
- which plans were executed already
- which slices are still only planned
- which later plans supersede or extend earlier ones

The cleanup should preserve this history without requiring every historical file to be rewritten.

## Core Design Rule

Status should be tracked centrally, not embedded into every spec and plan file.

This avoids:

- rewriting historical design docs
- duplicating status information across many files
- future drift between per-file annotations

The first cleanup slice should therefore introduce one central index file instead of editing every existing file header.

## Scope

This cleanup slice should do three things:

1. track the currently untracked plan files in git
2. add a central navigation/index document for `docs/superpowers`
3. map specs and plans to current implementation status in that index

It should not:

- rewrite the historical spec content
- merge or delete old spec/plan files
- introduce a heavy metadata system
- touch code outside documentation cleanup

## Recommended File Structure

Add one central index file:

- [README.md](../README.md)

Keep the existing layout:

- `docs/superpowers/specs/`
- `docs/superpowers/plans/`

Do not reorganize those directories in this slice.

## Index Content

The new central index should have four sections:

### 1. Overview

Briefly explain what `docs/superpowers` contains:

- specs are design documents
- plans are executable implementation breakdowns
- statuses in the index reflect the current branch/repo state

### 2. Specs

List the existing specs grouped by theme, with a short description and one of:

- `implemented`
- `partially implemented`
- `planned`

### 3. Plans

List the implementation plans grouped by theme, with a short description and one of:

- `implemented`
- `partially implemented`
- `not started`

### 4. Current Themes

Group related work so the trail is readable, for example:

- shared clustering
- region contrasts
- region discovery
- global analysis
- plotting architecture

This section can point users to the most relevant spec/plan sequence for each theme.

## Status Semantics

The cleanup should use lightweight status labels only.

Recommended meanings:

- `implemented`
  The planned slice has been completed and corresponding code/docs are in the repo.

- `partially implemented`
  Some downstream or follow-on pieces exist, but the spec/plan family is not fully exhausted.

- `planned`
  The design exists, but implementation for that slice has not yet landed.

- `not started`
  A plan exists but has not been executed.

These labels should be applied only in the central index.

## Treatment of Existing Untracked Plans

The current untracked files under `docs/superpowers/plans` should be committed rather than removed.

Rationale:

- they correspond to real design/implementation work
- they are useful historical context
- they make the executed development trail auditable

This cleanup slice should not try to prune them aggressively.

## Minimal Editing Principle

This slice should prefer:

- adding one new index
- committing the existing plan files as-is

over:

- editing every plan/spec header
- retrofitting status fields into older docs
- renaming the existing files

That keeps the cleanup low-risk and easy to maintain.

## Implementation Notes

The index should reflect the current codebase state on the branch where the cleanup is executed.

That means it should call out, at minimum, the slices that have clearly landed:

- shared clustering foundations
- region contrasts foundations
- global analysis foundations
- region discovery foundations
- discovery to clustering workflow
- discovery to clustering to contrasts workflow
- cluster occupancy region contrasts
- plotting axis architecture
- region contrasts plotting

If a broader feature area still has follow-on work remaining, mark the family as `partially implemented` rather than overstating completion.

## Tests / Verification

This is docs-only cleanup, so verification can stay lightweight:

- confirm the new index references the actual existing files
- confirm the untracked plan files are now staged/committed
- confirm there are no broken file paths in the index text

## Recommended Build Order

1. Add `docs/superpowers/README.md`.
2. Populate it with grouped spec/plan links and lightweight status labels.
3. Add the currently untracked plan docs to git.
4. Verify the index paths match the actual repo contents.
5. Commit the docs cleanup slice.

