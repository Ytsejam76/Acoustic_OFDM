# AGENTS.md

Repository-specific instructions for Codex and related coding agents.

## General

- Preserve separation of concerns aggressively.
- Keep the CLI thin. The CLI is a tool to run and test the modem, not the place for modem logic.
- Keep modem subsystems in separate modules when they represent distinct concerns.
- Wake-up logic and coarse search/synchronization belong in the modem library, but in separate modules.
- Equalization belongs in its own module.
- Debug/dump helpers used for inspection should live in dedicated modules, not mixed into core DSP logic.

## Modem architecture

- Prefer explicit structure over heuristic entanglement.
- Keep baseband DSP, passband conversion, synchronization, equalization, wake/search, and diagnostics separable.
- When extending the equalizer, prefer composable features/configuration over growing mode enums.
- Favor builders and explicit feature flags for equalizer configuration.

## Documentation

- Public types and public configuration enums must be documented properly.
- Add field-level documentation to public structs when the field meaning is not trivial.
- Add rationale comments before important methods, especially in DSP/equalization/sync code.
- Equations in comments are welcome and should be written in a documentation-friendly math style.
- Keep crate/module naming precise. Avoid names that are misleadingly narrow or ambiguous.

## Rust style

- Always run `cargo fmt` after Rust code edits. Do not ask whether to run it; just run it.
- Keep naming clear and unambiguous. Avoid names that collide conceptually with standard Rust traits when a clearer name exists.
- Prefer small, explicit helpers over opaque monolithic functions.

## Git and history

- Prefer `git mv` for renames instead of delete-and-recreate renames.
- Commit messages must use this format:
  - `<main_module>: <summary>`
  - a blank line
  - a short, idiomatic body explaining what changed and why
- Every line in the commit message must be at most 80 columns.
- Keep the subject line at or below 80 columns.
- Wrap the body at 80 columns; rephrase instead of letting lines exceed the limit.
- When asked to keep a commit focused, do not mix unrelated changes into it.

## README

- Keep the README aligned with the actual working scripts and current recommended workflow.
- Prefer shorter, current quickstart sections over long stale descriptions.
- If a directory contains a `README.md`, keep it up to date with the current
  development in that part of the codebase.
