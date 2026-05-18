# rajiveplus Context

This context defines the language used to describe RaJIVE+ multi-view decomposition work. It is a glossary for the domain concepts, not an implementation specification.

## Language

**Data Block**:
One measured data layer in a multi-view study.
_Avoid_: view, assay, matrix

**Joint Component**:
The shared structure estimated across multiple **Data Blocks**.
_Avoid_: shared part

**Individual Component**:
The structure estimated as specific to one **Data Block**.
_Avoid_: block-specific part, Indiv

**Residual Component**:
The remaining structure in a **Data Block** after the **Joint Component** and **Individual Component**.
_Avoid_: noise, Resid

**Observed Mask**:
A logical indicator of which entries in each **Data Block** are observed and available for native missing-data fitting.
_Avoid_: missing mask

## Relationships

- A RaJIVE+ analysis compares two or more **Data Blocks**.
- A **Data Block** can have observations absent for some study samples while other **Data Blocks** still contain those samples.
- Each **Data Block** is decomposed into a **Joint Component**, an **Individual Component**, and a **Residual Component**.
- An **Observed Mask** belongs to exactly one **Data Block**.

## Example Dialogue

> **Dev:** "Should we call mRNA and proteomics views in the refactoring docs?"
> **Domain expert:** "No. Call each one a **Data Block**; 'view' can appear informally, but **Data Block** is the canonical term."
>
> **Dev:** "Is the third part of a decomposition noise?"
> **Domain expert:** "Call it the **Residual Component**. It is what remains after joint and individual structure, not necessarily pure random noise."
>
> **Dev:** "Can public tables abbreviate components as Indiv and Resid?"
> **Domain expert:** "No. Use the full names **Individual** and **Residual** consistently in code-facing and documentation-facing labels."
>
> **Dev:** "In native missing mode, does `TRUE` in the mask mean missing?"
> **Domain expert:** "No. The mask is an **Observed Mask**: `TRUE` means the entry is observed and used for fitting."

## Flagged Ambiguities

- "view", "assay", and "matrix" were used informally for the same concept; resolved: the canonical term is **Data Block**.
- "noise" and "residual" were both used for the third decomposition role; resolved: the canonical term is **Residual Component**.
- "Indiv" and "Resid" were used as public shorthand labels; resolved: public code-facing and documentation-facing labels use **Individual** and **Residual**.
- "mask" was ambiguous because some codebases use `TRUE` to mean missing; resolved: RaJIVE+ native missing-data masks are **Observed Masks** where `TRUE` means observed.
