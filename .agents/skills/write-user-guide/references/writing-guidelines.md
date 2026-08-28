# Writing Guidelines — User Guide Content

## Voice and tone

- **Second person** ("you") for instructions.
- **Active voice** preferred over passive.
- **Present tense** for descriptions, imperative for instructions.
- Keep sentences short (under 25 words when possible).
- One idea per paragraph.
- **No emoji.** Never use emoji in generated prose.
- **No em dashes.** Use commas, semicolons, colons, or
  parentheses instead. Restructure the sentence if needed.
- **Always introduce a list.** Every bullet or numbered list
  must be preceded by a sentence or phrase that introduces it.

## Structure

### Opening paragraph

Every page starts with 2-3 sentences that answer:

1. What does this page cover?
2. Why should the reader care?

Do not start with "This page explains..." — jump straight into the
value proposition.

**Good:**

> Your docstrings are the single biggest input to your API reference.
> A well-written docstring becomes a polished reference page with
> almost no extra effort.

**Avoid:**

> This page will explain how to write docstrings for your package.

### Sections

- Use `##` for major topics, `###` for details.
- Keep sections focused — one concept per `##`.
- Lead each section with a sentence summarizing what follows.
- End long pages with a "Next steps" or "Summary" section.

### Code examples

- Show the minimal example that demonstrates the concept.
- Use realistic values (not `foo`, `bar`, `baz`).
- Prefer executable cells (`{python}`) when the output helps.
- Use `{.yaml filename="great-docs.yml"}` for config examples.
- Include comments only when the code is not self-explanatory.

### Tables

Use tables for reference-style information (options, parameters,
comparison matrices). Keep columns to 3-4 maximum.

### Data tables and DataFrames

When a page shows sample data or the output of a DataFrame
transformation, prefer the built-in table widgets over raw
`print()` output:

- `tbl_preview(df)` or `{{< tbl-preview file="..." >}}` for a
  compact head/tail view.
- `tbl_explorer(df)` or `{{< tbl-explorer file="..." >}}` for an
  interactive table with sorting, filtering, and pagination.

Both accept Pandas DataFrames, Polars DataFrames, and file paths
(CSV, TSV, Parquet, Arrow, JSONL). Use `tbl-preview` by default;
reach for `tbl-explorer` when readers benefit from interacting
with the data.

### Callouts

| Type          | Use for                                      |
| ------------- | -------------------------------------------- |
| `.callout-tip`      | Best practices, shortcuts, pro tips    |
| `.callout-note`     | Additional context, version notes      |
| `.callout-warning`  | Common mistakes, breaking changes      |
| `.callout-important`| Critical information, must-reads       |
| `.callout-caution`  | Destructive actions, data loss risks   |

Use callouts sparingly — more than 2-3 per page dilutes their impact.

## Common mistakes

1. **Walls of text.** Break prose with code blocks, tables, or
   callouts every 3-4 paragraphs.
2. **Missing frontmatter.** Every `.qmd` needs at least `title`.
3. **Deep heading nesting.** If you reach `####`, the page probably
   needs splitting into two pages.
4. **Stale cross-references.** After renaming a file, search the
   entire `user_guide/` directory for old references.
5. **Overly long pages.** If a page exceeds ~800 words, consider
   splitting it. Each page should cover one focused topic.
6. **Orphaned lists.** A list that appears without any introductory
   text feels abrupt. Always write a lead-in sentence or phrase
   before the first bullet or number.
7. **Em dashes in prose.** Quarto renders them fine, but the
   project style avoids them. Use commas, semicolons, colons, or
   parentheses instead.
