---
name: write-user-guide
description: >
  Write and maintain narrative user-guide pages for a Great Docs site.
  Covers page creation, QMD frontmatter, section grouping, sidebar
  ordering, callouts, executable code cells, cross-references, and
  content guidelines. Use when adding, reorganizing, or improving
  user-guide content.
license: MIT
compatibility: Requires Great Docs >=0.8, Quarto CLI installed.
metadata:
  author: rich-iannone
  version: "1.0"
  tags:
    - documentation
    - user-guide
    - quarto
    - content-authoring
---

# Write User Guide

Skill for authoring user-guide pages in a Great Docs documentation
site. User guides provide narrative documentation (tutorials,
conceptual explanations, and task walkthroughs) that complement
the auto-generated API reference.

## Quick start

```bash
mkdir -p user_guide
cat > user_guide/00-introduction.qmd << 'EOF'
---
title: "Introduction"
guide-section: "Getting Started"
tags: [Getting Started]
---

# Introduction

Welcome to the project. This guide walks you through...
EOF

great-docs build
```

## Skill directory structure

```
skills/write-user-guide/
├── SKILL.md
└── references/
    ├── page-anatomy.md
    └── writing-guidelines.md
```

## When to use this skill

| Need                         | Action                                         |
| ---------------------------- | ---------------------------------------------- |
| Add a new guide page         | Create `user_guide/NN-topic.qmd`               |
| Reorder pages                | Rename numeric prefixes                        |
| Group pages into sections    | Set `guide-section` in frontmatter             |
| Add an interactive example   | Use `{python}` code cells in the `.qmd`        |
| Cross-reference another page | Use `[text](../user-guide/page.qmd)` links     |
| Embed a callout              | Use `:::{.callout-tip}` / `:::{.callout-note}` |
| Add images                   | Place in `assets/` and reference from QMD      |

## Core concepts

### File naming convention

Every page in `user_guide/` must have a two-digit numeric prefix
that controls sidebar ordering:

```
user_guide/
├── 00-introduction.qmd     # appears first
├── 01-installation.qmd
├── 02-quickstart.qmd
├── 03-authoring-qmd-files.qmd
└── 04-writing-docstrings.qmd
```

Great Docs strips the prefix for clean URLs:
`00-introduction.qmd` → `user-guide/introduction.html`.

### QMD frontmatter

Every user-guide page starts with YAML frontmatter:

```yaml
---
title: "Writing Docstrings"
guide-section: "Getting Started"
tags: [API, Content]
---
```

**Required keys:**

| Key     | Description                    |
| ------- | ------------------------------ |
| `title` | Page heading and sidebar label |

**Optional keys:**

| Key             | Description                                        |
| --------------- | -------------------------------------------------- |
| `guide-section` | Group pages under a sidebar section header         |
| `tags`          | Content tags for discoverability                   |
| `bread-crumbs`  | Set `false` to hide breadcrumb navigation          |
| `status`        | Page status badge: `experimental`, `new`, `stable` |

### Guide sections

Pages with the same `guide-section` value are grouped together in
the sidebar under a collapsible section heading:

```
Getting Started
├── Introduction
├── Installation
└── Quick Start
Site Content
├── Authoring QMD Files
├── Writing Docstrings
└── User Guides
```

If no `guide-section` is set, the page appears at the top level.

### Page body structure

A well-structured page follows this outline:

```markdown
# Page Title

Opening paragraph: 2-3 sentences explaining what this page covers
and why the reader cares.

## First Major Section

Narrative prose. Keep paragraphs short (3-5 sentences).

### Subsection

More detail. Use tables, code blocks, and callouts to break up text.

## Second Major Section

...
```

**Guidelines:**

- Start every page with a single `#` heading matching the `title`.
- Use `##` for major sections, `###` for subsections.
- Keep the hierarchy flat; avoid `####` if possible.
- Lead each section with a sentence explaining what follows.
- End with a summary or "next steps" when appropriate.

### Callouts

Quarto callouts highlight important information:

```markdown
:::{.callout-tip}

## Pro tip

You can combine `guide-section` with `tags` to make pages
discoverable from multiple angles.
:::

:::{.callout-warning}

## Watch out

Renaming a page file changes its URL. Update any cross-references.
:::

:::{.callout-note}
This feature requires Great Docs 0.8 or later.
:::
```

Available types: `note`, `tip`, `warning`, `caution`, `important`.

### Executable code cells

Embed live Python examples that run during the build:

````markdown
```{python}
import great_docs
gd = great_docs.GreatDocs()
print(gd.project_path)
```
````

**`{python}` vs `{.python}`**: this distinction is critical.

- `{python}` (no dot) creates an **executable** code cell. Quarto
  runs it through the Jupyter kernel during the build and captures
  the output.
- `{.python}` (with a dot) creates a **display-only** code block.
  Quarto syntax-highlights it but never executes it.

Use `{python}` when the output matters (tables, plots, printed
values). Use `{.python}` for illustrative snippets where execution
is unnecessary or undesirable.

Additional cell-level controls:

- Use `#| eval: false` to show code without executing it.
- Use `#| echo: false` to show only the output.

### Table previews and explorers

When a page involves sample datasets or transformed DataFrames,
use the built-in table widgets instead of raw `print()` output.

**Shortcodes** (for static data files in `assets/data/`):

```markdown
{{< tbl-preview file="assets/data/students.csv" >}}
{{< tbl-explorer file="assets/data/students.csv" >}}
```

**Python API** (for DataFrames produced in executable cells):

````markdown
```{python}
from great_docs import tbl_preview, tbl_explorer

tbl_preview(df)           # compact head/tail preview
tbl_explorer(df)          # interactive: sort, filter, paginate
```
````

Use `tbl-preview` (or `tbl_preview()`) for a quick glance at a
dataset. Use `tbl-explorer` (or `tbl_explorer()`) when readers
need to sort, search, or paginate the data. Both accept Pandas
DataFrames, Polars DataFrames, and file paths (CSV, TSV, Parquet,
Arrow, JSONL).

### Cross-references

Link to other pages in the site:

```markdown
See the [Configuration](../user-guide/configuration.qmd) page.
See the [API reference](../reference/GreatDocs.qmd) page.
```

Use relative paths from the rendered output location
(`great-docs/user-guide/`), not the source.

### Images and assets

Place images in `assets/` at the project root:

```markdown
![Architecture diagram](../assets/architecture.png)
```

Great Docs copies the `assets/` directory into the build
automatically.

## Workflows

### Adding a new page

```
Task Progress:
- [ ] Step 1: Choose a filename
- [ ] Step 2: Write frontmatter
- [ ] Step 3: Write content
- [ ] Step 4: Build and preview
```

**Step 1**: Pick the next numeric prefix. If the last file is
`08-user-guides.qmd`, name your file `09-new-topic.qmd`.

**Step 2**: Add frontmatter with `title`, `guide-section`, and
optional `tags`.

**Step 3**: Write the body using the page structure guidelines
above.

**Step 4**: Run `great-docs build && great-docs preview` and check
the sidebar ordering and rendered content.

### Reorganizing existing pages

1. Rename files to adjust numeric prefixes.
2. Update `guide-section` values to regroup.
3. Rebuild. Great Docs regenerates the sidebar automatically.
4. Check for broken cross-references.

### Converting a README into a guide page

1. Copy the README content into a new `.qmd` file.
2. Add frontmatter with `title` and `guide-section`.
3. Replace any GitHub-flavored Markdown extensions with Quarto
   equivalents (e.g., `> [!NOTE]` → `:::{.callout-note}`).
4. Rebuild and verify.

## Gotchas

1. **Numeric prefixes control ordering by default.** Without them,
   pages sort alphabetically. Alternatively, you can omit prefixes
   and define an explicit page order in `great-docs.yml`.
2. **Don't skip numbers.** Gaps are fine (`01`, `03`, `05`) but
   large jumps make it hard to insert pages later.
3. **Title must match the `#` heading.** If `title: "Foo"` but
   the body starts with `# Bar`, the sidebar says "Foo" but the
   page says "Bar".
4. **`guide-section` is case-sensitive.** `"Getting Started"` and
   `"getting started"` create separate sections.
5. **Don't nest directories.** All pages must be directly in
   `user_guide/`, not in subdirectories.
6. **Hyphens, not underscores, in filenames.** Great Docs converts
   underscores to hyphens in URLs, so `my_page.qmd` becomes
   `my-page.html`.
7. **Lists need a blank line before them.** A bullet or numbered
   list that immediately follows a paragraph (no blank line) will
   not be parsed as a list by Quarto; it renders as plain text.
