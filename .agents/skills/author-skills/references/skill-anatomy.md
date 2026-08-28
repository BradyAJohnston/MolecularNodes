# Skill Anatomy

A SKILL.md file has two parts: YAML frontmatter and a Markdown
body. This reference covers every supported field.

## Frontmatter

```yaml
---
name: my-package
description: >
  One-paragraph description of the skill. Maximum 1024 characters.
  Should explain what the skill covers and when to use it.
license: MIT
compatibility: Requires Python >=3.11.
metadata:
  author: github-handle
  version: "1.0"
  tags:
    - python
    - relevant-topic
---
```

### Frontmatter keys

| Key             | Required | Type     | Max length | Description                                        |
| --------------- | -------- | -------- | ---------- | -------------------------------------------------- |
| `name`          | yes      | `str`    | 64 chars   | Skill identifier; lowercase with hyphens (e.g., `my-skill`) |
| `description`   | yes      | `str`    | 1024 chars | What the skill covers and when to use it            |
| `license`       | no       | `str`    |            | SPDX license identifier (e.g., `MIT`, `Apache-2.0`) |
| `compatibility` | no       | `str`    |            | Runtime requirements (Python version, CLI tools)    |
| `metadata`      | no       | `object` |            | Container for `author`, `version`, and `tags`       |

### Metadata sub-keys

| Key       | Type     | Description                                    |
| --------- | -------- | ---------------------------------------------- |
| `author`  | `str`    | GitHub handle or name of the skill author      |
| `version` | `str`    | Semantic version of the skill content          |
| `tags`    | `list`   | Searchable tags for skill discovery            |

## Body sections

The body is standard Markdown. The following sections are
conventional for Agent Skills, though none are strictly required.

### Recommended section order

| Section             | Heading level | Purpose                                          |
| ------------------- | ------------- | ------------------------------------------------ |
| Package title       | `#`           | Package name as the top-level heading            |
| Opening paragraph   |               | 2-3 sentences: what, when, why                   |
| Installation        | `##`          | `pip install` command                            |
| When to use what    | `##`          | Decision table mapping tasks to API calls        |
| API overview        | `##`          | Section-by-section API listing with summaries    |
| Gotchas             | `##`          | Numbered list of common pitfalls                 |
| Best practices      | `##`          | Bullet list of recommended patterns              |
| Resources           | `##`          | Links to docs, llms.txt, source code             |

### Decision table format

The "When to use what" table maps user needs to specific API
calls. Each row should be actionable:

```markdown
| Need                    | Use                  |
| ----------------------- | -------------------- |
| Create a table from CSV | `GT(pd.read_csv())` |
| Format as percentage    | `fmt_percent()`      |
| Add a footnote          | `tab_footnote()`     |
```

### API overview format

List functions with one-line summaries. Group by API reference
section:

```markdown
## API overview

### Table creation

- `GT(data)`: Create a display table from a DataFrame.
- `GT.as_raw_html()`: Render to an HTML string.

### Formatting

- `fmt_number()`: Format numeric columns.
- `fmt_currency()`: Format as currency values.
```

### Gotchas format

Use a numbered list. Lead each item with a bold label:

```markdown
## Gotchas

1. **Immutable objects.** Methods return new GT objects.
2. **Column selectors are strings.** Pass column names as
   strings, not bare identifiers.
```

## Reference files

Companion files live in a `references/` subdirectory alongside
the SKILL.md:

```
skills/my-package/
├── SKILL.md
└── references/
    ├── style-guide.md
    ├── pattern-library.md
    └── checklist.md
```

Reference files are plain Markdown (no frontmatter required).
They are copied into `.well-known/agent-skills/<name>/references/`
during the build, making them available to agents that install
the skill.

Use reference files for detailed content that would make the
main SKILL.md too long: style examples, checklists, decision
matrices, and extended pattern libraries.
