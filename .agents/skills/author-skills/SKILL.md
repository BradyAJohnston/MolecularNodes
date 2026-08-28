---
name: author-skills
description: >
  Author, configure, and distribute Agent Skills for a Great Docs
  site. Covers the three scenarios: automatic skill generation,
  adding a single hand-written skill, and distributing multiple
  named skills with a switcher page. Use when creating, editing,
  or configuring SKILL.md files for package documentation.
license: MIT
compatibility: Requires Great Docs >=0.8, Quarto CLI installed.
metadata:
  author: rich-iannone
  version: "1.0"
  tags:
    - skills
    - agent-skills
    - distribution
    - configuration
---

# Author Skills

Skill for creating and distributing Agent Skills through a Great
Docs documentation site. An Agent Skill is a structured Markdown
file (`SKILL.md`) that gives AI coding agents context about a
package, its API, and best practices.

Great Docs supports three scenarios for skills, each with a
different level of effort and control.

## The three scenarios

### Scenario 1: Automatic generation

Great Docs generates a `skill.md` automatically from your
package metadata, API reference sections, and `great-docs.yml`
config. This is the zero-effort default.

**How it works**: during `great-docs build`, the tool reads your
`pyproject.toml` (for name, description, license, Python version),
your API reference sections, and any `skill.*` config keys. It
produces a skill file with installation instructions, an API
overview, gotchas, best practices, and resource links.

**Config options that enrich the auto-generated skill:**

```yaml
# great-docs.yml
skill:
  enabled: true # default; set false to skip
  well_known: true # publish to .well-known/agent-skills/
  gotchas:
    - "Column expressions are lazy; call .collect() to materialize."
    - "Method chaining returns a new object; the original is not mutated."
  best_practices:
    - "Prefer method chaining over intermediate variables."
  decision_table:
    - need: "Create a table from a dict"
      use: "GT(data)"
    - need: "Format currency values"
      use: "fmt_currency()"
  extra_body: "skills/extra-skill-content.md"
```

**When to use this scenario**: you want a skill page with minimal
effort and are satisfied with auto-generated content derived from
your API reference.

### Scenario 2: Single hand-written skill

Write a `SKILL.md` by hand (or with LLM assistance) and point
`great-docs.yml` at it. Great Docs copies it into the build and
generates the Skills page from it.

**Directory layout:**

```
project-root/
├── skills/
│   └── my-package/
│       ├── SKILL.md
│       └── references/
│           ├── api-patterns.md
│           └── gotchas.md
└── great-docs.yml
```

**Config:**

```yaml
# great-docs.yml
skill:
  file: skills/my-package/SKILL.md
```

Great Docs copies the file to `<docs>/skill.md`, places it
under `.well-known/agent-skills/<name>/SKILL.md`, and generates
the Skills page with install instructions and the full skill
content.

**When to use this scenario**: you want full editorial control
over the skill content, including custom sections, curated
examples, and reference files that go beyond what automatic
generation can produce.

### Scenario 3: Multiple named skills

Distribute several focused skills from a single site. Each skill
gets its own panel on the Skills page (with a switcher bar) and
its own entry in `.well-known/agent-skills/index.json` for
`npx skills add` discovery.

**Directory layout:**

```
project-root/
├── skills/
│   ├── my-package/
│   │   ├── SKILL.md
│   │   └── references/
│   ├── write-guides/
│   │   ├── SKILL.md
│   │   └── references/
│   └── review-code/
│       ├── SKILL.md
│       └── references/
└── great-docs.yml
```

**Config:**

```yaml
# great-docs.yml
skill:
  skills:
    - name: my-package
      file: skills/my-package/SKILL.md
    - name: write-guides
      file: skills/write-guides/SKILL.md
    - name: review-code
      file: skills/review-code/SKILL.md
```

The first entry becomes the primary `skill.md` at the site root.
All entries appear in the switcher bar and are individually
installable via `npx skills add`.

**When to use this scenario**: your project has distinct
focus-area skills (e.g., one for general usage, one for
contributing, one for configuration) and you want users to
install exactly what they need.

## Quick start

Create a hand-written skill in three steps:

````bash
mkdir -p skills/my-package/references

cat > skills/my-package/SKILL.md << 'SKILL_EOF'
---
name: my-package
description: >
  Use the my-package Python library. Covers installation, core API,
  and common patterns.
license: MIT
compatibility: Requires Python >=3.11.
metadata:
  author: your-name
  version: "1.0"
  tags:
    - python
    - my-package
---

# my-package

Short description of what the package does and when to use it.

## Installation

```bash
pip install my-package
````

## When to use what

| Need            | Use              |
| --------------- | ---------------- |
| Create a widget | `Widget()`       |
| Style a widget  | `widget.style()` |

## API overview

### Core

- `Widget()`: Create a new widget.
- `Widget.style()`: Apply styling to a widget.

## Gotchas

1. **Immutable by default.** Methods return new objects.
2. **Lazy evaluation.** Call `.render()` to produce output.

## Best practices

- Prefer method chaining over intermediate variables.
- Use type hints in all public API calls.
  SKILL_EOF

````

Then add one line to `great-docs.yml`:

```yaml
skill:
  file: skills/my-package/SKILL.md
````

Run `great-docs build` and verify the Skills page.

## Skill directory structure

```
skills/author-skills/
├── SKILL.md
└── references/
    ├── skill-anatomy.md
    └── config-reference.md
```

## Writing a SKILL.md

### Frontmatter

Every SKILL.md starts with YAML frontmatter:

```yaml
---
name: my-package
description: >
  One-paragraph description of the skill (max 1024 characters).
license: MIT
compatibility: Requires Python >=3.11.
metadata:
  author: your-github-handle
  version: "1.0"
  tags:
    - python
    - relevant-topic
---
```

The required and optional keys are documented in the
`references/skill-anatomy.md` companion file.

### Body structure

A well-structured skill follows this outline:

```
# Package Name
Opening paragraph (what, when, why).

## Installation
pip install command.

## When to use what
Decision table mapping needs to API calls.

## API overview
Sections matching your API reference, with one-line summaries.

## Gotchas
Numbered list of common pitfalls.

## Best practices
Bullet list of recommended patterns.

## Resources
Links to docs, llms.txt, and source code.
```

### Reference files

Place companion Markdown files in a `references/` subdirectory
alongside the SKILL.md. These files are copied into
`.well-known/agent-skills/<name>/references/` and are available
to agents that install the skill.

Use reference files for content that is too detailed for the
main SKILL.md: style guides, checklists, pattern libraries,
and decision matrices.

### Writing tips

Follow these guidelines when writing skill content:

- Write for an LLM reader, not a human. Be explicit about
  distinctions that a model might confuse (e.g., `{python}` vs
  `{.python}`).
- Use concrete examples over abstract descriptions.
- Keep the decision table ("When to use what") actionable:
  each row should map a task to a specific function or method.
- Avoid emoji, em dashes, and starting a list without
  introductory text.
- Keep the SKILL.md under 2000 lines. Move detailed reference
  material into companion files.

## Workflows

### Adding a skill to an existing site

1. Create the `skills/<name>/` directory with a `SKILL.md`.
2. Optionally add a `references/` subdirectory with companion
   files.
3. Set `skill.file` in `great-docs.yml` to point at the SKILL.md.
4. Run `great-docs build` and check the Skills page.

### Converting from automatic to hand-written

1. Run `great-docs build` to generate the automatic `skill.md`.
2. Copy the generated file from `<docs>/skill.md` into
   `skills/<name>/SKILL.md`.
3. Edit the content: add custom sections, rewrite descriptions,
   curate the decision table.
4. Set `skill.file` in `great-docs.yml`.
5. Rebuild.

### Adding a second skill (switching to multi-skill mode)

1. Create a second `skills/<name>/` directory with its SKILL.md.
2. Replace the `skill.file` key with `skill.skills` in
   `great-docs.yml`:

   ```yaml
   skill:
     skills:
       - name: original-skill
         file: skills/original-skill/SKILL.md
       - name: new-skill
         file: skills/new-skill/SKILL.md
   ```

3. Rebuild. The Skills page now shows a switcher bar with both
   skills.

### Testing skill installation

Verify your skill is discoverable:

```bash
# Check the published skill
great-docs check-skill

# List all skills in the .well-known manifest
great-docs list-skills

# Test npx installation (after deploying)
npx skills add <site-url>
```

## Discovery and distribution

Great Docs publishes skills at two well-known paths for
auto-discovery by `npx skills add` and other agent tooling:

```
<site>/.well-known/agent-skills/<name>/SKILL.md
<site>/.well-known/agent-skills/index.json
```

The `index.json` manifest lists all skills with their name,
description, and relative path. In multi-skill mode, each skill
gets its own entry.

The root-level `<site>/skill.md` always contains the primary
skill (the first entry in `skill.skills`, or the only skill).

## Gotchas

1. **`skill.skills` overrides `skill.file`.** If both are set,
   `skill.skills` takes precedence and `skill.file` is ignored.
2. **The `name` field must be unique.** Each skill in
   `skill.skills` needs a distinct `name`. Duplicates cause the
   later entry to overwrite the earlier one in `.well-known/`.
3. **Reference files must be in `references/`.** Only files under
   the `references/` subdirectory (relative to the SKILL.md) are
   copied into `.well-known/`. Files elsewhere are ignored.
4. **Use hyphens in skill names, not underscores.** The `name`
   field becomes a URL path segment
   (`.well-known/agent-skills/<name>/`), and hyphens are the
   standard web convention for URL slugs. Match the directory
   name: `skills/my-skill/` with `name: my-skill`.
5. **The first skill in `skill.skills` is special.** It becomes
   the primary `skill.md` at the site root and is the default
   panel shown on the Skills page.
6. **Set `skill.enabled: false` to disable entirely.** This
   suppresses both automatic generation and hand-written skill
   processing.
