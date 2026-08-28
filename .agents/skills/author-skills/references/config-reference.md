# Config Reference -- Skill Settings

All skill-related settings live under the `skill:` key in
`great-docs.yml`. This reference covers every option.

## Top-level skill keys

```yaml
skill:
  enabled: true
  file: skills/my-package/SKILL.md
  well_known: true
  gotchas: []
  best_practices: []
  decision_table: []
  extra_body: null
  skills: []
```

### Key reference

| Key              | Type          | Default | Description                                             |
| ---------------- | ------------- | ------- | ------------------------------------------------------- |
| `enabled`        | `bool`        | `true`  | Enable or disable skill generation entirely             |
| `file`           | `str \| null` | `null`  | Path to a hand-written SKILL.md (relative to project root) |
| `well_known`     | `bool`        | `true`  | Publish skills to `.well-known/agent-skills/` for discovery |
| `gotchas`        | `list[str]`   | `[]`    | Gotcha strings appended to auto-generated skill         |
| `best_practices` | `list[str]`   | `[]`    | Best-practice strings appended to auto-generated skill  |
| `decision_table` | `list[dict]`  | `[]`    | Manual "When to use what" rows for auto-generated skill |
| `extra_body`     | `str \| null` | `null`  | Path to extra Markdown appended to auto-generated body  |
| `skills`         | `list[dict]`  | `[]`    | Multi-skill entries (overrides `file` when non-empty)   |

## Scenario configs

### Scenario 1: Automatic generation (default)

No `skill:` config needed. Great Docs generates a skill from
package metadata and API sections.

To enrich the auto-generated skill, add optional keys:

```yaml
skill:
  gotchas:
    - "Column selectors are strings, not bare identifiers."
  best_practices:
    - "Prefer method chaining over intermediate variables."
  decision_table:
    - need: "Create a table"
      use: "GT(data)"
    - need: "Format numbers"
      use: "fmt_number()"
  extra_body: "skills/extra-content.md"
```

### Scenario 2: Single hand-written skill

```yaml
skill:
  file: skills/my-package/SKILL.md
```

When `file` is set, automatic generation is skipped. The
referenced SKILL.md is copied verbatim to `<docs>/skill.md`.

### Scenario 3: Multiple named skills

```yaml
skill:
  skills:
    - name: my-package
      file: skills/my-package/SKILL.md
    - name: write-guides
      file: skills/write-guides/SKILL.md
    - name: review-code
      file: skills/review-code/SKILL.md
```

Each entry requires two keys:

| Key    | Type  | Description                                      |
| ------ | ----- | ------------------------------------------------ |
| `name` | `str` | Unique skill identifier (used in URLs and paths) |
| `file` | `str` | Path to the SKILL.md (relative to project root)  |

## Precedence rules

Great Docs evaluates skill config in this order:

1. If `skill.enabled` is `false`, no skill is generated.
2. If `skill.skills` is non-empty, multi-skill mode is used.
   The `skill.file` key is ignored.
3. If `skill.file` is set, the hand-written file is used.
4. If a curated skill exists at `skills/<package-name>/SKILL.md`,
   it is used automatically (no config needed).
5. Otherwise, a skill is auto-generated from package metadata.

## Discovery output

After a build, skills are published at these paths:

```
<docs>/
├── skill.md                                    # primary skill
├── skills.qmd                                  # rendered Skills page
└── .well-known/
    └── agent-skills/
        ├── index.json                          # discovery manifest
        ├── my-package/
        │   ├── SKILL.md
        │   └── references/
        │       └── ...
        └── write-guides/
            ├── SKILL.md
            └── references/
                └── ...
```

The `index.json` manifest lists each skill with its `name`,
`description`, and relative path to the SKILL.md. Agent tooling
like `npx skills add` reads this manifest to offer skill
installation.
