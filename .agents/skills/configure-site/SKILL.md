---
name: configure-site
description: >
  Configure a Great Docs documentation site through great-docs.yml.
  Covers theming (navbar gradients, content glow, dark mode), hero
  sections, logos, announcements, sidebar options, page tags, page
  status badges, SEO, analytics, and deployment settings. Use when
  customizing site appearance, enabling features, or tuning build
  behavior.
license: MIT
compatibility: Requires Great Docs >=0.8, Quarto CLI installed.
metadata:
  author: rich-iannone
  version: "1.0"
  tags:
    - documentation
    - configuration
    - theming
    - deployment
---

# Configure Site

Skill for customizing a Great Docs documentation site through
`great-docs.yml`. All configuration is centralized in this single
YAML file at the project root.

## Quick start

```bash
# Generate a starter config
great-docs init

# Or generate a full template with all options
great-docs config > great-docs.yml

# Edit, then rebuild
great-docs build && great-docs preview
```

## Skill directory structure

```
skills/configure-site/
├── SKILL.md
└── references/
    ├── theming-options.md
    └── feature-matrix.md
```

## When to use this skill

| Need                      | Config key                              |
| ------------------------- | --------------------------------------- |
| Add a gradient navbar     | `navbar_style: sky`                     |
| Set a solid navbar color  | `navbar_color: "#1a1a2e"`               |
| Add content area glow     | `content_style: lilac`                  |
| Add a logo                | `logo: assets/logo.svg`                 |
| Add light/dark logos      | `logo: {light: ..., dark: ...}`         |
| Set a favicon             | `favicon: assets/favicon.svg`           |
| Enable hero section       | `hero: true` or `hero: {enabled: true}` |
| Add starfield animation   | `hero: {starfield: true}`               |
| Show announcement banner  | `announcement: {content: "...", ...}`   |
| Enable dark mode toggle   | `dark_mode_toggle: true`                |
| Enable page tags          | `tags: {enabled: true}`                 |
| Enable page status badges | `page_status: true`                     |
| Configure sidebar filter  | `sidebar_filter: {enabled: true}`       |
| Set display name          | `display_name: "My Package"`            |
| Exclude symbols from docs | `exclude: [PrivateClass]`               |
| Enable CLI documentation  | `cli: {enabled: true, module: ...}`     |
| Set docstring parser      | `parser: numpy`                         |
| Use static analysis       | `dynamic: false`                        |
| Enable link checker       | Built-in: `great-docs check-links`      |
| Add analytics             | `include_in_header: [{text: ...}]`      |
| Enable multi-version docs | `versions: [{label: ..., tag: ...}]`    |

## Core concepts

### Configuration file location

`great-docs.yml` must be at the project root, alongside
`pyproject.toml`. Great Docs reads it automatically on every
command.

### Configuration categories

The config file is organized into logical groups:

1. **Package metadata** — `module`, `display_name`, `parser`,
   `dynamic`, `exclude`
2. **GitHub integration** — `repo`, `github_style`, `source`
3. **Navigation & theming** — `navbar_style`, `navbar_color`,
   `content_style`, `dark_mode_toggle`, `nav_icons`
4. **Branding** — `logo`, `favicon`, `hero`, `announcement`
5. **Content** — `sections`, `cli`, `tags`, `page_status`,
   `sidebar_filter`
6. **Build & publish** — `versions`, `skill`, `changelog`, `seo`,
   `include_in_header`
7. **Author metadata** — `authors`, `funding`

### Theming presets

Great Docs ships with named gradient presets for the navbar and
content area. See
[references/theming-options.md](references/theming-options.md) for
the full palette.

```yaml
navbar_style: sky       # gradient navbar
content_style: lilac    # content glow

# Or limit content glow to the homepage only:
content_style:
  preset: lilac
  pages: homepage
```

Available presets: `sky`, `peach`, `prism`, `lilac`, `slate`,
`honey`, `dusk`, `mint`.

### Logo configuration

```yaml
# Single logo for both themes
logo: assets/logo.svg

# Separate light/dark logos
logo:
  light: assets/logo-light.svg
  dark: assets/logo-dark.svg
```

### Hero section

The hero is the large landing area on the homepage:

```yaml
hero:
  enabled: true
  tagline: "Your tagline here."
  starfield: true # animated starfield canvas
```

Or use the shorthand: `hero: true` (uses package description as
the tagline).

### Announcement banner

```yaml
announcement:
  content: "v2 is out! <a href='...'>Read more</a>"
  type: info # info, warning, success, danger
  style: mint # gradient preset (optional)
  dismissable: true # show close button
```

Shorthand: `announcement: "Your message here"`.

### Page tags

```yaml
tags:
  enabled: true
  icons:
    Getting Started: rocket
    Configuration: cog
    API: file-code
```

Pages opt in via frontmatter: `tags: [Getting Started, API]`.

### Page status badges

```yaml
page_status: true
```

Pages set their status in frontmatter:

```yaml
---
title: "New Feature"
status: experimental # experimental, new, stable, deprecated
---
```

### Navigation icons

Add Lucide icons to navbar and sidebar labels:

```yaml
nav_icons:
  navbar:
    User Guide: book-open
    Reference: code
  sidebar:
    Installation: download
    Quick Start: rocket
```

### Multi-version documentation

```yaml
versions:
  - label: "2.0 (dev)"
    tag: dev
    version: "2.0.0"
    prerelease: true
  - label: "1.0"
    tag: "1.0"
    git_ref: v1.0.0
    latest: true
  - label: "0.9"
    tag: "0.9"
    git_ref: v0.9
    eol: true # end-of-life badge
```

### SEO

```yaml
seo:
  sitemap: true
  canonical:
    base_url: https://your-package.readthedocs.io/
  meta:
    description: "Custom meta description"
```

### Excluding symbols

```yaml
exclude:
  - _InternalClass
  - deprecated_function
  - Config # re-exported third-party type
```

## Workflows

### Setting up a new site's appearance

```
Task Progress:
- [ ] Step 1: Choose a navbar style
- [ ] Step 2: Set content glow
- [ ] Step 3: Add logo and favicon
- [ ] Step 4: Configure hero section
- [ ] Step 5: Build and preview
```

**Step 1**: Pick a `navbar_style` preset. Preview each by building.

**Step 2**: Set `content_style` to a complementary preset.

**Step 3**: Add SVG or PNG logo files to `assets/` and configure
the `logo` key. Set `favicon` similarly.

**Step 4**: Enable `hero` with a tagline. Add `starfield: true`
for the animated background.

**Step 5**: Run `great-docs build && great-docs preview`.

### Enabling a feature

1. Find the relevant config key (use the table above or
   [references/feature-matrix.md](references/feature-matrix.md)).
2. Add it to `great-docs.yml`.
3. Rebuild with `great-docs build`.
4. Check the rendered output in the browser.

### Migrating from a minimal config

If you started with `great-docs init` and want to add more features:

1. Run `great-docs config` to see all available options.
2. Copy the sections you want into your existing `great-docs.yml`.
3. Customize values.
4. Rebuild.

## Gotchas

1. **YAML indentation matters.** Use 2-space indentation. Tabs
   cause parse errors.
2. **`module` is the import name.** For PyPI package `py-shiny`,
   set `module: shiny`, not `module: py-shiny`.
3. **Preset names are case-sensitive.** Use lowercase: `sky`, not
   `Sky`.
4. **`hero: true` is shorthand.** For full control, use the dict
   form with `enabled`, `tagline`, and `starfield` keys.
5. **Logo paths are relative to the project root.** Not relative
   to `great-docs.yml`.
6. **`content_style` on homepage only.** Use the dict form with
   `pages: homepage` to avoid the glow on every page.
7. **Changes require rebuild.** Config changes are not hot-reloaded.
   Always run `great-docs build` after editing `great-docs.yml`.
8. **Version `git_ref` must exist.** The tag or branch must exist
   in the Git repository or the versioned build will fail.
