# Page Anatomy — User Guide Pages

A Great Docs user-guide page is a Quarto Markdown (`.qmd`) file with
three layers: frontmatter, body content, and optional executable code.

## Frontmatter block

```yaml
---
title: "Page Title"
guide-section: "Section Name"
tags: [Tag1, Tag2]
bread-crumbs: true
status: new
---
```

### Frontmatter keys

| Key                | Required | Type     | Default        | Description                                         |
| ------------------ | -------- | -------- | -------------- | --------------------------------------------------- |
| `title`            | yes      | `str`    |                | Page title shown in sidebar and heading             |
| `guide-section`    | no       | `str`    |                | Groups pages under a sidebar section                |
| `tags`             | no       | `list`   | `[]`           | Content tags for the page-tags widget               |
| `tag-location`     | no       | `str`    | (site default) | Override tag badge position: `top` or `bottom`      |
| `status`           | no       | `str`    |                | Status badge: `experimental`, `new`, `stable`       |
| `bread-crumbs`     | no       | `bool`   | `true`         | Show or hide breadcrumb navigation                  |
| `description`      | no       | `str`    |                | Page description for metadata and social cards      |
| `image`            | no       | `str`    |                | Social-card / Open Graph image path                 |
| `toc`              | no       | `bool`   | `true`         | Show or hide the table of contents                  |
| `toc-depth`        | no       | `int`    | `2`            | Heading depth for TOC entries                       |
| `page-navigation`  | no       | `bool`   | `true`         | Show or hide prev/next navigation links             |
| `freeze`           | no       | `bool`   | `false`        | Cache code cell output across builds                |
| `layout`           | no       | `str`    | `passthrough`  | Page layout mode: `passthrough` or `raw`            |
| `body-classes`     | no       | `str`    |                | Extra CSS classes added to the page body             |
| `code-fold`        | no       | `bool`   | `false`        | Collapse code blocks by default (page-wide)         |
| `code-summary`     | no       | `str`    | `"Code"`       | Label for collapsed code blocks (page-wide)         |

## Body content

The body uses standard Quarto Markdown:

```
# Heading 1          — page title (should match `title`)
## Heading 2         — major section
### Heading 3        — subsection
#### Heading 4       — rarely needed; avoid if possible

**Bold**, *italic*, `inline code`

- Bullet list
- Another item

1. Numbered list
2. Another item

| Col A | Col B |
| ----- | ----- |
| val   | val   |

> Blockquote

:::{.callout-note}
Callout body
:::

[Link text](url)
![Alt text](image-path)
```

## Executable code cells

~~~
```{.python}
#| label: example
#| echo: true
#| eval: true
import mypackage
mypackage.hello()
```
~~~

### Cell options

| Option         | Default | Description                                      |
| -------------- | ------- | ------------------------------------------------ |
| `label`        |         | Unique identifier for the cell (for references)  |
| `eval`         | `true`  | Execute the cell during build                    |
| `echo`         | `true`  | Show the source code in output                   |
| `output`       | `true`  | Show the cell output                             |
| `warning`      | `true`  | Show warnings emitted during execution           |
| `code-fold`    | `false` | Collapse the code block (reader can expand)      |
| `code-summary` | `Code`  | Label on the collapsed code toggle               |
| `fig-cap`      |         | Caption for figure output                        |
| `fig-width`    |         | Width of figure output (inches)                  |
| `fig-height`   |         | Height of figure output (inches)                 |
| `tbl-cap`      |         | Caption for table output                         |
| `classes`      |         | Extra CSS classes on the output container        |

## File location

```
project-root/
├── user_guide/
│   ├── 00-introduction.qmd
│   ├── 01-installation.qmd
│   └── ...
├── great-docs.yml
└── pyproject.toml
```

The numeric prefix (`00-`, `01-`, ...) controls sidebar ordering.
Great Docs strips it for clean URLs.
