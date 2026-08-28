# Theming Options — Great Docs

## Navbar gradient presets

Set with `navbar_style` in `great-docs.yml`.

| Preset   | Description                              |
| -------- | ---------------------------------------- |
| `sky`    | Blue gradient — professional, calm       |
| `peach`  | Warm peach/coral — friendly, inviting    |
| `prism`  | Multi-color rainbow — vibrant, playful   |
| `lilac`  | Purple/lavender — elegant, creative      |
| `slate`  | Gray/steel — minimal, serious            |
| `honey`  | Gold/amber — warm, earthy                |
| `dusk`   | Deep blue/purple — dramatic, modern      |
| `mint`   | Green/teal — fresh, clean                |

## Content area glow presets

Set with `content_style` in `great-docs.yml`. Uses the same palette
names as navbar presets.

```yaml
# All pages
content_style: lilac

# Homepage only
content_style:
  preset: lilac
  pages: homepage
```

## Solid navbar color

Use `navbar_color` instead of `navbar_style` for a flat color:

```yaml
# Same color for both themes
navbar_color: "#1a1a2e"

# Different colors per theme
navbar_color:
  light: "#ffffff"
  dark: "#1a1a2e"
```

`navbar_style` and `navbar_color` are mutually exclusive. If both
are set, `navbar_style` wins.

## Dark mode

```yaml
dark_mode_toggle: true   # show toggle switch in navbar (default)
dark_mode_toggle: false  # hide toggle, use system preference only
```

All gradient presets and content glows have dark-mode variants that
activate automatically.

## Custom CSS

For fine-grained control, add a custom SCSS/CSS file:

```yaml
include_in_header:
  - text: |
      <link rel="stylesheet" href="custom.css">
```

Place `custom.css` in the project root or `assets/` directory.
Override Great Docs CSS variables for targeted changes:

```css
/* custom.css */
:root {
  --gd-navbar-bg: #2d3748;
  --gd-navbar-text: #e2e8f0;
}
```

## Combining options

Typical aesthetic combinations:

| Style         | Navbar     | Content    | Hero       |
| ------------- | ---------- | ---------- | ---------- |
| Professional  | `sky`      | `sky`      | no starfield|
| Playful       | `prism`    | `peach`    | starfield   |
| Minimal       | `slate`    | none       | no starfield|
| Elegant       | `dusk`     | `lilac`    | starfield   |
| Nature        | `mint`     | `mint`     | no starfield|
