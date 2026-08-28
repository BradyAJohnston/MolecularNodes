# Feature Matrix — great-docs.yml

Quick reference for every toggleable feature and its config key.

## Features enabled by default

| Feature                | Config key             | Default | How to disable       |
| ---------------------- | ---------------------- | ------- | -------------------- |
| API reference          | (always on)            | —       | Cannot be disabled   |
| Dark mode toggle       | `dark_mode_toggle`     | `true`  | Set to `false`       |
| Skill generation       | `skill.enabled`        | `true`  | Set to `false`       |
| Changelog              | `changelog.enabled`    | `true`  | Set to `false`       |
| Source links           | `source.enabled`       | `true`  | Set to `false`       |
| Attribution footer     | `attribution`          | `true`  | Set to `false`       |
| `.well-known` skills   | `skill.well_known`     | `true`  | Set to `false`       |
| Dynamic introspection  | `dynamic`              | `true`  | Set to `false`       |

## Features disabled by default

| Feature                | Config key              | How to enable                       |
| ---------------------- | ----------------------- | ----------------------------------- |
| Navbar gradient        | `navbar_style`          | Set to a preset name                |
| Solid navbar color     | `navbar_color`          | Set to a hex color                  |
| Content glow           | `content_style`         | Set to a preset name                |
| Hero section           | `hero.enabled`          | Set to `true`                       |
| Announcement banner    | `announcement`          | Set content string or dict          |
| CLI documentation      | `cli.enabled`           | Set to `true`                       |
| Page tags              | `tags.enabled`          | Set to `true`                       |
| Page status badges     | `page_status`           | Set to `true`                       |
| Sidebar filter         | `sidebar_filter.enabled`| Set to `true`                       |
| Page timestamps        | `site.show_dates`       | Set to `true`                       |
| Back-to-top button     | `back_to_top`           | Set to `true`                       |
| Keyboard navigation    | `keyboard_nav`          | Set to `true`                       |
| Navigation icons       | `nav_icons`             | Add navbar/sidebar icon mappings    |
| Multi-version docs     | `versions`              | Add version entries list            |
| SEO sitemap            | `seo.sitemap`           | Set to `true`                       |
| Social cards           | `seo.social_cards`      | Set to `true`                       |
| Custom sections        | `sections`              | Add section entries list            |
| Authors sidebar        | `authors`               | Add author entries list             |
| Funding info           | `funding`               | Add funding dict                    |
| Logo                   | `logo`                  | Set path or light/dark dict         |
| Favicon                | `favicon`               | Set path                            |
| Display name           | `display_name`          | Set string                          |
| Analytics              | `include_in_header`     | Add script tags                     |

## Feature dependencies

Some features require other features or external setup:

| Feature              | Requires                                   |
| -------------------- | ------------------------------------------ |
| CLI docs             | Click-based CLI module                     |
| Multi-version docs   | Git tags/branches for each version         |
| Source links         | GitHub repository URL (auto-detected)      |
| Changelog            | GitHub Releases on the repository          |
| Navigation icons     | Lucide icon names                          |
| Page status badges   | `status` key in page frontmatter           |
| Page tags            | `tags` key in page frontmatter             |
| SEO sitemap          | `seo.canonical.base_url` must be set       |
