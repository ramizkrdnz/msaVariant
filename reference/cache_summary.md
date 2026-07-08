# Summarize what's currently cached

Lists all locally-cached per-gene annotation files, with their size and
date of caching. Useful for inspecting disk use after querying many
genes.

## Usage

``` r
cache_summary()
```

## Value

A \`data.frame\` with columns \`gene\`, \`size_kb\`, \`cached_on\`,
sorted by size descending. Returns an empty data.frame if the cache is
empty.
