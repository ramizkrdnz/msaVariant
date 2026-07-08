# Clear the msaVariant download cache

Removes locally-cached annotation files. Use this if you suspect a
cached file is stale or corrupt, or if you want to free disk space.

## Usage

``` r
clear_cache(gene = NULL)
```

## Arguments

- gene:

  If supplied, only that gene's cached file is removed. If \`NULL\`
  (default), the entire cache is removed.

## Value

Invisibly, the number of files deleted.
