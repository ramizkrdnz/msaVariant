# Show the location of the msaVariant cache

Show the location of the msaVariant cache

## Usage

``` r
cache_location()
```

## Value

The cache directory path (character).

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
cache_location()
#> Cache directory does not exist yet (no annotations have been downloaded).
#> [1] "/tmp/RtmpUwZ6Qa/msaVariant_cache_1a0c39a0e068"
```
