# Function to merge similar clusters

This function will iteratively merge clusters that result similar in the
centroid. The merging will stop as soon as the cosine similarity between
all pairs of clusters is below \`cutoff\`

## Usage

``` r
merge_clusters(x, cutoff = 0.8)
```

## Arguments

- x:

  bascule object.

- cutoff:

  Minimum value of similarity to merge two clusters.

## Value

bascule object.
