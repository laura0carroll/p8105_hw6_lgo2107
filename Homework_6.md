Homework 6
================
Laura O’Carroll
2023-11-28

## Problem 2

First we will download the data.

``` r
weather_df = 
  rnoaa::meteo_pull_monitors(
    c("USW00094728"),
    var = c("PRCP", "TMIN", "TMAX"), 
    date_min = "2022-01-01",
    date_max = "2022-12-31") |>
  mutate(
    name = recode(id, USW00094728 = "CentralPark_NY"),
    tmin = tmin / 10,
    tmax = tmax / 10) |>
  select(name, id, everything())
```

    ## using cached file: /Users/lauraocarroll/Library/Caches/org.R-project.R/R/rnoaa/noaa_ghcnd/USW00094728.dly

    ## date created (size, mb): 2023-11-28 20:52:59.726532 (8.544)

    ## file min/max dates: 1869-01-01 / 2023-11-30

Let’s draw 5000 bootstrap samples and test the distribution of two
quantities estimated from these data:

- r̂2
- log(β̂1∗β̂2)

First make the function.

``` r
boot_sample = function(df) {
  sample_frac(df, replace = TRUE)
}
```

Then let’s pull 5000 samples.

``` r
boot_straps = 
  tibble(strap_number = 1:5000) |> 
  mutate(
    strap_sample = map(strap_number, \(i) boot_sample(df = weather_df))
  )

boot_straps
```

    ## # A tibble: 5,000 × 2
    ##    strap_number strap_sample      
    ##           <int> <list>            
    ##  1            1 <tibble [365 × 6]>
    ##  2            2 <tibble [365 × 6]>
    ##  3            3 <tibble [365 × 6]>
    ##  4            4 <tibble [365 × 6]>
    ##  5            5 <tibble [365 × 6]>
    ##  6            6 <tibble [365 × 6]>
    ##  7            7 <tibble [365 × 6]>
    ##  8            8 <tibble [365 × 6]>
    ##  9            9 <tibble [365 × 6]>
    ## 10           10 <tibble [365 × 6]>
    ## # ℹ 4,990 more rows

Now let’s get estimates for the two quantities of interest for each of
oursamples.

``` r
boot_results = 
  boot_straps |> 
  mutate(
    models = map(strap_sample, \(df) lm(tmax ~ tmin + prcp, data = df)),
    results = map(models, broom::tidy)
  ) |> 
  select(strap_number, models, results) |> 
  unnest(results) |> 
  select(strap_number, models, term, estimate, std.error) |>  
  mutate(
    r2 = map(models, broom::glance)
  ) |> 
  unnest(r2)

compact_boot_results =
  boot_results |> 
  select(strap_number, models, term, estimate, r.squared) |> 
  pivot_wider(names_from = term,
              values_from = estimate) |> 
  rename(beta_0 = "(Intercept)",
         beta_1 = "tmin",
         beta_2 = "prcp") |> 
  mutate(
    log_b1_b2 = log(beta_1 * beta_2)
  )
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `log_b1_b2 = log(beta_1 * beta_2)`.
    ## Caused by warning in `log()`:
    ## ! NaNs produced

``` r
compact_boot_results
```

    ## # A tibble: 5,000 × 7
    ##    strap_number models r.squared beta_0 beta_1    beta_2 log_b1_b2
    ##           <int> <list>     <dbl>  <dbl>  <dbl>     <dbl>     <dbl>
    ##  1            1 <lm>       0.910   7.96  1.03   0.00452      -5.38
    ##  2            2 <lm>       0.944   7.37  1.06  -0.00707     NaN   
    ##  3            3 <lm>       0.925   7.92  1.04  -0.00701     NaN   
    ##  4            4 <lm>       0.911   8.10  1.00  -0.000395    NaN   
    ##  5            5 <lm>       0.899   8.16  0.990  0.00330      -5.72
    ##  6            6 <lm>       0.909   7.94  1.02  -0.00309     NaN   
    ##  7            7 <lm>       0.918   7.87  0.998  0.00136      -6.60
    ##  8            8 <lm>       0.906   8.42  1.00  -0.000283    NaN   
    ##  9            9 <lm>       0.925   7.28  1.05   0.000331     -7.96
    ## 10           10 <lm>       0.919   8.15  1.01  -0.00201     NaN   
    ## # ℹ 4,990 more rows
