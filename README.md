# BKZ Algorithm

Optimized implementation of BKZ lattice reduction algorithm using `ntl-5.5.2` and `fplll`.

```
git clone https://github.com/diptadas/bkz-algorithm.git
cd bkz-algorithm
docker run -it -v $(pwd):/bkz -w /bkz diptadas/bkz-utils
```

```
g++ main.cpp -lntl
./a.out
```

[Lattice challenge](https://www.latticechallenge.org/svp-challenge/halloffame.php) submissions:

| Dimension | Seed | Norm |
|:---------:|:----:|:----:|
|     94    |  650 | 2473 |
|     91    | 1203 | 2424 |
|     80    |  426 | 2163 |
|     69    | 2675 | 1896 |
|     69    |  555 | 1972 |
|     64    | 2793 | 1880 |
|     59    | 2335 | 1801 |
|     59    | 1650 | 1817 |
