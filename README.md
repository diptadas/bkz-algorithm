# BKZ Algorithm

## Run BKZ algorithm inside Docker

- Go to docker [playground](https://labs.play-with-docker.com/), login and create an instance.

- Clone the source repository.

```
$ git clone https://github.com/diptadas/bkz-algorithm.git && cd bkz-algorithm
Cloning into 'bkz-algorithm'...
remote: Enumerating objects: 3, done.
remote: Counting objects: 100% (3/3), done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 21 (delta 0), reused 1 (delta 0), pack-reused 18
Unpacking objects: 100% (21/21), done.
```

- Run the pre-built docker image and mount the source files.

```
$ docker run -it -v $(pwd):/bkz -w /bkz diptadas/bkz-utils
```

- Compile and run the source.

```
/bkz # g++ main.cpp -lntl && ./a.out
1. BKZ Pruned Enumeration
2. RBKZ
3. RBKZ DP
4. BKZ Randomize Input
Option: 4
Dim(n): 59
Blocksize(k): 30
Seed: 1650


Base: 2056.260198
Loop: 1 Norm: 1817.233062
Loop: 2 Norm: 1817.233062


Shortest norm after BKZ: 1817.233062
Shortest vector:
[-195 -85 -396 -269 -37 178 342 146 -195 357 -68 -313 270 -113 75 -42 22 -449 -132 341 -109 675 329 49 -127 -53 -100 -89 69 -167 -182 180 -289 -60 -252 209 -159 53 -485 7 204 -345 104 81 -276 466 142 35 -82 -69 -194 -452 -160 57 42 -128 -212 57 384]
```

## [Lattice challenge](https://www.latticechallenge.org/svp-challenge/halloffame.php) submissions

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
