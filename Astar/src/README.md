# $A\ ^* \ Path\ Finding$

## 1. Simulation Result 

![sim_result](https://github.com/Lkaho/Motion_Planning_For_Mobile_Robot/blob/main/Astar/src/sim_result.png)

## 2. Benchmark

5 random maps for each heuristic, target is the same, which locates at the right-down corner with z

equal to 0. In the following Table: 

- h1 is Manhattan; 
- h2 is Euclidean; 
- h3 is Diagonal Heuristic; 

TB is Tie Breaker. h1 might overestimate the path cost, it is not admissible but greedier. It shows less visited nodes. h2 and h3 are admissible and h3 is closer to the real path cost, which shows better

performance compared to h2.

|                      | MAP 1 | MAP 2 | MAP 3 | MAP 4 | MAP 5 | AVG   |
| -------------------- | :---: | :---: | :---: | ----- | ----- | ----- |
| time / ms h1         | 0.136 | 0.363 | 0.233 | 0.221 | 0.534 | 0.297 |
| Visited node h1      |  29   |  27   |  28   | 26    | 40    | 30    |
| time / ms h2         | 1.499 | 0.392 | 3.686 | 5.897 | 4.658 | 3.226 |
| Visited node h2      |  998  |  237  |  771  | 505   | 2218  | 946   |
| time / ms h3         | 0.864 | 0.600 | 0.341 | 0.497 | 0.642 | 0.589 |
| Visited node h3      |  70   |  93   |  33   | 130   | 30    | 71    |
| time / ms h3 + TB    | 0.376 | 0.108 | 0.101 | 0.287 | 0.083 | 0.191 |
| Visited node h3 + TB |  27   |  26   |  40   | 27    | 27    | 29    |

