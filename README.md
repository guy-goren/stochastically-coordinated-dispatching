## Stochastic Coordination in Heterogeneous Load Balancing Systems

This repository contains the implementation of our stochastic coordination dispatching policies (SCD) from our paper in PODC 2021 (Full version: https://arxiv.org/abs/2105.09389, YouTube: https://www.youtube.com/watch?v=JdbLsF0nEPY).

### Abstract

Current-day data centers and high-volume cloud services employ a broad set of heterogeneous servers. In such settings, client requests typically arrive at multiple entry points, and dispatching them to servers is an urgent distributed systems problem. This paper presents an efficient solution to the load balancing problem in such systems that improves on and overcomes problems of previous solutions. The load balancing problem is formulated as a stochastic optimization problem, and an efficient algorithmic solution is obtained based on a subtle mathematical analysis of the problem. Finally, extensive evaluation of the solution on simulated data shows that it outperforms previous solutions. Moreover, the resulting dispatching policy can be computed very efficiently, making the solution practically viable.

### SCD

The dispatcher class that implements SCD is "FastHeterogenousDispatcherSplittableFullServerState". This is the fast O(nlogn) implementation. The slower O(n^2) implementation is "SlowHeterogenousDispatcherSplittableFullServerState".

The dispatcher class that implements "Tidal Water Filling" (TWF) from our DISC'20 paper (Conference version: https://drops.dagstuhl.de/opus/volltexte/2020/13092/, Full version: https://arxiv.org/abs/2008.00793, YouTube: https://www.youtube.com/watch?v=lURNVYFUvJ4) that achieves stochastic coordination in homogeneous systems is "DispatcherTWFSplittableFullServerState".

